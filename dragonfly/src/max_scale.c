#include "max_internal.h"

void update_scale(struct max_data *common) {
	int d ;
	struct iterate *iter = common->iter ;
	struct params *param = iter->par ;
	struct detector *det = iter->det ;
	
	if (det[0].with_bg) {
		update_scale_bg(common) ;
	}
	else {
		for (d = 0 ; d < iter->tot_num_data ; ++d)
		if (!iter->blacklist[d])
			iter->scale[d] = iter->fcounts[d] / common->psum_d[d] ;
	}

	print_max_time("scale", "", param->verbosity > 1 && param->rank == 0) ;
}

void gradient_d(struct max_data *common, uint8_t *mask, double *scale, double *grad) {
	int d ;
	struct iterate *iter = common->iter ;
	struct dataset *frames = iter->dset ;
	struct detector *det = iter->det ;
	struct params *param = iter->par ;
	struct quaternion *quat = iter->quat ;
	struct model *mod = iter->mod ;
	
	for (d = 0 ; d < iter->tot_num_data ; ++d)
	if (mask[d] == 0)
		grad[d] = 0. ;
	
	#pragma omp parallel default(shared)
	{
		int r, d, t, detn, curr_d, pixel, rotind, mode, ind, dset ;
		double val ;
		double *view, **views = malloc(iter->num_det * sizeof(double*)) ;
		for (detn = 0 ; detn < iter->num_det ; ++detn)
			views[detn] = malloc(det[detn].num_pix * sizeof(double)) ;
		double *priv_grad = calloc(iter->tot_num_data, sizeof(double)) ;
		struct dataset *curr ;
		
		#pragma omp for schedule(static,1)
	 	for (r = 0 ; r < quat->num_rot_p ; ++r) {
			rotind = (r*param->num_proc + param->rank) / param->num_modes ;
			mode = (r*param->num_proc + param->rank) % param->num_modes ;
			if (rotind >= quat->num_rot) {
				mode = r*param->num_proc + param->rank - param->num_modes * (quat->num_rot - 1) ;
				rotind = 0 ;
			}
			for (detn = 0 ; detn < iter->num_det ; ++detn) {
				(*slice_gen)(&quat->quats[rotind*5], mode, views[detn], &det[detn], mod) ;
				for (t = 0 ; t < det[detn].num_pix ; ++t)
					views[detn][t] *= iter->rescale[detn] ;
			}
			curr = frames ;
			dset = 0 ;
			
			while (curr != NULL) {
				detn = iter->det_mapping[dset] ;
				view = views[detn] ;
				for (curr_d = 0 ; curr_d < curr->num_data ; ++curr_d) {
					d = curr->num_offset + curr_d ;
					
					if (mask[d] > 0)
						continue ;
					
					// For refinement, check if frame should be processed
					if (param->refine) {
						ind = -1 ;
						for (t = 0 ; t < iter->num_rel_quat[d] ; ++t)
						if (iter->quat_mapping[rotind] == iter->rel_quat[d][t]) {
							ind = t ;
							break ;
						}
						if (ind == -1)
							continue ;
					}
					
					// check if current frame has significant probability
					ind = -1 ;
					for (t = 0 ; t < common->num_prob[d] ; ++t)
					if (r*param->num_proc + param->rank == common->place_prob[d][t]) {
						ind = t ;
						break ;
					}
					if (ind == -1)
						continue ;
					
					// For each pixel with one photon
					for (t = 0 ; t < curr->ones[curr_d] ; ++t) {
						pixel = curr->place_ones[curr->ones_accum[curr_d] + t] ;
						if (det[detn].raw_mask[pixel] < 1) {
							val = view[pixel] * scale[d] + iter->bgscale[d] * det[detn].background[pixel] ;
							priv_grad[d] += common->prob[d][ind] * view[pixel] / val ;
						}
					}
					
					// For each pixel with count_multi photons
					for (t = 0 ; t < curr->multi[curr_d] ; ++t) {
						pixel = curr->place_multi[curr->multi_accum[curr_d] + t] ;
						if (det[detn].raw_mask[pixel] < 1) {
							val = view[pixel] * scale[d] + iter->bgscale[d] * det[detn].background[pixel] ;
							priv_grad[d] += common->prob[d][ind] * curr->count_multi[curr->multi_accum[curr_d] + t] * view[pixel] / val ;
						}
					}
				}
				
				dset++ ;
				curr = curr->next ;
			}
		}
		
		#pragma omp critical(grad)
		{
			for (d = 0 ; d < iter->tot_num_data ; ++d)
				grad[d] += priv_grad[d] ;
		}
		
		free(priv_grad) ;
		for (detn = 0 ; detn < iter->num_det ; ++detn)
			free(views[detn]) ;
		free(views) ;
	}
	
	MPI_Allreduce(MPI_IN_PLACE, grad, iter->tot_num_data, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
	
	for (d = 0 ; d < iter->tot_num_data ; ++d)
	if (mask[d] == 0)
		grad[d] -= common->psum_d[d] ;
}

void update_scale_bg(struct max_data *common) {
	int d, i, num_mask ;
	struct iterate *iter = common->iter ;
	struct params *param = iter->par ;
	double *scale_old = calloc(iter->tot_num_data, sizeof(double)) ;
	double *scale_new = calloc(iter->tot_num_data, sizeof(double)) ;
	double *scale_mid = calloc(iter->tot_num_data, sizeof(double)) ;
	double *scale_latest = calloc(iter->tot_num_data, sizeof(double)) ;
	double *Gd_old = calloc(iter->tot_num_data, sizeof(double)) ;
	double *Gd_new = calloc(iter->tot_num_data, sizeof(double)) ;
	double *Gd_mid = calloc(iter->tot_num_data, sizeof(double)) ;
	double *Gd_latest = calloc(iter->tot_num_data, sizeof(double)) ;
	uint8_t *mask = calloc(iter->tot_num_data, sizeof(uint8_t)) ;
	for (d = 0 ; d < iter->tot_num_data ; ++d) {
		if (iter->blacklist[d]) {
			mask[d] = 255 ;
			continue ;
		}
		
		if (param->iteration > 1 || param->known_scale)
			scale_new[d] = 4. * iter->scale[d] ;
		else
			scale_new[d] = 32. ;
	}
	
	// Set search bounds
	// 	Calculate G(0)
	gradient_d(common, mask, scale_old, Gd_old) ;
	for (d = 0 ; d < iter->tot_num_data ; ++d)
	if ((mask[d] != 255) & (Gd_old[d] <= 0.)) {
		// The constrained optimum is at the non-negative boundary.
		iter->scale[d] = 0.001 ;
		mask[d] = 255 ;
	}
	
	// 	Calculate phi_max and G(phi_max)
	for (i = 0 ; i < 5 ; ++i) {
		gradient_d(common, mask, scale_new, Gd_new) ;
		num_mask = 0 ;
		
		for (d = 0 ; d < iter->tot_num_data ; ++d) {
			if (mask[d] != 0) {
				num_mask++ ;
			}
			else if (Gd_new[d] < 0.) {
				mask[d] = 192 ;
				num_mask++ ;
			}
			else {
				scale_new[d] *= 4. ;
			}
		}
		
		if (num_mask > ((double) 0.99*iter->tot_num_data))
			break ;
	}
	if (i == 5)
		fprintf(stderr, "WARNING: Could not find search bounds for %d/%d frames\n", iter->tot_num_data - num_mask, iter->tot_num_data) ;
	
	// Bounded root finding using Ridder's method
	for (i = 0 ; i < 50 ; ++i) {
		// Phase 1: Convergence check and midpoint computation
		num_mask = 0 ;
		for (d = 0 ; d < iter->tot_num_data ; ++d) {
			if (mask[d] == 255) {
				num_mask++ ;
			}
			else if (fabs(scale_old[d] - scale_new[d]) < 1.e-5) {
				mask[d] = 255 ;
				iter->scale[d] = 0.5 * (scale_old[d] + scale_new[d]) ;
				num_mask++ ;
			}
			else if (mask[d] == 192) {
				mask[d] = 0 ;
				scale_mid[d] = 0.5 * (scale_old[d] + scale_new[d]) ;
			}
			else {
				scale_mid[d] = 0.5 * (scale_old[d] + scale_new[d]) ;
			}
		}
		if (num_mask == iter->tot_num_data)
			break ;
		
		// Phase 2: Evaluate gradient at midpoint
		gradient_d(common, mask, scale_mid, Gd_mid) ;
		
		// Phase 3: Check midpoint convergence and compute Ridder point
		for (d = 0 ; d < iter->tot_num_data ; ++d) {
			if (mask[d] == 255) continue ;
			if (fabs(Gd_mid[d]) < 1.e-3) {
				iter->scale[d] = scale_mid[d] ;
				mask[d] = 255 ;
				num_mask++ ;
				continue ;
			}
			double s = sqrt(Gd_mid[d]*Gd_mid[d] - Gd_old[d]*Gd_new[d]) ;
			if (s < 1.e-15) {
				scale_latest[d] = scale_mid[d] ;
			}
			else {
				scale_latest[d] = scale_mid[d] + (scale_mid[d] - scale_old[d]) * copysign(1.0, Gd_old[d]) * Gd_mid[d] / s ;
				double lo = fmin(scale_old[d], scale_new[d]) ;
				double hi = fmax(scale_old[d], scale_new[d]) ;
				if (scale_latest[d] <= lo || scale_latest[d] >= hi)
					scale_latest[d] = scale_mid[d] ;
			}
		}
		if (num_mask == iter->tot_num_data)
			break ;
		
		// Phase 4: Evaluate gradient at Ridder point
		gradient_d(common, mask, scale_latest, Gd_latest) ;
		
		// Phase 5: Check Ridder convergence and update brackets
		for (d = 0 ; d < iter->tot_num_data ; ++d)
		if (mask[d] < 255) {
			if (fabs(Gd_latest[d]) < 1.e-3) {
				iter->scale[d] = scale_latest[d] ;
				mask[d] = 255 ;
				num_mask++ ;
			}
			else if (Gd_mid[d] * Gd_latest[d] < 0.) {
				// Tightest bracket: between midpoint and Ridder point
				scale_old[d] = scale_mid[d] ;
				Gd_old[d] = Gd_mid[d] ;
				scale_new[d] = scale_latest[d] ;
				Gd_new[d] = Gd_latest[d] ;
			}
			else if (Gd_old[d] * Gd_latest[d] < 0.) {
				scale_new[d] = scale_latest[d] ;
				Gd_new[d] = Gd_latest[d] ;
			}
			else {
				scale_old[d] = scale_latest[d] ;
				Gd_old[d] = Gd_latest[d] ;
			}
		}

		if (num_mask == iter->tot_num_data)
			break ;
	}
	if ((i == 50) & (param->rank == 0))
		fprintf(stderr, "WARNING: scale optimization did not converge for %d/%d frames\n", iter->tot_num_data-num_mask, iter->tot_num_data) ;
	
	MPI_Bcast(iter->scale, iter->tot_num_data, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
	
	// Free memory
	free(scale_old) ; free(scale_new) ; free(scale_mid) ; free(scale_latest) ;
	free(Gd_old) ; free(Gd_new) ; free(Gd_mid) ; free(Gd_latest) ;
	free(mask) ;
	char tag[128] ;
	sprintf(tag, "(%d iterations)", i) ;
	print_max_time("scale", tag, param->rank == 0) ;
}
