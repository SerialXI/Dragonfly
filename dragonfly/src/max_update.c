#include "max_internal.h"

void update_tomogram(int r, struct max_data *priv, struct max_data *common) {
	double scalemax ;
	struct iterate *iter = common->iter ;
	struct params *param = iter->par ;
	struct detector *det = iter->det ;
	
	scalemax = calc_psum_r(r, priv, common) ;
	
	if (det[0].with_bg && param->need_scaling)
		update_tomogram_bg(r, scalemax, priv, common) ;
	else
		update_tomogram_nobg(r, priv, common) ;

}

void merge_tomogram(int r, struct max_data *priv) {
	int detn, mode, rotind ;
	struct iterate *iter = priv->iter ;
	struct model *mod = iter->mod ;
	struct detector *det = iter->det ;
	struct quaternion *quat = iter->quat ;
	struct params *param = iter->par ;
	
	rotind = (r*param->num_proc + param->rank) / param->num_modes ;
	mode = (r*param->num_proc + param->rank) % param->num_modes ;
	if (rotind >= quat->num_rot) {
		mode = r*param->num_proc + param->rank - param->num_modes * (quat->num_rot - 1) ;
		rotind = 0 ;
	}
	
	// If no data frame has any probability for this orientation, don't merge
	for (detn = 0 ; detn < iter->num_det ; ++detn)
	if (priv->psum_r[detn] > 0.)
		(*slice_merge)(&quat->quats[rotind*5], mode, priv->all_views[detn], priv->model, priv->weight, mod->size, &det[detn]) ;
	
	if ((param->verbosity > 2) && ((r*param->num_proc + param->rank)%(quat->num_rot * param->num_modes / 10) == 0))
		fprintf(stderr, "\t\tFinished r = %d/%d\n", r*param->num_proc + param->rank, quat->num_rot * param->num_modes + param->nonrot_modes) ;
}

double calc_psum_r(int r, struct max_data *priv, struct max_data *common) {
	int dset = 0, d, curr_d, detn, rotind, mode, t, ind ;
	struct iterate *iter = common->iter ;
	struct dataset *frames = iter->dset ;
	struct model *mod = iter->mod ;
	struct quaternion *quat = iter->quat ;
	struct params *param = iter->par ;
	double temp, scalemax ;
	struct dataset *curr = frames ;
	double **prob = priv->prob ;
	int **place_prob = priv->place_prob, *num_prob = priv->num_prob ;
	
	scalemax = -DBL_MAX ;
	memset(priv->psum_r, 0, (iter->num_det)*sizeof(double)) ;
	rotind = (r*param->num_proc + param->rank) / param->num_modes ;
	mode = (r*param->num_proc + param->rank) % param->num_modes ;
	if (rotind >= quat->num_rot) {
		mode = r*param->num_proc + param->rank - param->num_modes * (quat->num_rot - 1) ;
		rotind = 0 ;
	}
	
	while (curr != NULL) {
		detn = iter->det_mapping[dset] ;
		for (curr_d = 0 ; curr_d < curr->num_data ; ++curr_d) {
			// Calculate frame number in full list
			d = curr->num_offset + curr_d ;

			// check if frame is blacklisted
			if (iter->blacklist[d])
				continue ;
			
			// For refinement, check if frame should be processed
			if (param->refine && mode < param->num_modes) {
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
			for (t = 0 ; t < num_prob[d] ; ++t)
			if (r*param->num_proc + param->rank == place_prob[d][t]) {
				ind = t ;
				break ;
			}
			if (ind == -1)
				continue ;
			
			// Exponentiate log-likelihood and normalize to get probabilities
			temp = prob[d][ind] ;
			if (curr->ftype < 2)
				prob[d][ind] = exp(iter->beta[d] * (prob[d][ind] - common->max_exp[d])) / common->p_norm[d] ; 
			else
				prob[d][ind] = exp(iter->beta[d] * (prob[d][ind] - common->max_exp[d]) / 2. / param->sigmasq) / common->p_norm[d] ;
			
			//if (param->need_scaling)
			//	priv->likelihood[d] += prob[d][ind] * (temp - iter->sum_fact[d] + frames->count[d]*log(iter->scale[d])) ;
			//else
			priv->likelihood[d] += prob[d][ind] * (temp - iter->sum_fact[d]) ;
			
			// Calculate denominator for update rule
			if (param->need_scaling) {
				priv->psum_r[detn] += prob[d][ind] * iter->scale[d] ;
				
				// Calculate denominator for scale factor update rule
				if (param->update_scale)
					priv->psum_d[d] += prob[d][ind] * common->u[detn][r] * iter->rescale[detn] ;
			}
			else
				priv->psum_r[detn] += prob[d][ind] ; 
			
			// Skip if probability is very low (saves time)
			if (prob[d][ind] < PROB_MIN)
				continue ;
			
			// If multiple modes, calculate occupancy of frame into each mode
			if (mod->num_modes > 1)
				priv->quat_norm[d*mod->num_modes + mode] += prob[d][ind] ;
			
			// Calculate mutual information of probability distribution
			if (mode >= param->num_modes)
				priv->info[d] += prob[d][ind] * log(prob[d][ind] / mod->num_modes) ;
			else
				priv->info[d] += prob[d][ind] * log(prob[d][ind] / quat->quats[rotind*5 + 4] * mod->num_modes) ;
			
			if (param->need_scaling && iter->scale[d] > scalemax)
				scalemax = iter->scale[d] ;
		}
		dset++ ;
		curr = curr->next ;
	}
	
	return scalemax ;
}

void update_tomogram_nobg(int r, struct max_data *priv, struct max_data *common) {
	int dset = 0, t, d, curr_d, pixel, detn, ind, rotind ;
	double *view ;
	struct iterate *iter = common->iter ;
	struct dataset *frames = iter->dset ;
	struct detector *det = iter->det ;
	struct params *param = iter->par ;
	struct dataset *curr = frames ;
	double **prob = priv->prob ;
	int **place_prob = priv->place_prob ;
	
	for (detn = 0 ; detn < iter->num_det ; ++detn) 
		memset(priv->all_views[detn], 0, det[detn].num_pix*sizeof(double)) ;
	rotind = (r*param->num_proc + param->rank) / param->num_modes ;
	
	while (curr != NULL) {
		// Calculate slice for current detector
		detn = iter->det_mapping[dset] ;
		view = priv->all_views[detn] ;
		
		for (curr_d = 0 ; curr_d < curr->num_data ; ++curr_d) {
			// Calculate frame number in full list
			d = curr->num_offset + curr_d ;
			
			// check if frame is blacklisted
			if (iter->blacklist[d])
				continue ;
			
			// For refinement, check if frame should be processed
			if (param->refine && rotind < iter->quat->num_rot) {
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
			for (t = 0 ; t < priv->num_prob[d] ; ++t)
			if (r*param->num_proc + param->rank == place_prob[d][t]) {
				ind = t ;
				break ;
			}
			if (ind == -1)
				continue ;
			
			// Skip if probability is very low (saves time)
			if (!(prob[d][ind] > PROB_MIN))
				continue ;
			
			if (curr->ftype == 0) {
				// For all pixels with one photon
				for (t = 0 ; t < curr->ones[curr_d] ; ++t) {
					pixel = curr->place_ones[curr->ones_accum[curr_d] + t] ;
					if (det[detn].raw_mask[pixel] < 2)
						view[pixel] += prob[d][ind] ;
				}
				
				// For all pixels with count_multi photons
				for (t = 0 ; t < curr->multi[curr_d] ; ++t) {
					pixel = curr->place_multi[curr->multi_accum[curr_d] + t] ;
					if (det[detn].raw_mask[pixel] < 2)
						view[pixel] += curr->count_multi[curr->multi_accum[curr_d] + t] * prob[d][ind] ;
				}
			}
			else if (curr->ftype == 1) {
				for (t = 0 ; t < curr->num_pix ; ++t)
					view[t] += curr->int_frames[curr_d*curr->num_pix + t] * prob[d][ind] ;
			}
			else if (curr->ftype == 2) { // Gaussian EMC update without scaling
				for (t = 0 ; t < curr->num_pix ; ++t)
				if (det[detn].raw_mask[t] < 2)
					view[t] += curr->frames[curr_d*curr->num_pix + t] * prob[d][ind] ;
			}
		}
		
		curr = curr->next ;
		dset++ ;
	}
	
	for (detn = 0 ; detn < iter->num_det ; ++detn) 
	for (t = 0 ; t < det[detn].num_pix ; ++t)
	if (priv->psum_r[detn] > 0.)
		priv->all_views[detn][t] /= priv->psum_r[detn] ;
}

void gradient_rt(int r, struct max_data *priv, double **views, double **gradients) {
	int dset = 0, t, d, curr_d, pixel, detn, ind, rotind ;
	double val, *grad, *view ;
	struct iterate *iter = priv->iter ;
	struct dataset *frames = iter->dset ;
	struct detector *det = iter->det ;
	struct params *param = iter->par ;
	struct dataset *curr = frames ;
	int **place_prob = priv->place_prob ;
	double **prob = priv->prob ;
	
	// Initialization:
	// 	grad = -sum_d P_dr * phi_d
	// 	mask = 1
	for (detn = 0 ; detn < iter->num_det ; ++detn) {
		grad = gradients[detn] ;
		view = views[detn] ;
		
		for (t = 0 ; t < det[detn].num_pix ; ++t) {
			if (priv->mask[detn][t] < 128) {
				grad[t] = - priv->psum_r[detn] ;
				priv->mask[detn][t] = 1 ;
			}
			else if (priv->mask[detn][t] == 160) {
				grad[t] = -1e100 ;
			}
			else {
				grad[t] = 0. ;
			}
		}
	}
	rotind = (r*param->num_proc + param->rank) / param->num_modes ;
	
	while (curr != NULL) {
		detn = iter->det_mapping[dset] ;
		grad = gradients[detn] ;
		view = views[detn] ;
		
		for (curr_d = 0 ; curr_d < curr->num_data ; ++curr_d) {
			// Calculate frame number in full list
			d = curr->num_offset + curr_d ;
			
			// check if frame is blacklisted
			if (iter->blacklist[d])
				continue ;
			
			// For refinement, check if frame should be processed
			if (param->refine && rotind < iter->quat->num_rot) {
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
			for (t = 0 ; t < priv->num_prob[d] ; ++t)
			if (r*param->num_proc + param->rank == place_prob[d][t]) {
				ind = t ;
				break ;
			}
			if (ind == -1)
				continue ;
			
			// Skip if probability is very low (saves time)
			if (!(prob[d][ind] > PROB_MIN))
				continue ;
			
			// Currently only working with type-0 data
			// For each pixel with one photon
			for (t = 0 ; t < curr->ones[curr_d] ; ++t) {
				pixel = curr->place_ones[curr->ones_accum[curr_d] + t] ;
				if (priv->mask[detn][pixel] < 128) {
					val = view[pixel] * iter->scale[d] + iter->bgscale[d] * det[detn].background[pixel] ;
					grad[pixel] += prob[d][ind] * iter->scale[d] / val ;
					priv->mask[detn][pixel] = 0 ;
				}
				else if (priv->mask[detn][pixel] == 128) {
					grad[pixel] += prob[d][ind] ;
				}
				else if (priv->mask[detn][pixel] == 160) {
					grad[pixel] = fmax(grad[pixel], iter->scale[d] / iter->bgscale[d]) ;
				}
			}
			
			// For each pixel with count_multi photons
			for (t = 0 ; t < curr->multi[curr_d] ; ++t) {
				pixel = curr->place_multi[curr->multi_accum[curr_d] + t] ;
				if (priv->mask[detn][pixel] < 128) {
					val = view[pixel] * iter->scale[d] + iter->bgscale[d] * det[detn].background[pixel] ;
					grad[pixel] += prob[d][ind] * curr->count_multi[curr->multi_accum[curr_d] + t] * iter->scale[d] / val ;
					priv->mask[detn][pixel] = 0 ;
				}
				else if (priv->mask[detn][pixel] == 128) {
					grad[pixel] += prob[d][ind] * curr->count_multi[curr->multi_accum[curr_d] + t] ;
				}
				else if (priv->mask[detn][pixel] == 160) {
					grad[pixel] = fmax(grad[pixel], iter->scale[d] / iter->bgscale[d]) ;
				}
			}
		}
		
		dset++ ;
		curr = curr->next ;
	}
}

void update_tomogram_bg(int r, double scalemax, struct max_data *priv, struct max_data *common) {
	int i, t, detn ;
	int nmask, tot_num_pix = 0 ;
	double val ;
	double grad_tol = 1.e-4 ;
	
	struct iterate *iter = common->iter ;
	// struct dataset *frames = iter->dset ; // unused - accessed via gradient_rt
	struct detector *det = iter->det ;
	struct params *param = iter->par ;
	double mean_count = 0. ;
	for (detn = 0 ; detn < iter->num_det ; ++detn)
		mean_count += iter->mean_count[detn] ;
	mean_count /= iter->num_det ;
	
	nmask = 0 ;
	for (detn = 0 ; detn < iter->num_det ; ++detn) {
		memset(priv->all_views[detn], 0, det[detn].num_pix*sizeof(double)) ;
		memset(priv->W_old[detn], 0, det[detn].num_pix*sizeof(double)) ;
		memset(priv->W_new[detn], 0, det[detn].num_pix*sizeof(double)) ;
		memset(priv->W_latest[detn], 0, det[detn].num_pix*sizeof(double)) ;
		tot_num_pix += det[detn].num_pix ;
		
		for (t = 0 ; t < det[detn].num_pix ; ++t) {
			if (det[detn].raw_mask[t] < 2) {
				if (det[detn].background[t] / scalemax < 1.e-2 * det[detn].powder[t]) {
					priv->mask[detn][t] = 128 ;
				}
				else {
					priv->W_new[detn][t] = 1.e-8 - det[detn].background[t] / scalemax ;
					priv->mask[detn][t] = 0 ;
				}
			}
			else {
				priv->mask[detn][t] = 255 ;
			}
		}
	}
	
	// mask values:
	//   0 : Still to be optimized
	// 128 : Negligible background (can do standard update)
	// 160 : Signal to gradient_rt() to calculate scale_max rather than gradient
	// 192 : Had to adjust search window to have opposite signs
	// 255 : Optimal pixel
	// Set search bounds
	// Calculate G(0) and check sign
	gradient_rt(r, priv, priv->W_old, priv->G_old) ;
	for (detn = 0 ; detn < iter->num_det ; ++detn)
	for (t = 0 ; t < det[detn].num_pix ; ++t) {
		if (fabs(priv->G_old[detn][t]) < grad_tol) {
			priv->all_views[detn][t] = 0. ;
			priv->mask[detn][t] = 255 ;
		}
		else if (priv->G_old[detn][t] > 0.) {
			val = mean_count / det[detn].num_pix ; // Average photons/pixel
			priv->W_new[detn][t] = det[detn].powder[t] > val ? det[detn].powder[t] : val ; // Check against powder sum value
			priv->W_new[detn][t] *= 2. ; // Double
		}
	}
	
	// Find far end of search window
	for (i = 0 ; i < 10 ; ++i) {
		nmask = 0 ;
		gradient_rt(r, priv, priv->W_new, priv->G_new) ;
		
		for (detn = 0 ; detn < iter->num_det ; ++detn)
		for (t = 0 ; t < det[detn].num_pix ; ++t) {
			if (priv->mask[detn][t] == 0) {
				if (priv->G_old[detn][t] < 0. && priv->G_new[detn][t] < 0.) {      // If both negative, calculate scale_max
					priv->mask[detn][t] = 160 ; // Forces gradient_rt to calculate scale_max
				}
				else if (priv->G_old[detn][t] > 0. && priv->G_new[detn][t] > 0.) { // If both positive, double W_new
					priv->W_new[detn][t] *= 2 ;
				}
				else {                                                             // If opposite signs, stop updating
					priv->mask[detn][t] = 192 ;
					nmask++ ;
				}
			}
			else if (priv->mask[detn][t] == 160) { // Set W_new to be just above minimum background
				priv->W_new[detn][t] = 1.e-8 - det[detn].background[t] / priv->G_new[detn][t] ;
				priv->mask[detn][t] = 0 ;
			}
			else {
				nmask++ ;
			}
		}
		
		if (nmask == tot_num_pix)
			break ;
	}
	if (i == 10 && nmask/((double)tot_num_pix) < 0.9)
		fprintf(stderr, "%.5d bad search bounds, %d/%d\n", r*param->num_proc + param->rank, nmask, tot_num_pix) ;
	
	// Bounded root-finding using Ridder's method
	for (i = 0 ; i < 50 ; ++i) {
		// Phase 1: Convergence checks and midpoint computation
		nmask = 0 ;
		for (detn = 0 ; detn < iter->num_det ; ++detn)
		for (t = 0 ; t < det[detn].num_pix ; ++t) {
			if (priv->mask[detn][t] == 255) {
				nmask++ ;
			}
			else if (priv->mask[detn][t] == 1) {
				priv->mask[detn][t] = 255 ;
				priv->all_views[detn][t] = 0. ;
				nmask++ ;
			}
			else if (priv->mask[detn][t] == 128) {
				priv->all_views[detn][t] = priv->G_old[detn][t] / priv->psum_r[detn] ;
				priv->mask[detn][t] = 255 ;
				nmask++ ;
			}
			else if (priv->mask[detn][t] == 192) {
				priv->mask[detn][t] = 0 ;
				priv->W_mid[detn][t] = 0.5 * (priv->W_old[detn][t] + priv->W_new[detn][t]) ;
			}
			else if (priv->mask[detn][t] == 0) {
				priv->W_mid[detn][t] = 0.5 * (priv->W_old[detn][t] + priv->W_new[detn][t]) ;
			}
		}
		if (nmask == tot_num_pix)
			break ;
		
		// Phase 2: Evaluate gradient at midpoint
		gradient_rt(r, priv, priv->W_mid, priv->G_mid) ;
		
		// Phase 3: Handle no-photon pixels, check midpoint convergence, compute Ridder point
		for (detn = 0 ; detn < iter->num_det ; ++detn)
		for (t = 0 ; t < det[detn].num_pix ; ++t) {
			if (priv->mask[detn][t] == 1) {
				priv->mask[detn][t] = 255 ;
				priv->all_views[detn][t] = 0. ;
				nmask++ ;
			}
			else if (priv->mask[detn][t] == 0) {
				if (fabs(priv->G_mid[detn][t]) < grad_tol) {
					priv->all_views[detn][t] = priv->W_mid[detn][t] ;
					priv->mask[detn][t] = 255 ;
					nmask++ ;
				}
				else {
					double s = sqrt(priv->G_mid[detn][t]*priv->G_mid[detn][t] - priv->G_old[detn][t]*priv->G_new[detn][t]) ;
					if (s < 1.e-15) {
						priv->W_latest[detn][t] = priv->W_mid[detn][t] ;
					}
					else {
						priv->W_latest[detn][t] = priv->W_mid[detn][t] + (priv->W_mid[detn][t] - priv->W_old[detn][t]) * copysign(1.0, priv->G_old[detn][t]) * priv->G_mid[detn][t] / s ;
						double lo = fmin(priv->W_old[detn][t], priv->W_new[detn][t]) ;
						double hi = fmax(priv->W_old[detn][t], priv->W_new[detn][t]) ;
						if (priv->W_latest[detn][t] <= lo || priv->W_latest[detn][t] >= hi)
							priv->W_latest[detn][t] = priv->W_mid[detn][t] ;
					}
				}
			}
		}
		if (nmask == tot_num_pix)
			break ;
		
		// Phase 4: Evaluate gradient at Ridder point
		gradient_rt(r, priv, priv->W_latest, priv->G_latest) ;
		
		// Phase 5: Check Ridder convergence and update brackets
		for (detn = 0 ; detn < iter->num_det ; ++detn)
		for (t = 0 ; t < det[detn].num_pix ; ++t) {
			if (priv->mask[detn][t] == 1) {
				priv->mask[detn][t] = 255 ;
				priv->all_views[detn][t] = 0. ;
				nmask++ ;
			}
			else if (priv->mask[detn][t] == 0) {
				if (fabs(priv->G_latest[detn][t]) < grad_tol) {
					priv->all_views[detn][t] = priv->W_latest[detn][t] ;
					priv->mask[detn][t] = 255 ;
					nmask++ ;
				}
				else if (priv->G_mid[detn][t] * priv->G_latest[detn][t] < 0.) {
					// Tightest bracket: between midpoint and Ridder point
					priv->W_old[detn][t] = priv->W_mid[detn][t] ;
					priv->G_old[detn][t] = priv->G_mid[detn][t] ;
					priv->W_new[detn][t] = priv->W_latest[detn][t] ;
					priv->G_new[detn][t] = priv->G_latest[detn][t] ;
				}
				else if (priv->G_old[detn][t] * priv->G_latest[detn][t] < 0.) {
					priv->W_new[detn][t] = priv->W_latest[detn][t] ;
					priv->G_new[detn][t] = priv->G_latest[detn][t] ;
				}
				else {
					priv->W_old[detn][t] = priv->W_latest[detn][t] ;
					priv->G_old[detn][t] = priv->G_latest[detn][t] ;
				}
			}
		}
		
		if (nmask == tot_num_pix)
			break ;
	}
	if (i == 50 && nmask/((double)tot_num_pix) < 0.9)
		fprintf(stderr, "%.5d not converged, %d/%d\n", r*param->num_proc + param->rank, nmask, tot_num_pix) ;
}
