#include "max_internal.h"

void calculate_rescale(struct max_data *data) {
	int detn ;
	struct iterate *iter = data->iter ;
	struct model *mod = iter->mod ;
	struct detector *det = iter->det ;
	struct quaternion *quat = iter->quat ;
	struct params *param = iter->par ;
	double *total = calloc(iter->num_det, sizeof(double)) ;
	char res_string[1024] = {'\0'}  ;
	
	// Calculate rescale factor by calculating mean model value over detector
	#pragma omp parallel default(shared)
	{
		int r, t, detn, mode, rotind ;
		double *priv_total = calloc(iter->num_det, sizeof(double)) ;
		double **views = malloc(iter->num_det * sizeof(double*)) ;
		for (detn = 0 ; detn < iter->num_det ; ++detn)
			views[detn] = malloc(det[detn].num_pix * sizeof(double)) ;
		
		#pragma omp for schedule(static,1)
		for (r = 0 ; r < quat->num_rot_p ; ++r) {
			rotind = (r*param->num_proc + param->rank) / param->num_modes ;
			mode = (r*param->num_proc + param->rank) % param->num_modes ;
			if (rotind >= quat->num_rot) {
				mode = r*param->num_proc + param->rank - param->num_modes * (quat->num_rot - 1) ;
				rotind = 0 ;
			}
			//fprintf(stderr, "%d: %.3d - %.2d %.2d\n", omp_get_thread_num(), r, rotind, mode) ;
			
			for (detn = 0 ; detn < iter->num_det ; ++detn) {
				// Second argument being 0. tells slice_gen to generate un-rescaled tomograms
				(*slice_gen)(&quat->quats[rotind*5], mode, views[detn], &det[detn], mod) ;
				
				for (t = 0 ; t < det[detn].num_pix ; ++t)
				if (det[detn].raw_mask[t] < 1)
					data->u[detn][r] += views[detn][t] ;
				
				priv_total[detn] += quat->quats[rotind*5 + 4] * data->u[detn][r] ;
			}
		}
		
		#pragma omp critical(total)
		{
			for (detn = 0 ; detn < iter->num_det ; ++detn)
				total[detn] += priv_total[detn] ;
		}
		for (detn = 0 ; detn < iter->num_det ; ++detn)
			free(views[detn]) ;
		free(views) ;
		free(priv_total) ;
	}
	
	MPI_Allreduce(MPI_IN_PLACE, total, iter->num_det, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
	
	sprintf(res_string, "(=") ;
	for (detn = 0 ; detn < iter->num_det ; ++detn) {
		iter->rescale[detn] = iter->mean_count[detn] / total[detn] * mod->num_modes ;
		sprintf(res_string + strlen(res_string), " %.6e", iter->rescale[detn]) ;
	}
	sprintf(res_string + strlen(res_string), ")") ;
	free(total) ;
	print_max_time("rescale", res_string, param->verbosity > 1 && param->rank == 0) ;
}

void calculate_prob(int r, struct max_data *priv, struct max_data *common) {
	int dset = 0, t, d, curr_d, pixel, mode, rotind, detn, old_detn = -1, ind ;
	struct iterate *iter = common->iter ;
	struct dataset *frames = iter->dset ;
	struct model *mod = iter->mod ;
	struct detector *det = iter->det ;
	struct quaternion *quat = iter->quat ;
	struct params *param = iter->par ;
	
	struct dataset *curr = frames ;
	double ipred, pval, *view ;
	int *num_prob = priv->num_prob, **place_prob = priv->place_prob ;
	double **prob = priv->prob ;
	
	rotind = (r*param->num_proc + param->rank) / param->num_modes ;
	mode = (r*param->num_proc + param->rank) % param->num_modes ;
	if (rotind >= quat->num_rot) {
		mode = r*param->num_proc + param->rank - param->num_modes * (quat->num_rot - 1) ;
		rotind = 0 ;
	}
	
	// Linked list of data sets from different files
	while (curr != NULL) {
		// Calculate slice for current detector
		detn = iter->det_mapping[dset] ;
		view = priv->all_views[detn] ;
		if (detn != old_detn) {
			(*slice_gen)(&quat->quats[rotind*5], mode, view, &det[detn], mod) ;
			if (!(det[0].with_bg && param->need_scaling)) {
				for (t = 0 ; t < det[detn].num_pix ; ++t)
					view[t] = log(view[t]) ;
			}
		}
		old_detn = detn ;
		
		// For each frame in data set
		for (curr_d = 0 ; curr_d < curr->num_data ; ++curr_d) {
			// Calculate frame number in full list
			d = curr->num_offset + curr_d ;
			
			// check if frame is blacklisted
			if (iter->blacklist[d])
				continue ;
			
			// For refinement, check if frame should be processed
			if (param->refine) {
				ind = -1 ;
				for (t = 0 ; t < iter->num_rel_quat[d] ; ++t)
				if (iter->quat_mapping[r*param->num_proc+param->rank] == iter->rel_quat[d][t]) {
					ind = t ;
					break ;
				}
				if (ind == -1)
					continue ;
			}
			
			if (curr->ftype < 2) {
				// need_scaling is for if we want to assume variable incident intensity
				if (param->need_scaling && (param->iteration > 1 || param->known_scale))
					pval = log(quat->quats[rotind*5 + 4]) - common->u[detn][r] * iter->scale[d] ;
				else
					pval = log(quat->quats[rotind*5 + 4]) - common->u[detn][r] * iter->rescale[detn] ;
			}
			else {
				pval = 0. ;
			}
			
			if (curr->ftype == SPARSE) {
				if (det[0].with_bg && param->need_scaling) {
					// For each pixel with one photon
					for (t = 0 ; t < curr->ones[curr_d] ; ++t) {
						pixel = curr->place_ones[curr->ones_accum[curr_d] + t] ;
						if (det[detn].raw_mask[pixel] < 1) {
							ipred = view[pixel] * iter->scale[d] + iter->bgscale[d] * det[detn].background[pixel] ;
							if (ipred < 0)
								ipred = DBL_MIN ;
							pval += log(ipred) ;
						}
					}
					
					// For each pixel with count_multi photons
					for (t = 0 ; t < curr->multi[curr_d] ; ++t) {
						pixel = curr->place_multi[curr->multi_accum[curr_d] + t] ;
						if (det[detn].raw_mask[pixel] < 1) {
							ipred = view[pixel] * iter->scale[d] + iter->bgscale[d] * det[detn].background[pixel] ;
							if (ipred < 0)
								ipred = DBL_MIN ;
							pval += curr->count_multi[curr->multi_accum[curr_d] + t] * log(ipred) ;
						}
					}
				}
				else {
					// For each pixel with one photon
					for (t = 0 ; t < curr->ones[curr_d] ; ++t) {
						pixel = curr->place_ones[curr->ones_accum[curr_d] + t] ;
						if (det[detn].raw_mask[pixel] < 1)
							pval += view[pixel] ;
					}
					
					// For each pixel with count_multi photons
					for (t = 0 ; t < curr->multi[curr_d] ; ++t) {
						pixel = curr->place_multi[curr->multi_accum[curr_d] + t] ;
						if (det[detn].raw_mask[pixel] < 1)
							pval += curr->count_multi[curr->multi_accum[curr_d] + t] * view[pixel] ;
					}
				}
			}
			else if (curr->ftype == DENSE_INT) {
				for (t = 0 ; t < det[detn].num_pix ; ++t)
				if (det[detn].raw_mask[t] < 1)
					pval += curr->int_frames[curr_d*curr->num_pix + t] * view[t] ;
			}
			else if (curr->ftype == DENSE_DOUBLE) { // Gaussian EMC for double precision data without scaling
				for (t = 0 ; t < det[detn].num_pix ; ++t)
				if (det[detn].raw_mask[t] < 1)
					pval -= pow(curr->frames[curr_d*curr->num_pix + t] - view[t]*iter->rescale[detn], 2.) ;
			}
			
			// Only save value in prob array if it is significant
			if (pval + PDIFF_THRESH/iter->beta[d] > priv->max_exp_p[d]) {
				prob[d][num_prob[d]] = pval ;
				place_prob[d][num_prob[d]] = r*param->num_proc + param->rank ;
				num_prob[d]++ ;
				
				// If num_prob is a power of two, expand array
				if (num_prob[d] >= 4 && (num_prob[d] & (num_prob[d] - 1)) == 0) {
					prob[d] = realloc(prob[d], num_prob[d] * 2 * sizeof(double)) ;
					place_prob[d] = realloc(place_prob[d], num_prob[d] * 2 * sizeof(int)) ;
				}
			}
			
			// Note maximum log-likelihood for each frame among 'r's tested by this MPI rank and OMP rank
			// Recalculate sparse array with new maximum
			if (pval > priv->max_exp_p[d]) {
				priv->max_exp_p[d] = pval ;
				priv->rmax[d] = r*param->num_proc + param->rank ;
				num_prob[d] = resparsify(prob[d], place_prob[d], num_prob[d], pval - PDIFF_THRESH/iter->beta[d]) ;
			}
		}
		
		curr = curr->next ;
		dset++ ;
	}
	
	if ((param->verbosity > 2) && ((r*param->num_proc + param->rank)%(quat->num_rot * param->num_modes / 10) == 0))
		fprintf(stderr, "\t\tFinished r = %d/%d\n", r*param->num_proc + param->rank, quat->num_rot * param->num_modes + param->nonrot_modes) ;
	print_max_time("prob", "", param->verbosity > 1 && (r == quat->num_rot_p-1) && (param->rank == 0)) ;
}

void normalize_prob(struct max_data *priv, struct max_data *common) {
	int r, d, omp_rank = omp_get_thread_num() ;
	struct iterate *iter = common->iter ;
	struct dataset *frames = iter->dset ;
	struct params *param = iter->par ;
	double *priv_norm = calloc(iter->tot_num_data, sizeof(double)) ;
	
	// Calculate max_log_prob over all OpenMP ranks (and the r for that maximum)
	#pragma omp critical(maxexp)
	{
		for (d = 0 ; d < iter->tot_num_data ; ++d)
		if (priv->max_exp_p[d] > common->max_exp_p[d]) {
			common->max_exp_p[d] = priv->max_exp_p[d] ;
			common->rmax[d] = priv->rmax[d] ;
		}
	}
	#pragma omp barrier
	
	if (omp_rank == 0) {
		MPI_Allreduce(common->max_exp_p, common->max_exp, iter->tot_num_data, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD) ;
		
		// Determine 'r' for which log-likelihood is maximum
		for (d = 0 ; d < iter->tot_num_data ; ++d)
		if (common->max_exp[d] != common->max_exp_p[d] || common->max_exp_p[d] == MAX_EXP_START)
			common->rmax[d] = -1 ;
		
		MPI_Allreduce(MPI_IN_PLACE, common->rmax, iter->tot_num_data, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;
	}
	#pragma omp barrier
	
	for (d = 0 ; d < iter->tot_num_data ; ++d)
	for (r = 0 ; r < priv->num_prob[d] ; ++r) 
	if (frames->ftype < 2)
		priv_norm[d] += exp(iter->beta[d] * (priv->prob[d][r] - common->max_exp[d])) ;
	else
		priv_norm[d] += exp(iter->beta[d] * (priv->prob[d][r] - common->max_exp[d]) / 2. / param->sigmasq) ;
	
	#pragma omp critical(psum)
	{
		for (d = 0 ; d < iter->tot_num_data ; ++d)
			common->p_norm[d] += priv_norm[d] ;
	}
	#pragma omp barrier
	
	if (omp_rank == 0)
		MPI_Allreduce(MPI_IN_PLACE, common->p_norm, iter->tot_num_data, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
	#pragma omp barrier
	
	free(priv_norm) ;
	print_max_time("norm", "", param->verbosity > 1 && param->rank == 0 && omp_rank == 0) ;
}

void combine_information_omp(struct max_data *priv, struct max_data *common) {
	int d, r, omp_rank = omp_get_thread_num() ;
	int nthreads = omp_get_num_threads() ;
	struct iterate *iter = common->iter ;
	struct model *mod = iter->mod ;
	struct detector *det = iter->det ;
	struct params *param = iter->par ;
	long x ;
	
	print_max_time("update", "", param->verbosity > 1 && param->rank == 0 && omp_rank == 0) ;
	 
	#pragma omp critical(model)
	{
		for (x = 0 ; x < mod->num_modes * mod->vol ; ++x) {
			mod->model2[x] += priv->model[x] ;
			mod->inter_weight[x] += priv->weight[x] ;
		}
	}
	
	#pragma omp critical(like_info)
	{
		for (d = 0 ; d < iter->tot_num_data ; ++d) {
			common->likelihood[d] += priv->likelihood[d] ;
			common->info[d] += priv->info[d] ;
		}
	}
	
	// Calculate common probabilities if needed for update_scale_bg or to save
	if ((param->need_scaling && param->update_scale && det[0].with_bg) || param->save_prob) {
		// Calculate offsets to combine sparse probabilities for each OpenMP rank
		#pragma omp critical(offset_prob)
		{
			for (d = 0 ; d < iter->tot_num_data ; ++d) {
				common->num_prob[d] += priv->num_prob[d] ;
				for (r = omp_rank + 1 ; r < nthreads ; ++r)
					common->offset_prob[d*nthreads + r] += priv->num_prob[d] ;
			}
		}
		#pragma omp barrier
			
		// Allocate common prob arrays
		#pragma omp for schedule(static,1)
		for (d = 0 ; d < iter->tot_num_data ; ++d) {
			common->prob[d] = malloc(common->num_prob[d] * sizeof(double)) ;
			common->place_prob[d] = malloc(common->num_prob[d] * sizeof(int)) ;
		}
		
		// Populate common->prob array for all d
		for (d = 0 ; d < iter->tot_num_data ; ++d)
		for (r = 0 ; r < priv->num_prob[d] ; ++r) {
			common->prob[d][common->offset_prob[d*nthreads + omp_rank] + r] = priv->prob[d][r] ;
			common->place_prob[d][common->offset_prob[d*nthreads + omp_rank] + r] = priv->place_prob[d][r] ;
		}
		#pragma omp barrier
		
		// Sparsify probs based on PROB_MIN threshold
		#pragma omp for schedule(static,1)
		for (d = 0 ; d < iter->tot_num_data ; ++d)
			common->num_prob[d] = resparsify(common->prob[d], common->place_prob[d], common->num_prob[d], PROB_MIN) ;
	}
	
	if (param->need_scaling && param->update_scale) {
		#pragma omp critical(scale)
		{
			for (d = 0 ; d < iter->tot_num_data ; ++d)
			if (!iter->blacklist[d])
				common->psum_d[d] += priv->psum_d[d] ;
		}
	}
	
	if (mod->num_modes > 1) {
		#pragma omp critical(quat_norm)
		for (d = 0 ; d < iter->tot_num_data * mod->num_modes ; ++d)
			common->quat_norm[d] += priv->quat_norm[d] ;
	}
}

double combine_information_mpi(struct max_data *common) {
	int d ;
	struct iterate *iter = common->iter ;
	struct model *mod = iter->mod ;
	struct detector *det = iter->det ;
	struct params *param = iter->par ;
	double avg_likelihood = 0. ;
	iter->mutual_info = 0. ;
	
	// Combine 3D volumes from all MPI ranks
	if (param->rank) {
		MPI_Reduce(mod->model2, mod->model2, mod->num_modes*mod->vol, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) ;
		MPI_Reduce(mod->inter_weight, mod->inter_weight, mod->num_modes*mod->vol, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) ;
	}
	else {
		MPI_Reduce(MPI_IN_PLACE, mod->model2, mod->num_modes*mod->vol, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) ;
		MPI_Reduce(MPI_IN_PLACE, mod->inter_weight, mod->num_modes*mod->vol, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) ;
	}
	
	// Combine mutual info and likelihood from all MPI ranks
	MPI_Allreduce(MPI_IN_PLACE, common->likelihood, iter->tot_num_data, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
	MPI_Allreduce(MPI_IN_PLACE, common->info, iter->tot_num_data, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
	
	for (d = 0 ; d < iter->tot_num_data ; ++d)
	if (!iter->blacklist[d]) {
		iter->mutual_info += common->info[d] ;
		avg_likelihood += common->likelihood[d] ;
	}
	
	iter->mutual_info /= (iter->tot_num_data - iter->num_blacklist) ;
	avg_likelihood /= (iter->tot_num_data - iter->num_blacklist) ;
	
	// Calculate common probabilities if needed for update_scale_bg or to save
	if ((param->need_scaling && param->update_scale && det[0].with_bg) || param->save_prob) {
		int p, q, tot_num_prob ;
		int *num_prob_p = calloc(iter->tot_num_data * param->num_proc, sizeof(int)) ;
		int *displ_prob_p = NULL ;
		for (d = 0 ; d < iter->tot_num_data ; ++d)
			num_prob_p[d*param->num_proc + param->rank] = common->num_prob[d] ;
		if (param->rank) {
			MPI_Reduce(num_prob_p, NULL, iter->tot_num_data*param->num_proc, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD) ;
			free(num_prob_p) ;
		}
		else {
			MPI_Reduce(MPI_IN_PLACE, num_prob_p, iter->tot_num_data*param->num_proc, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD) ;
			displ_prob_p = malloc(param->num_proc * sizeof(int)) ;
		}
		for (d = 0 ; d < iter->tot_num_data ; ++d) {
			if (param->rank) {
				MPI_Gatherv(common->place_prob[d], common->num_prob[d], MPI_INT, NULL, 0, 0, MPI_INT, 0, MPI_COMM_WORLD) ;
				MPI_Gatherv(common->prob[d], common->num_prob[d], MPI_DOUBLE, NULL, 0, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
			}
			else {
				tot_num_prob = 0 ;
				memset(displ_prob_p, 0, param->num_proc*sizeof(int)) ;
				for (p = 0 ; p < param->num_proc ; ++p) {
					tot_num_prob += num_prob_p[d*param->num_proc + p] ;
					for (q = 0 ; q < p ; ++q)
						displ_prob_p[p] += num_prob_p[d*param->num_proc + q] ;
				}
				common->prob[d] = realloc(common->prob[d], tot_num_prob * sizeof(double)) ;
				common->place_prob[d] = realloc(common->place_prob[d], tot_num_prob * sizeof(int)) ;
				
				MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT, common->place_prob[d], &num_prob_p[d*param->num_proc], displ_prob_p, MPI_INT, 0, MPI_COMM_WORLD) ;
				MPI_Gatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, common->prob[d], &num_prob_p[d*param->num_proc], displ_prob_p, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
				common->num_prob[d] = tot_num_prob ;
			}
		}
		if (!param->rank) {
			free(num_prob_p) ;
			free(displ_prob_p) ;
		}
	}
	
	// Combine scale factor information from all MPI ranks
	if (param->need_scaling && param->update_scale)
		MPI_Allreduce(MPI_IN_PLACE, common->psum_d, iter->tot_num_data, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
	
	if (mod->num_modes > 1)
		MPI_Allreduce(MPI_IN_PLACE, common->quat_norm, iter->tot_num_data*mod->num_modes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
	
	print_max_time("sync", "", param->verbosity > 1 && param->rank == 0) ;
	return avg_likelihood ;
}

int resparsify(double *vals, int *pos, int num_vals, double thresh) {
	int i, j ;
	
	for (i = 0 ; i < num_vals ; ++i) {
		if (vals[i] <= thresh) {
			num_vals-- ;
			for (j = i ; j < num_vals ; ++j) {
				vals[j] = vals[j+1] ;
				pos[j] = pos[j+1] ;
				pos[j+1] = -1 ;
			}
			i-- ;
		}
	}
	
	//vals = realloc(vals, num_vals*sizeof(double)) ;
	//pos = realloc(pos, num_vals*sizeof(int)) ;
	
	return num_vals ;
}
