#include "max_internal.h"

struct timeval tm1, tm2 ;

void (*slice_gen)(double*, int, double*, struct detector*, struct model*) ;
void (*slice_merge)(double*, int, double*, double*, double*, long, struct detector*) ;

double maximize(struct max_data *common_data) {
	double avg_likelihood ;
	struct iterate *iter = common_data->iter ;
	if (iter == NULL) {
		fprintf(stderr, "No iterate in max_data!\n") ;
		return -1. ;
	}
	
	struct quaternion *quat = iter->quat ;
	struct params *param = iter->par ;

	gettimeofday(&tm1, NULL) ;
	common_data->within_openmp = 0 ;
	free_max_data(common_data) ;
	
	allocate_memory(common_data) ;
	calculate_rescale(common_data) ;

	#pragma omp parallel default(shared)
	{
		int r ;
		struct max_data *priv_data = malloc(sizeof(struct max_data)) ;
		
		priv_data->within_openmp = 1 ;
		priv_data->iter = iter ;
		allocate_memory(priv_data) ;
		
		#pragma omp for schedule(static,1)
		for (r = 0 ; r < quat->num_rot_p ; ++r)
			calculate_prob(r, priv_data, common_data) ;
		
		normalize_prob(priv_data, common_data) ;
		
		#pragma omp for schedule(static,1)
		for (r = 0 ; r < quat->num_rot_p ; ++r) {
			update_tomogram(r, priv_data, common_data) ;
			merge_tomogram(r, priv_data) ;
		}
		
		combine_information_omp(priv_data, common_data) ;
		
		free_max_data(priv_data) ;
		free(priv_data) ;
	}

	avg_likelihood = combine_information_mpi(common_data) ;
	if (param->need_scaling && param->update_scale)
		update_scale(common_data) ;
	
	return avg_likelihood ;
}

void allocate_memory(struct max_data *data) {
	int detn, d ;
	struct iterate *iter = data->iter ;
	struct model *mod = iter->mod ;
	struct detector *det = iter->det ;
	struct quaternion *quat = iter->quat ;
	struct params *param = iter->par ;
	
	// Both private and common
	data->rmax = calloc(iter->tot_num_data, sizeof(int)) ;
	data->info = calloc(iter->tot_num_data, sizeof(double)) ;
	data->likelihood = calloc(iter->tot_num_data, sizeof(double)) ;
	data->max_exp_p = malloc(iter->tot_num_data * sizeof(double)) ;
	for (d = 0 ; d < iter->tot_num_data ; ++d)
		data->max_exp_p[d] = MAX_EXP_START ;
	if (mod->num_modes > 1)
		data->quat_norm = calloc(mod->num_modes * iter->tot_num_data, sizeof(double)) ;
	
	data->prob = malloc(iter->tot_num_data * sizeof(double*)) ;
	data->place_prob = malloc(iter->tot_num_data * sizeof(int*)) ;
	for (d = 0 ; d < iter->tot_num_data ; ++d) {
		data->prob[d] = NULL ;
		data->place_prob[d] = NULL ;
	}
	data->num_prob = calloc(iter->tot_num_data, sizeof(int)) ;
	if (param->need_scaling && param->update_scale)
		data->psum_d = calloc(iter->tot_num_data, sizeof(double)) ;
		
	if (!data->within_openmp) { // common_data
		data->u = malloc(iter->num_det * sizeof(double*)) ;
		for (detn = 0 ; detn < iter->num_det ; ++detn)
			data->u[detn] = calloc(quat->num_rot_p, sizeof(double)) ;
		data->max_exp = calloc(iter->tot_num_data, sizeof(double)) ;
		data->p_norm = calloc(iter->tot_num_data, sizeof(double)) ;
		data->offset_prob = calloc(iter->tot_num_data * omp_get_max_threads(), sizeof(int)) ;
		
		memset(mod->model2, 0, mod->num_modes*mod->vol*sizeof(double)) ;
		memset(mod->inter_weight, 0, mod->num_modes*mod->vol*sizeof(double)) ;
		print_max_time("alloc", "", param->verbosity > 1 && param->rank == 0) ;
	}
	else { // priv_data
		data->all_views = malloc(iter->num_det * sizeof(double*)) ;
		for (d = 0 ; d < iter->num_det ; ++d)
			data->all_views[d] = malloc(det[d].num_pix * sizeof(double)) ;
		
		data->model = calloc(mod->num_modes*mod->vol, sizeof(double)) ;
		data->weight = calloc(mod->num_modes*mod->vol, sizeof(double)) ;
		
		data->psum_r = calloc(iter->num_det, sizeof(double)) ;
		
		for (d = 0 ; d < iter->tot_num_data ; ++d) {
			data->prob[d] = malloc(4 * sizeof(double)) ;
			data->place_prob[d] = malloc(4 * sizeof(int)) ;
		}
		
		// Only for background-aware update
		if (det[0].with_bg && param->need_scaling) {
			data->mask = malloc(iter->num_det * sizeof(uint8_t*)) ;
			data->G_old = malloc(iter->num_det * sizeof(double*)) ;
			data->G_new = malloc(iter->num_det * sizeof(double*)) ;
			data->G_mid = malloc(iter->num_det * sizeof(double*)) ;
			data->G_latest = malloc(iter->num_det * sizeof(double*)) ;
			data->W_old = malloc(iter->num_det * sizeof(double*)) ;
			data->W_new = malloc(iter->num_det * sizeof(double*)) ;
			data->W_mid = malloc(iter->num_det * sizeof(double*)) ;
			data->W_latest = malloc(iter->num_det * sizeof(double*)) ;
			for (detn = 0 ; detn < iter->num_det ; ++detn) {
				data->mask[detn] = calloc(det[detn].num_pix, sizeof(uint8_t)) ;
				data->G_old[detn] = calloc(det[detn].num_pix, sizeof(double)) ;
				data->G_new[detn] = calloc(det[detn].num_pix, sizeof(double)) ;
				data->G_mid[detn] = calloc(det[detn].num_pix, sizeof(double)) ;
				data->G_latest[detn] = calloc(det[detn].num_pix, sizeof(double)) ;
				data->W_old[detn] = calloc(det[detn].num_pix, sizeof(double)) ;
				data->W_new[detn] = calloc(det[detn].num_pix, sizeof(double)) ;
				data->W_mid[detn] = calloc(det[detn].num_pix, sizeof(double)) ;
				data->W_latest[detn] = calloc(det[detn].num_pix, sizeof(double)) ;
			}
		}
	}
}

void free_max_data(struct max_data *data) {
	int detn, d ;
	if (data == NULL || data->iter == NULL)
		return ;

	struct iterate *iter = data->iter ;
	struct params *param = iter->par ;
	struct detector *det = iter->det ;

	free(data->max_exp_p) ;
	free(data->info) ;
	free(data->likelihood) ;
	free(data->rmax) ;
	free(data->quat_norm) ;
	if (data->prob != NULL) {
		for (d = 0 ; d < iter->tot_num_data ; ++d) {
			free(data->prob[d]) ;
			if (data->place_prob != NULL)
				free(data->place_prob[d]) ;
		}
	}
	free(data->prob) ;
	free(data->place_prob) ;
	free(data->num_prob) ;
	free(data->psum_d) ;
	
	if (!data->within_openmp) {
		free(data->max_exp) ;
		if (data->u != NULL) {
			for (detn = 0 ; detn < iter->num_det ; ++detn)
				free(data->u[detn]) ;
		}
		free(data->u) ;
		free(data->p_norm) ;
		free(data->offset_prob) ;
	}
	else {
		if (data->all_views != NULL) {
			for (detn = 0 ; detn < iter->num_det ; ++detn)
				free(data->all_views[detn]) ;
		}
		free(data->all_views) ;
		free(data->model) ;
		free(data->weight) ;
		free(data->psum_r) ;
		if (det[0].with_bg && param->need_scaling) {
			for (detn = 0 ; detn < iter->num_det ; ++detn) {
				if (data->mask != NULL) free(data->mask[detn]) ;
				if (data->G_old != NULL) free(data->G_old[detn]) ;
				if (data->G_new != NULL) free(data->G_new[detn]) ;
				if (data->G_mid != NULL) free(data->G_mid[detn]) ;
				if (data->G_latest != NULL) free(data->G_latest[detn]) ;
				if (data->W_old != NULL) free(data->W_old[detn]) ;
				if (data->W_new != NULL) free(data->W_new[detn]) ;
				if (data->W_mid != NULL) free(data->W_mid[detn]) ;
				if (data->W_latest != NULL) free(data->W_latest[detn]) ;
			}
			free(data->mask) ;
			free(data->G_old) ;
			free(data->G_new) ;
			free(data->G_mid) ;
			free(data->G_latest) ;
			free(data->W_old) ;
			free(data->W_new) ;
			free(data->W_mid) ;
			free(data->W_latest) ;
		}
	}

	data->model = NULL ;
	data->weight = NULL ;
	data->all_views = NULL ;
	data->psum_r = NULL ;
	data->mask = NULL ;
	data->G_old = NULL ;
	data->G_new = NULL ;
	data->G_mid = NULL ;
	data->G_latest = NULL ;
	data->W_old = NULL ;
	data->W_new = NULL ;
	data->W_mid = NULL ;
	data->W_latest = NULL ;
	data->max_exp_p = NULL ;
	data->info = NULL ;
	data->likelihood = NULL ;
	data->rmax = NULL ;
	data->quat_norm = NULL ;
	data->prob = NULL ;
	data->place_prob = NULL ;
	data->num_prob = NULL ;
	data->psum_d = NULL ;
	data->max_exp = NULL ;
	data->u = NULL ;
	data->p_norm = NULL ;
	data->offset_prob = NULL ;
}

void print_max_time(char *pre_tag, char *post_tag, int flag) {
	if (!flag)
		return ;
	
	double diff ;
	double time_1 = tm1.tv_sec + tm1.tv_usec*1.e-6 ;
	double time_2 = tm2.tv_sec + tm2.tv_usec*1.e-6 ;
	
	if (time_1 > time_2) {
		gettimeofday(&tm2, NULL) ;
		diff = tm2.tv_sec + tm2.tv_usec*1.e-6 - time_1 ;
	}
	else {
		gettimeofday(&tm1, NULL) ;
		diff = tm1.tv_sec + tm1.tv_usec*1.e-6 - time_2 ;
	}
	
	fprintf(stderr, "\t%s\t%f %s\n", pre_tag, diff, post_tag) ;
}
