#include "iterate.h"

void calc_frame_counts(struct iterate *self, uint8_t *blacklist) {
	long t, d, detn, dset_num = 0 ;
	struct dataset *curr = self->dset ;
	struct detector *cdet ;
	int *num_data = calloc(self->num_det, sizeof(int)) ;
	if (blacklist == NULL) {
		fprintf(stderr, "Need to generate an empty blacklist before calculating frame counts\n") ;
		free(num_data) ;
		return ;
	}
	
	if (self->fcounts == NULL)
		self->fcounts = malloc(self->tot_num_data*sizeof(int)) ;
	if (self->mean_count == NULL)
		self->mean_count = malloc(self->num_det*sizeof(double)) ;
	if (self->rescale == NULL)
		self->rescale = malloc(self->num_det*sizeof(double)) ;

	memset(self->fcounts, 0, self->tot_num_data*sizeof(int)) ;
	memset(self->mean_count, 0, self->num_det*sizeof(double)) ;
	memset(self->rescale, 0, self->num_det*sizeof(double)) ;
	
	while (curr != NULL) {
		cdet = curr->det ;
		detn = self->det_mapping[dset_num] ;
		
		if (curr->ftype == SPARSE) {
			for (d = 0 ; d < curr->num_data ; ++d) {
				if (blacklist[curr->num_offset + d])
					continue ;
				for (t = 0 ; t < curr->ones[d] ; ++t)
				if (cdet->raw_mask[curr->place_ones[curr->ones_accum[d] + t]] < 1)
					self->fcounts[curr->num_offset + d]++ ;
				
				for (t = 0 ; t < curr->multi[d] ; ++t)
				if (cdet->raw_mask[curr->place_multi[curr->multi_accum[d] + t]] < 1)
					self->fcounts[curr->num_offset + d] += curr->count_multi[curr->multi_accum[d] + t] ;
				
				self->mean_count[detn] += self->fcounts[curr->num_offset + d] ;
				num_data[detn]++ ;
			}
		}
		else if (curr->ftype == DENSE_INT) {
			for (d = 0 ; d < curr->num_data ; ++d) {
				if (blacklist[curr->num_offset + d])
					continue ;
				for (t = 0 ; t < curr->num_pix ; ++t)
					self->fcounts[curr->num_offset + d] += curr->int_frames[d*curr->num_pix + t] ;
				
				self->mean_count[detn] += self->fcounts[curr->num_offset + d] ;
				num_data[detn]++ ;
			}
		}
		else if (curr->ftype == DENSE_DOUBLE) {
			for (d = 0 ; d < curr->num_data ; ++d) {
				if (blacklist[curr->num_offset + d])
					continue ;
				for (t = 0 ; t < curr->num_pix ; ++t)
					self->fcounts[curr->num_offset + d] += curr->frames[d*curr->num_pix + t] ;
				
				self->mean_count[detn] += self->fcounts[curr->num_offset + d] ;
				num_data[detn]++ ;
			}
		}
		
		dset_num++ ;
		curr = curr->next ;
	}
	
	for (d = 0 ; d < self->num_det ; ++d)
		self->mean_count[d] /= num_data[d] ;
	free(num_data) ;
}

void calc_frame_counts_partial(struct iterate *self, int rank, int num_proc, int *num_data, uint8_t *blacklist) {
	long t, d, detn, dset_num = 0 ;
	struct dataset *curr = self->dset ;
	struct detector *cdet ;
	if (blacklist == NULL) {
		fprintf(stderr, "Need to generate an empty blacklist before calculating frame counts\n") ;
		return ;
	}
	
	if (self->fcounts == NULL)
		self->fcounts = malloc(self->tot_num_data*sizeof(int)) ;
	if (self->mean_count == NULL)
		self->mean_count = malloc(self->num_det*sizeof(double)) ;
	if (self->rescale == NULL)
		self->rescale = malloc(self->num_det*sizeof(double)) ;

	memset(self->fcounts, 0, self->tot_num_data*sizeof(int)) ;
	memset(self->mean_count, 0, self->num_det*sizeof(double)) ;
	memset(self->rescale, 0, self->num_det*sizeof(double)) ;
	memset(num_data, 0, self->num_det*sizeof(int)) ;
	
	while (curr != NULL) {
		if (dset_num % num_proc != rank) {
			dset_num++ ;
			curr = curr->next ;
			continue ;
		}
		cdet = curr->det ;
		detn = self->det_mapping[dset_num] ;
		
		if (curr->ftype == SPARSE) {
			for (d = 0 ; d < curr->num_data ; ++d) {
				if (blacklist[curr->num_offset + d])
					continue ;
				for (t = 0 ; t < curr->ones[d] ; ++t)
				if (cdet->raw_mask[curr->place_ones[curr->ones_accum[d] + t]] < 1)
					self->fcounts[curr->num_offset + d]++ ;
				
				for (t = 0 ; t < curr->multi[d] ; ++t)
				if (cdet->raw_mask[curr->place_multi[curr->multi_accum[d] + t]] < 1)
					self->fcounts[curr->num_offset + d] += curr->count_multi[curr->multi_accum[d] + t] ;
				
				self->mean_count[detn] += self->fcounts[curr->num_offset + d] ;
				num_data[detn]++ ;
			}
		}
		else if (curr->ftype == DENSE_INT) {
			for (d = 0 ; d < curr->num_data ; ++d) {
				if (blacklist[curr->num_offset + d])
					continue ;
				for (t = 0 ; t < curr->num_pix ; ++t)
					self->fcounts[curr->num_offset + d] += curr->int_frames[d*curr->num_pix + t] ;
				
				self->mean_count[detn] += self->fcounts[curr->num_offset + d] ;
				num_data[detn]++ ;
			}
		}
		else if (curr->ftype == DENSE_DOUBLE) {
			for (d = 0 ; d < curr->num_data ; ++d) {
				if (blacklist[curr->num_offset + d])
					continue ;
				for (t = 0 ; t < curr->num_pix ; ++t)
					self->fcounts[curr->num_offset + d] += curr->frames[d*curr->num_pix + t] ;
				
				self->mean_count[detn] += self->fcounts[curr->num_offset + d] ;
				num_data[detn]++ ;
			}
		}
		
		dset_num++ ;
		curr = curr->next ;
	}
}

void calc_beta(double start, struct iterate *self) {
	int d ;

	self->beta_start = malloc(self->tot_num_data * sizeof(double)) ;
	self->beta = malloc(self->tot_num_data * sizeof(double)) ;
	
	if (!self->par->need_scaling && start < 0)
		start = exp(-6.5 * pow(self->mean_count[0] * 1.e-5, 0.15)) ; // Empirical
	
	if (start > 0) {
		for (d = 0 ; d < self->tot_num_data ; ++d)
			self->beta_start[d] = start ;
	}
	else {
		for (d = 0 ; d < self->tot_num_data ; ++d)
			self->beta_start[d] = exp(-6.5 * pow(self->fcounts[d] * 1.e-5, 0.15)) ; // Empirical
	}
}

void calc_sum_fact(struct iterate *self, uint8_t *blacklist) {
	int dset_num = 0, d, t ;
	struct dataset *curr = self->dset ;
	struct detector *det ;
	if (blacklist == NULL) {
		fprintf(stderr, "Need to generate an empty blacklist before calculating sum_fact\n") ;
		return ;
	}
	
	if (self->sum_fact == NULL)
		self->sum_fact = malloc(self->tot_num_data*sizeof(double)) ;
	memset(self->sum_fact, 0, self->tot_num_data*sizeof(double)) ;
	
	while (curr != NULL) {
		det = curr->det ;
		
		if (curr->ftype == SPARSE) {
			for (d = 0 ; d < curr->num_data ; ++d) {
				if (blacklist[d + curr->num_offset])
					continue ;
				for (t = 0 ; t < curr->multi[d] ; ++t)
				if (det->raw_mask[curr->place_multi[curr->multi_accum[d] + t]] < 1)
					self->sum_fact[curr->num_offset+d] += gsl_sf_lnfact(curr->count_multi[curr->multi_accum[d] + t]) ;
			}
		}
		else if (curr->ftype == DENSE_INT) {
			for (d = 0 ; d < curr->num_data ; ++d) {
				if (blacklist[d + curr->num_offset])
					continue ;
				for (t = 0 ; t < curr->num_pix ; ++t)
				if (det->raw_mask[t] < 1)
					self->sum_fact[curr->num_offset+d] += gsl_sf_lnfact(curr->int_frames[d*curr->num_pix + t]) ;
			}
		}
		else if (curr->ftype == DENSE_DOUBLE) {
			for (d = 0 ; d < curr->num_data ; ++d)
				self->sum_fact[curr->num_offset+d] = 0. ;
		}
		
		dset_num++ ;
		curr = curr->next ;
	}
}

void calc_sum_fact_partial(struct iterate *self, int rank, int num_proc, uint8_t *blacklist) {
	int dset_num = 0, d, t ;
	struct dataset *curr = self->dset ;
	struct detector *det ;
	if (blacklist == NULL) {
		fprintf(stderr, "Need to generate an empty blacklist before calculating sum_fact\n") ;
		return ;
	}
	
	if (self->sum_fact == NULL)
		self->sum_fact = malloc(self->tot_num_data*sizeof(double)) ;
	memset(self->sum_fact, 0, self->tot_num_data*sizeof(double)) ;
	
	while (curr != NULL) {
		if (dset_num % num_proc != rank) {
			dset_num++ ;
			curr = curr->next ;
			continue ;
		}
		det = curr->det ;
		
		if (curr->ftype == SPARSE) {
			for (d = 0 ; d < curr->num_data ; ++d) {
				if (blacklist[d + curr->num_offset])
					continue ;
				for (t = 0 ; t < curr->multi[d] ; ++t)
				if (det->raw_mask[curr->place_multi[curr->multi_accum[d] + t]] < 1)
					self->sum_fact[curr->num_offset+d] += gsl_sf_lnfact(curr->count_multi[curr->multi_accum[d] + t]) ;
			}
		}
		else if (curr->ftype == DENSE_INT) {
			for (d = 0 ; d < curr->num_data ; ++d) {
				if (blacklist[d + curr->num_offset])
					continue ;
				for (t = 0 ; t < curr->num_pix ; ++t)
				if (det->raw_mask[t] < 1)
					self->sum_fact[curr->num_offset+d] += gsl_sf_lnfact(curr->int_frames[d*curr->num_pix + t]) ;
			}
		}
		else if (curr->ftype == DENSE_DOUBLE) {
			for (d = 0 ; d < curr->num_data ; ++d)
				self->sum_fact[curr->num_offset+d] = 0. ;
		}
		
		dset_num++ ;
		curr = curr->next ;
	}
}

void calc_powder(struct iterate *self, uint8_t *blacklist) {
	int dset_num = 0, detn, d, t, pixel ;
	struct dataset *curr = self->dset ;
	struct detector *det ;
	int *nframes = calloc(self->num_det, sizeof(int)) ;
	if (blacklist == NULL) {
		fprintf(stderr, "Need to generate an empty blacklist before calculating powder\n") ;
		free(nframes) ;
		return ;
	}
	
	// Only allocate powder if detector has background
	if (self->det[0].with_bg) {
		for (detn = 0 ; detn < self->num_det ; ++detn) {
			if (self->det[detn].powder != NULL)
				free(self->det[detn].powder) ;
			self->det[detn].powder = calloc(self->det[detn].num_pix, sizeof(double)) ;
		}
		
		while (curr != NULL) {
			detn = self->det_mapping[dset_num] ;
			det = &self->det[detn] ;
			
			if (curr->ftype == SPARSE) {
				for (d = 0 ; d < curr->num_data ; ++d) {
					if (blacklist[curr->num_offset + d])
						continue ;
					nframes[detn]++ ;
					for (t = 0 ; t < curr->ones[d] ; ++t) {
						pixel = curr->place_ones[curr->ones_accum[d] + t] ;
						if (det->raw_mask[pixel] < 1)
							det->powder[pixel]++ ;
					}
					for (t = 0 ; t < curr->multi[d] ; ++t) {
						pixel = curr->place_multi[curr->multi_accum[d] + t] ;
						if (det->raw_mask[pixel] < 1)
							det->powder[pixel] += curr->count_multi[curr->multi_accum[d] + t] ;
					}
				}
			}
			else if (curr->ftype == DENSE_INT) {
				for (d = 0 ; d < curr->num_data ; ++d) {
					if (blacklist[curr->num_offset + d])
						continue ;
					nframes[detn]++ ;
					for (t = 0 ; t < curr->num_pix ; ++t)
						if (det->raw_mask[t] < 1)
							det->powder[t] += curr->int_frames[d*curr->num_pix + t] ;
				}
			}
			else if (curr->ftype == DENSE_DOUBLE) {
				for (d = 0 ; d < curr->num_data ; ++d) {
					if (blacklist[curr->num_offset + d])
						continue ;
					nframes[detn]++ ;
					for (t = 0 ; t < curr->num_pix ; ++t)
						if (det->raw_mask[t] < 1)
							det->powder[t] += curr->frames[d*curr->num_pix + t] ;
				}
			}
			
			dset_num++ ;
			curr = curr->next ;
		}
		
		for (detn = 0 ; detn < self->num_det ; ++detn)
			for (t = 0 ; t < self->det[detn].num_pix ; ++t)
				self->det[detn].powder[t] /= nframes[detn] ;
	}
	
	free(nframes) ;
}

void calc_powder_partial(struct iterate *self, int rank, int num_proc, int *nframes, uint8_t *blacklist) {
	int dset_num = 0, detn, d, t, pixel ;
	struct dataset *curr = self->dset ;
	struct detector *det ;
	if (blacklist == NULL) {
		fprintf(stderr, "Need to generate an empty blacklist before calculating powder\n") ;
		return ;
	}
	memset(nframes, 0, self->num_det*sizeof(int)) ;
	
	// Only allocate powder if detector has background
	if (self->det[0].with_bg) {
		for (detn = 0 ; detn < self->num_det ; ++detn) {
			if (self->det[detn].powder != NULL)
				free(self->det[detn].powder) ;
			self->det[detn].powder = calloc(self->det[detn].num_pix, sizeof(double)) ;
		}
		
		while (curr != NULL) {
			if (dset_num % num_proc != rank) {
				dset_num++ ;
				curr = curr->next ;
				continue ;
			}
			detn = self->det_mapping[dset_num] ;
			det = &self->det[detn] ;
			
			if (curr->ftype == SPARSE) {
				for (d = 0 ; d < curr->num_data ; ++d) {
					if (blacklist[curr->num_offset + d])
						continue ;
					nframes[detn]++ ;
					for (t = 0 ; t < curr->ones[d] ; ++t) {
						pixel = curr->place_ones[curr->ones_accum[d] + t] ;
						if (det->raw_mask[pixel] < 1)
							det->powder[pixel]++ ;
					}
					for (t = 0 ; t < curr->multi[d] ; ++t) {
						pixel = curr->place_multi[curr->multi_accum[d] + t] ;
						if (det->raw_mask[pixel] < 1)
							det->powder[pixel] += curr->count_multi[curr->multi_accum[d] + t] ;
					}
				}
			}
			else if (curr->ftype == DENSE_INT) {
				for (d = 0 ; d < curr->num_data ; ++d) {
					if (blacklist[curr->num_offset + d])
						continue ;
					nframes[detn]++ ;
					for (t = 0 ; t < curr->num_pix ; ++t)
						if (det->raw_mask[t] < 1)
							det->powder[t] += curr->int_frames[d*curr->num_pix + t] ;
				}
			}
			else if (curr->ftype == DENSE_DOUBLE) {
				for (d = 0 ; d < curr->num_data ; ++d) {
					if (blacklist[curr->num_offset + d])
						continue ;
					nframes[detn]++ ;
					for (t = 0 ; t < curr->num_pix ; ++t)
						if (det->raw_mask[t] < 1)
							det->powder[t] += curr->frames[d*curr->num_pix + t] ;
				}
			}
			
			dset_num++ ;
			curr = curr->next ;
		}
	}
}
