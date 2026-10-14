#include "iterate.h"

void calc_frame_counts(struct iterate *self) {
	long t, d, detn, dset_num = 0 ;
	struct dataset *curr = self->dset ;
	struct detector *cdet ;
	int *num_data = calloc(self->num_det, sizeof(int)) ;
	if (self->static_blacklist == NULL) {
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
				if (self->static_blacklist[curr->num_offset + d])
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
				if (self->static_blacklist[curr->num_offset + d])
					continue ;
				for (t = 0 ; t < curr->num_pix ; ++t)
					self->fcounts[curr->num_offset + d] += curr->int_frames[d*curr->num_pix + t] ;
				
				self->mean_count[detn] += self->fcounts[curr->num_offset + d] ;
				num_data[detn]++ ;
			}
		}
		else if (curr->ftype == DENSE_DOUBLE) {
			for (d = 0 ; d < curr->num_data ; ++d) {
				if (self->static_blacklist[curr->num_offset + d])
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

void calc_sum_fact(struct iterate *self) {
	int dset_num = 0, d, t ;
	struct dataset *curr = self->dset ;
	struct detector *det ;
	if (self->static_blacklist == NULL) {
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
				if (self->static_blacklist[d + curr->num_offset])
					continue ;
				for (t = 0 ; t < curr->multi[d] ; ++t)
				if (det->raw_mask[curr->place_multi[curr->multi_accum[d] + t]] < 1)
					self->sum_fact[curr->num_offset+d] += gsl_sf_lnfact(curr->count_multi[curr->multi_accum[d] + t]) ;
			}
		}
		else if (curr->ftype == DENSE_INT) {
			for (d = 0 ; d < curr->num_data ; ++d) {
				if (self->static_blacklist[d + curr->num_offset])
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

void calc_powder(struct iterate *self) {
	int dset_num = 0, detn, d, t, pixel ;
	struct dataset *curr = self->dset ;
	struct detector *det ;
	int *nframes = calloc(self->num_det, sizeof(int)) ;
	if (self->static_blacklist == NULL) {
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
					if (self->static_blacklist[curr->num_offset + d])
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
					if (self->static_blacklist[curr->num_offset + d])
						continue ;
					nframes[detn]++ ;
					for (t = 0 ; t < curr->num_pix ; ++t)
						if (det->raw_mask[t] < 1)
							det->powder[t] += curr->int_frames[d*curr->num_pix + t] ;
				}
			}
			else if (curr->ftype == DENSE_DOUBLE) {
				for (d = 0 ; d < curr->num_data ; ++d) {
					if (self->static_blacklist[curr->num_offset + d])
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

void free_rel_quat(struct iterate *iter) {
    int d ;
    for (d = 0 ; d < iter->tot_num_data ; ++d) {
        if (iter->rel_quat != NULL)
            free(iter->rel_quat[d]) ;
        if (iter->rel_prob != NULL)
            free(iter->rel_prob[d]) ;
    }
    free(iter->rel_quat) ;
    free(iter->rel_prob) ;
    free(iter->num_rel_quat) ;
    iter->rel_quat = NULL ;
    iter->rel_prob = NULL ;
    iter->num_rel_quat = NULL ;
}

static int read_rel_quat_h5(hid_t file, int num_rot_coarse, int parse_prob, struct iterate *iter) {
    int d, k, value, status ;
    hid_t dset, dspace, dtype ;
    hsize_t num_data ;
    const char *keys[] = {"params/num_div", "params/num_rot", "params/num_modes", "params/nonrot_modes"} ;
    int expected[] = {iter->par->coarse_div, num_rot_coarse, iter->par->num_modes, iter->par->nonrot_modes} ;

    if (H5Lexists(file, "params", H5P_DEFAULT) > 0)
    for (k = 0 ; k < 4 ; ++k) {
        if (H5Lexists(file, keys[k], H5P_DEFAULT) <= 0)
            continue ;
        dset = H5Dopen2(file, keys[k], H5P_DEFAULT) ;
        status = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) ;
        H5Dclose(dset) ;
        if (status < 0 || value != expected[k]) {
            fprintf(stderr, "Coarse checkpoint does not match refinement %s\n", keys[k]) ;
            return 1 ;
        }
    }

    int sparse = H5Lexists(file, "probabilities", H5P_DEFAULT) > 0 ;
    if (sparse) {
        dset = H5Dopen2(file, "probabilities/num_rot", H5P_DEFAULT) ;
        if (dset < 0)
            return 1 ;
        status = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) ;
        H5Dclose(dset) ;
        if (status < 0 || value != num_rot_coarse) {
            fprintf(stderr, "Rotation sampling in probs file does not match num_rot_coarse\n") ;
            return 1 ;
        }
    }

    const char *paths[] = {"probabilities/place", "probabilities/prob"} ;
    for (k = 0 ; k < (sparse && parse_prob ? 2 : 1) ; ++k) {
        dset = H5Dopen2(file, sparse ? paths[k] : "orientations", H5P_DEFAULT) ;
        if (dset < 0)
            return 1 ;
        dspace = H5Dget_space(dset) ;
        if (H5Sget_simple_extent_ndims(dspace) != 1 ||
                H5Sget_simple_extent_dims(dspace, &num_data, NULL) < 0 || num_data != (hsize_t) iter->tot_num_data) {
            fprintf(stderr, "Number of frames in coarse checkpoint does not match %d\n", iter->tot_num_data) ;
            H5Sclose(dspace) ;
            H5Dclose(dset) ;
            return 1 ;
        }

        if (!sparse) {
            int *buffer = malloc(num_data * sizeof(int)) ;
            status = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) ;
            if (status >= 0)
            for (d = 0 ; d < iter->tot_num_data ; ++d) {
                iter->num_rel_quat[d] = 1 ;
                iter->rel_quat[d] = malloc(sizeof(int)) ;
                iter->rel_quat[d][0] = buffer[d] ;
                if (parse_prob) {
                    iter->rel_prob[d] = malloc(sizeof(double)) ;
                    iter->rel_prob[d][0] = 1. ;
                }
            }
            free(buffer) ;
        }
        else {
            hvl_t *buffer = calloc(num_data, sizeof(hvl_t)) ;
            dtype = H5Tvlen_create(k == 0 ? H5T_NATIVE_INT : H5T_NATIVE_DOUBLE) ;
            status = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) ;
            if (status >= 0)
            for (d = 0 ; d < iter->tot_num_data ; ++d) {
                if (k == 0) {
                    iter->num_rel_quat[d] = buffer[d].len ;
                    if (buffer[d].len > 0) {
                        iter->rel_quat[d] = malloc(buffer[d].len * sizeof(int)) ;
                        memcpy(iter->rel_quat[d], buffer[d].p, buffer[d].len * sizeof(int)) ;
                    }
                }
                else {
                    if (buffer[d].len != (size_t) iter->num_rel_quat[d]) {
                        fprintf(stderr, "Probability count does not match orientations for frame %d\n", d) ;
                        status = -1 ;
                        break ;
                    }
                    if (buffer[d].len > 0) {
                        iter->rel_prob[d] = malloc(buffer[d].len * sizeof(double)) ;
                        memcpy(iter->rel_prob[d], buffer[d].p, buffer[d].len * sizeof(double)) ;
                    }
                }
            }
            H5Dvlen_reclaim(dtype, dspace, H5P_DEFAULT, buffer) ;
            free(buffer) ;
            H5Tclose(dtype) ;
        }
        H5Sclose(dspace) ;
        H5Dclose(dset) ;
        if (status < 0)
            return 1 ;
    }
    return 0 ;
}

static int read_rel_quat_binary(FILE *fp, int num_rot_coarse, int parse_prob, struct iterate *iter) {
    int d, num_data, num_rot ;
    if (fread(&num_data, sizeof(int), 1, fp) != 1 || fread(&num_rot, sizeof(int), 1, fp) != 1)
        return 1 ;
    if (num_data != iter->tot_num_data || num_rot != num_rot_coarse) {
        fprintf(stderr, "Frame count or rotation sampling in probs file does not match refinement\n") ;
        return 1 ;
    }
    if (fseek(fp, 1024, SEEK_SET) != 0 ||
            fread(iter->num_rel_quat, sizeof(int), num_data, fp) != (size_t) num_data)
        return 1 ;
    for (d = 0 ; d < num_data ; ++d) {
        int n = iter->num_rel_quat[d] ;
        if (n < 0)
            return 1 ;
        if (n == 0)
            continue ;
        iter->rel_quat[d] = malloc(n * sizeof(int)) ;
        if (fread(iter->rel_quat[d], sizeof(int), n, fp) != (size_t) n)
            return 1 ;
    }
    if (parse_prob)
    for (d = 0 ; d < num_data ; ++d) {
        int n = iter->num_rel_quat[d] ;
        if (n == 0)
            continue ;
        iter->rel_prob[d] = malloc(n * sizeof(double)) ;
        if (fread(iter->rel_prob[d], sizeof(double), n, fp) != (size_t) n)
            return 1 ;
    }
    return 0 ;
}

static int restrict_rel_quat(int num_rot_coarse, int parse_prob, struct iterate *iter) {
    int d, t, i, n, index, rot ;
    int num_modes = iter->par->num_modes ;
    for (d = 0 ; d < iter->tot_num_data ; ++d) {
        n = 0 ;
        for (t = 0 ; t < iter->num_rel_quat[d] ; ++t) {
            index = iter->rel_quat[d][t] ;
            if (index < 0 || index >= num_rot_coarse * num_modes)
                continue ;
            // Saved indices include modes; refinement restricts only orientations.
            rot = index / num_modes ;
            i = n ;
            if (num_modes > 1)
            for (i = 0 ; i < n ; ++i)
                if (iter->rel_quat[d][i] == rot)
                    break ;
            if (i == n) {
                iter->rel_quat[d][n] = rot ;
                if (parse_prob)
                    iter->rel_prob[d][n] = iter->rel_prob[d][t] ;
                n++ ;
            }
            else if (parse_prob)
                iter->rel_prob[d][i] += iter->rel_prob[d][t] ;
        }
        iter->num_rel_quat[d] = n ;
        if (n == 0 && !iter->static_blacklist[d] && iter->par->nonrot_modes == 0) {
            fprintf(stderr, "No coarse orientations for active frame %d\n", d) ;
            return 1 ;
        }
    }
    return 0 ;
}

int parse_rel_quat(const char *fname, int num_rot_coarse, int parse_prob, struct iterate *iter) {
    FILE *fp = fopen(fname, "rb") ;
    if (fp == NULL) {
        fprintf(stderr, "Could not open rel_quat file %s\n", fname) ;
        return 1 ;
    }
    // Keep the previous restrictions if the new checkpoint cannot be read.
    struct iterate data = *iter ;
    data.num_rel_quat = calloc(iter->tot_num_data, sizeof(int)) ;
    data.rel_quat = calloc(iter->tot_num_data, sizeof(int*)) ;
    data.rel_prob = parse_prob ? calloc(iter->tot_num_data, sizeof(double*)) : NULL ;
    int status ;
    if (H5Fis_hdf5(fname) > 0) {
        fclose(fp) ;
        hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT) ;
        status = file < 0 ? 1 : read_rel_quat_h5(file, num_rot_coarse, parse_prob, &data) ;
        if (file >= 0)
            H5Fclose(file) ;
    }
    else {
        status = read_rel_quat_binary(fp, num_rot_coarse, parse_prob, &data) ;
        fclose(fp) ;
    }
    if (status == 0)
        status = restrict_rel_quat(num_rot_coarse, parse_prob, &data) ;
    if (status != 0) {
        free_rel_quat(&data) ;
        return 1 ;
    }
    free_rel_quat(iter) ;
    iter->num_rel_quat = data.num_rel_quat ;
    iter->rel_quat = data.rel_quat ;
    iter->rel_prob = data.rel_prob ;
    return 0 ;
}
