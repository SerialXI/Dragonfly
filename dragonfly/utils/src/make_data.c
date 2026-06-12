#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <limits.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <sys/time.h>
#include <float.h>
#include <stdint.h>
#include <libgen.h>
#include <hdf5.h>
#include "../../src/detector.h"
#include "../../src/model.h"

#define NUM_AVE 5000
#define TESTING 0

void rand_quat(double quat[4], gsl_rng *rng) {
	int i ;
	double qq ;
	
	do {
		qq = 0. ;
		for (i = 0 ; i < 4 ; ++i) {
			quat[i] = gsl_rng_uniform(rng) -.5 ;
			qq += quat[i] * quat[i] ;
		}
	}
	while (qq > .25) ;
	
	qq = sqrt(qq) ;
	for (i = 0 ; i < 4 ; ++i)
		quat[i] /= qq ;
}

double rescale_intens(double fluence, double mean_count,
                      struct model *mod, struct detector *det) {
	int x ;
	double rescale = 0., intens_ave = 0. ;
	const gsl_rng_type *T = gsl_rng_default ;
	gsl_rng *rng = gsl_rng_alloc(T) ;
	unsigned long *seeds = malloc(omp_get_max_threads() * sizeof(unsigned long)) ;
	
	if (TESTING) {
		gsl_rng_set(rng, 0x5EED) ;
	}
	else {
		struct timeval tval ;
		gettimeofday(&tval, NULL) ;
		gsl_rng_set(rng, tval.tv_sec + tval.tv_usec) ;
	}
	for (x = 0 ; x < omp_get_max_threads() ; ++x)
		seeds[x] = gsl_rng_get(rng) ;
	
	#pragma omp parallel default(shared)
	{
		int d, t, rank = omp_get_thread_num() ;
		double quat[4] ;
		double *view = malloc(det->num_pix * sizeof(double)) ;
		gsl_rng *rng = gsl_rng_alloc(T) ;
		gsl_rng_set(rng, seeds[rank]) ;
		
		#pragma omp for schedule(static) reduction(+:intens_ave)
		for (d = 0 ; d < NUM_AVE ; ++d) {
			rand_quat(quat, rng) ;
			slice_gen3d(quat, 0, view, det, mod) ;
			
			for (t = 0 ; t < det->num_pix ; ++t){
				if (det->raw_mask[t] > 1)
					continue ;
				intens_ave += view[t] ;
			}
		}
		
		free(view) ;
		gsl_rng_free(rng) ;
	}
	
	free(seeds) ;
	gsl_rng_free(rng) ;
	intens_ave /= NUM_AVE ;
	
	if (fluence > 0. && mean_count > 0.) {
		fprintf(stderr, "ERROR: Only one of fluence and mean_count can be positive\n") ;
		return -1. ;
	}
	
	if (fluence < 0. && mean_count < 0.) {
		fprintf(stderr, "ERROR: Need to specify either fluence or mean_count for scaling\n") ;
		return -1. ;
	}
	
	if (fluence > 0.) {
		rescale = fluence*pow(2.81794e-9, 2) ;
		mean_count = rescale*intens_ave ;
		fprintf(stderr, "Target mean_count = %f for fluence = %.3e photons/um^2\n", mean_count, fluence) ;
	}
	else if (mean_count > 0.) {
		rescale = mean_count / intens_ave ;
		fprintf(stderr, "Using target mean_count of %f photons/frame\n", mean_count) ;
	}
	
	for (x = 0 ; x < mod->vol ; ++x)
		mod->model1[x] *= rescale ;
	
	return mean_count ;
}

double gen_and_save_dataset(long num_data, double mean_count, long do_gamma,
                            char *output_fname, long h5_output,
							char *likelihood_fname, char *scale_fname,
                            struct model *mod, struct detector *det) {
	long x, d, t ;
	double bg_count = 0., actual_mean_count = 0. ;
	int num_counts[2] ;
	const gsl_rng_type *T = gsl_rng_default ;
	gsl_rng *rng = gsl_rng_alloc(T) ;
	unsigned long *seeds = malloc(omp_get_max_threads() * sizeof(unsigned long)) ;
	int **place_ones, **place_multi, *ones, *multi, **count_multi ;
	double *likelihood, *scale_factors ;
	
	ones = calloc(num_data, sizeof(int)) ;
	multi = calloc(num_data, sizeof(int)) ;
	place_ones = malloc(num_data * sizeof(int*)) ;
	place_multi = malloc(num_data * sizeof(int*)) ;
	count_multi = malloc(num_data * sizeof(int*)) ;
	likelihood = calloc(num_data, sizeof(double)) ;
	scale_factors = malloc(num_data * sizeof(double)) ;
    if (det->with_bg) {
        for (t = 0 ; t < det->num_pix ; ++t)
            bg_count += det->background[t] ;
	}
	
	num_counts[1] = (mean_count + bg_count) > det->num_pix ? det->num_pix : (mean_count + bg_count) ;
	num_counts[0] = 10*num_counts[1] > det->num_pix ? det->num_pix : 10*num_counts[1];
	//fprintf(stderr, "Assuming maximum of %d and %d ones and multi pixels respectively.\n", num_counts[0], num_counts[1]) ;
	
	for (d = 0 ; d < num_data ; ++d) {
		place_ones[d] = malloc((size_t) num_counts[0] * sizeof(int)) ;
		place_multi[d] = malloc((size_t) num_counts[1] * sizeof(int)) ;
		count_multi[d] = malloc((size_t) num_counts[1] * sizeof(int)) ;
	}

	if (TESTING) {
		gsl_rng_set(rng, 0x5EED) ;
	}
	else {
		struct timeval tval ;
		gettimeofday(&tval, NULL) ;
		gsl_rng_set(rng, tval.tv_sec + tval.tv_usec) ;
	}
	for (x = 0 ; x < omp_get_max_threads() ; ++x)
		seeds[x] = gsl_rng_get(rng) ;
	
	#pragma omp parallel default(shared)
	{
        long d, t ;
		int photons, rank = omp_get_thread_num() ;
		int curr_counts[2] ;
		double scale = 1., quat[4], val ;
		double *view = malloc(det->num_pix * sizeof(double)) ;
		gsl_rng *rng = gsl_rng_alloc(T) ;
		gsl_rng_set(rng, seeds[rank]) ;
		
		#pragma omp for schedule(static,1) reduction(+:actual_mean_count)
		for (d = 0 ; d < num_data ; ++d) {
			rand_quat(quat, rng) ;
			slice_gen3d(quat, 0, view, det, mod) ;
			curr_counts[0] = num_counts[0] ;
			curr_counts[1] = num_counts[1] ;
			
			if (do_gamma)
				scale = gsl_ran_gamma(rng, 2., 0.5) ;
			
			if (scale > 0.) {
				for (t = 0 ; t < det->num_pix ; ++t) {
					if (det->raw_mask[t] > 1)
						continue ;
					
					val = view[t]*scale ;
					if (det->with_bg)
						val += det->background[t] ;
					photons = gsl_ran_poisson(rng, val) ;
					
					if (photons == 1) {
						place_ones[d][ones[d]++] = t ;
					}
					else if (photons > 1) {
						place_multi[d][multi[d]] = t ;
						count_multi[d][multi[d]++] = photons ;
						actual_mean_count += photons ;
					}
					
					if (likelihood_fname[0] != '\0') {
						if (photons == 0)
							likelihood[d] -= val ;
						else
							likelihood[d] += photons*log(val) - val - gsl_sf_lnfact(photons) ;
					}
					if (scale_fname[0] != '\0')
						scale_factors[d] = scale ;
					if (ones[d] >= curr_counts[0]) {
						curr_counts[0] *= 2 ;
						place_ones[d] = realloc(place_ones[d], curr_counts[0]*sizeof(int)) ;
					}
					if (multi[d] >= curr_counts[1]) {
						curr_counts[1] *= 2 ;
						place_multi[d] = realloc(place_multi[d], curr_counts[1]*sizeof(int)) ;
						count_multi[d] = realloc(count_multi[d], curr_counts[1]*sizeof(int)) ;
					}
				}
			}
			
			actual_mean_count += ones[d] ;
			
			if (rank == 0)
				fprintf(stderr, "\rFinished d = %ld", d) ;
		}
		
		free(view) ;
		gsl_rng_free(rng) ;
	}
	 
	free(seeds) ;
	gsl_rng_free(rng) ;
	fprintf(stderr, "\rFinished d = %ld\n", num_data) ;
	actual_mean_count /= num_data ;

	if (h5_output == 0) {
		int d, header[256] = {0} ;
		unsigned char emc_magic[8] = {137, 'E', 'M', 'C', '\r', '\n', 26, '\n'} ;
		header[0] = num_data ;
		header[1] = det->num_pix ;
		memcpy((char*) header + 12, emc_magic, 8) ;
		
		FILE *fp = fopen(output_fname, "wb") ;
		fwrite(header, sizeof(int), 256, fp) ;
		fwrite(ones, sizeof(int), num_data, fp) ;
		fwrite(multi, sizeof(int), num_data, fp) ;
		for (d = 0 ; d < num_data ; ++d)
			fwrite(place_ones[d], sizeof(int), ones[d], fp) ;
		for (d = 0 ; d < num_data ; ++d)
			fwrite(place_multi[d], sizeof(int), multi[d], fp) ;
		for (d = 0 ; d < num_data ; ++d)
			fwrite(count_multi[d], sizeof(int), multi[d], fp) ;
		fclose(fp) ;
		
		if (likelihood_fname[0] != '\0') {
			fp = fopen(likelihood_fname, "wb") ;
			fwrite(likelihood, sizeof(double), num_data, fp) ;
			fclose(fp) ;
		}
		if (scale_fname[0] != '\0') {
			fp = fopen(scale_fname, "w") ;
			for (d = 0 ; d < num_data ; ++d)
				fprintf(fp, "%13.10f\n", scale_factors[d]) ;
			fclose(fp) ;
		}
	}
	else {
		int d ;
		hid_t file, dset, dspace, dtype ;
		hsize_t dsize[1] = {1} ;
		hvl_t *po, *pm, *cm;
		
		file = H5Fcreate(output_fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT) ;
		dspace = H5Screate_simple(1, dsize, NULL) ;
		
		dset = H5Dcreate(file, "num_pix", H5T_STD_I32LE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
		H5Dwrite(dset, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(det->num_pix)) ;
		H5Dclose(dset) ;
		
		dsize[0] = num_data ;
		dspace = H5Screate_simple(1, dsize, NULL) ;
		dtype = H5Tvlen_create(H5T_STD_I32LE) ;
		po = malloc(num_data * sizeof(hvl_t)) ;
		pm = malloc(num_data * sizeof(hvl_t)) ;
		cm = malloc(num_data * sizeof(hvl_t)) ;
		for (d = 0 ; d < num_data ; ++d) {
			po[d].len = ones[d] ;
			po[d].p = place_ones[d] ;
			pm[d].len = multi[d] ;
			pm[d].p = place_multi[d] ;
			cm[d].len = multi[d] ;
			cm[d].p = count_multi[d] ;
		}
		dset = H5Dcreate(file, "place_ones", dtype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
		H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, po) ;
		H5Dclose(dset) ;
		free(po) ;
		
		dset = H5Dcreate(file, "place_multi", dtype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
		H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, pm) ;
		H5Dclose(dset) ;
		free(pm) ;
		
		dset = H5Dcreate(file, "count_multi", dtype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
		H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cm) ;
		H5Dclose(dset) ;
		free(cm) ;
		
		H5Sclose(dspace) ;
		H5Tclose(dtype) ;
		H5Fclose(file) ;
	}

	free(likelihood) ;
	free(scale_factors) ;
	
	free(ones) ;
	free(multi) ;
	for (d = 0 ; d < num_data ; ++d) {
		free(place_ones[d]) ;
		free(place_multi[d]) ;
		free(count_multi[d]) ;
	}
	free(place_ones) ;
	free(place_multi) ;
	free(count_multi) ;

	return actual_mean_count ;
}

