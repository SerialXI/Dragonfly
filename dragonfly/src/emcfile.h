#ifndef EMCFILE_H
#define EMCFILE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <hdf5.h>
#include "detector.h"

enum frame_type {SPARSE, DENSE_INT, DENSE_DOUBLE} ;

struct dataset {
	char *fname ;
    enum frame_type ftype ;
	int num_data, num_pix ;
	int has_binary_magic ;
	double mean_count ;
	struct detector *det ;

	// Linked list information
	struct dataset *next ;
	int num_offset ;

	// Sparse data
	int *ones, *multi ;
	int *ones_file, *multi_file ;
	int *place_ones, *place_multi, *count_multi ;
	long ones_total, multi_total ;
	long *ones_accum, *multi_accum ;
	long ones_total_file, multi_total_file ;
	long *ones_accum_file, *multi_accum_file ;
	int lazy ;

	// Dense data
	int *int_frames ;
	double *frames ;
} ;

int parse_dataset(char*, struct detector*, struct dataset*, int) ;
int load_active_frames(struct dataset*, uint8_t*) ;

#endif // EMCFILE_H
