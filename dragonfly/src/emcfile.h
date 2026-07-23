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
	struct detector *det ;

	// Linked list information
	struct dataset *next ;
	int num_offset ;

	// Sparse data
	int *ones, *multi ;
	int *place_ones, *place_multi, *count_multi ;
	long ones_total, multi_total ;
	long *ones_accum, *multi_accum ;

	// Dense data
	int *int_frames ;
	double *frames ;
} ;

int parse_dataset(char*, struct detector*, struct dataset*) ;

#endif // EMCFILE_H
