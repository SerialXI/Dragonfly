#ifndef ITERATE_H
#define ITERATE_H

#include <gsl/gsl_sf_gamma.h>
#include "detector.h"
#include "emcfile.h"
#include "quaternion.h"
#include "model.h"
#include "params.h"

struct iterate {
	struct detector *det ;
	struct model *mod ;
	struct dataset *dset ;
	struct quaternion *quat ;
	struct params *par ;

	// Parameters for each frame
	int tot_num_data, num_blacklist ;
	int *fcounts ;
	int *last_selected ;
	double *scale, *bgscale ;
	double *beta, *beta_start ;
	double *sum_fact ;
	uint8_t *blacklist, *static_blacklist ;

	// For refinement
	int *quat_mapping, **rel_quat, *num_rel_quat ;
	double **rel_prob ;

	// Parameters for each detector
	int num_det ; // Number of unique detectors
	int num_dfiles ; // Number of datasets in linked list
	int *det_mapping ; // Mapping to unique list
	double *rescale, *mean_count ;

	// Aggregate metrics
	double likelihood, mutual_info, rms_change ;
} ;

void calc_frame_counts(struct iterate*, uint8_t*) ;
void calc_frame_counts_partial(struct iterate*, int, int, int*, uint8_t*) ;
void calc_beta(double, struct iterate*) ;
void calc_sum_fact(struct iterate*, uint8_t*) ;
void calc_sum_fact_partial(struct iterate*, int, int, uint8_t*) ;
void calc_powder(struct iterate*, uint8_t*) ;
void calc_powder_partial(struct iterate*, int, int, int*, uint8_t*) ;

#endif // ITERATE_H
