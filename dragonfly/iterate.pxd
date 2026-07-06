from libc.stdint cimport uint8_t
from .detector cimport detector
from .model cimport model
from .emcfile cimport dataset
from .quaternion cimport quaternion
from .params cimport params

cdef extern from "src/iterate.h" nogil:
    struct iterate:
        detector *det
        model *mod
        dataset *dset
        quaternion *quat
        params *par

        # Parameters for each frame
        int tot_num_data, num_blacklist
        int *fcounts
        int *last_selected
        double *scale
        double *bgscale
        double *beta
        double *beta_start
        double *sum_fact
        uint8_t *blacklist
        uint8_t *static_blacklist

        # For refinement
        int *quat_mapping
        int **rel_quat
        int *num_rel_quat
        double **rel_prob
    
        # Parameters for each detector
        int num_det # Number of unique detectors
        int num_dfiles # Number of datasets in linked list
        int *det_mapping
        double *rescale
        double *mean_count

        # Aggregate metrics
        double likelihood, mutual_info, rms_change

    void calc_frame_counts(iterate*)
    void calc_frame_counts_partial(iterate*, int, int, int*)
    void calc_beta(double, iterate*)
    void calc_sum_fact(iterate*)
    void calc_sum_fact_partial(iterate*, int, int)
    void calc_powder(iterate*)
    void calc_powder_partial(iterate*, int, int, int*)

cdef class Iterate:
    cdef iterate *iter
