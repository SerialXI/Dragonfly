from libc.stdint cimport uint8_t
from .detector cimport detector

cdef extern from "src/emcfile.h" nogil:
    enum frame_type: SPARSE, DENSE_INT, DENSE_DOUBLE
    struct dataset:
        char *fname
        frame_type ftype
        int num_data, num_pix
        double mean_count
        detector *det

        # Linked list information
        dataset *next
        int num_offset

        # Sparse data
        int *ones
        int *multi
        int *ones_file
        int *multi_file
        int *place_ones
        int *place_multi
        int *count_multi
        long ones_total, multi_total
        long *ones_accum
        long *multi_accum
        long ones_total_file, multi_total_file
        long *ones_accum_file
        long *multi_accum_file
        int lazy

        # Dense data
        int *int_frames
        double *frames

    int parse_dataset(char*, detector*, dataset*, int)
    int load_active_frames(dataset*, uint8_t*)

cdef class CDataset:
    cdef dataset *dset
