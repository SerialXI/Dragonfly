from .detector cimport detector

cdef extern from "src/emcfile.h" nogil:
    enum frame_type: SPARSE, DENSE_INT, DENSE_DOUBLE
    struct dataset:
        char *fname
        frame_type ftype
        int num_data, num_pix
        int has_binary_magic
        detector *det

        # Linked list information
        dataset *next
        int num_offset

        # Sparse data
        int *ones
        int *multi
        int *place_ones
        int *place_multi
        int *count_multi
        long ones_total, multi_total
        long *ones_accum
        long *multi_accum

        # Dense data
        int *int_frames
        double *frames

    int parse_dataset(char*, detector*, dataset*)

cdef class CDataset:
    cdef dataset *dset
