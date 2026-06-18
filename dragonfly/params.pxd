cdef extern from "src/params.h" nogil:
    enum recon_type: RECON3D, RECON2D, RECONRZ
    struct params:
        int rank, num_proc
        int iteration, current_iter, start_iter, num_iter
        char *output_folder
        char *log_fname
        recon_type rtype
        int save_prob
        int verbosity
        
        # Algorithm parameters
        int beta_period, data_fraction_period, data_fraction_until
        int need_scaling, known_scale, update_scale, lazy_data
        double alpha, beta_jump, beta_factor, beta_config
        double data_fraction, coverage_bias
        #double *beta, *beta_start
        int friedel_sym # Symmetrization for 2D recon
        int axial_sym # N-fold symmetrization about Z-axis
        int refine, coarse_div, fine_div # If doing refinement

        # Radius refinement
        int radius_period
        double radius, radius_jump, oversampling
        
        # Gaussian EMC parameter
        double sigmasq

        # Mode information
        int num_modes, rot_per_mode, nonrot_modes

        # 2D beam-center shift search
        int num_shift_x, num_shift_y
        double shift_max_x, shift_max_y

        # Testing
        int fixed_seed

cdef class EMCParams:
    cdef params *par
