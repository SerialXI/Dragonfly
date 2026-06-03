#ifndef PARAMS_H
#define PARAMS_H

enum recon_type {RECON3D, RECON2D, RECONRZ} ;

struct params {
	int rank, num_proc ;
	int iteration, current_iter, start_iter, num_iter ;
	char *output_folder, *log_fname ;
	enum recon_type rtype ;
	int save_prob ;
	int verbosity ;
	
	// Algorithm parameters
	int beta_period, data_fraction_period ;
	int need_scaling, known_scale, update_scale, lazy_data ;
	double alpha, beta_jump, beta_factor, beta_config ;
	double data_fraction, coverage_bias ;
	//double *beta, *beta_start ;
	int friedel_sym ; // Symmetrization for 2D recon
	int axial_sym ; // N-fold symmetrization about Z-axis
	int refine, coarse_div, fine_div ; // If doing refinement

	// Radius refinement
	int radius_period ;
	double radius, radius_jump, oversampling ;
	
	// Gaussian EMC parameter
	double sigmasq ;

	// Mode information
	int num_modes, rot_per_mode, nonrot_modes ;
	
	// Testing
	int fixed_seed ;
} ;

#endif //PARAMS_H
