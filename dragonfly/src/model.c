#include "model.h"

#define TAU2 0.80901699437 // tau / 2
#define IAU2 0.30901699437 // 1 / (2 * tau)

/* Tri-linear interpolation:
 * Generates slice[t] from model3d[x] by interpolation using given quaternion
 * The locations of the pixels in slice[t] are given by det->qvals[t]
 * The logartihm of the rescaled slice is outputted unless rescale is set to 0.
 */
void slice_gen3d(double *quaternion, int mode, double *slice, struct detector *det, struct model *mod) {
	long t, i, j, x, y, z, size = mod->size ;
	double tx, ty, tz, fx, fy, fz, cx, cy, cz ;
	double rot_pix[3], rot[3][3] = {{0}} ;
	double *model = &mod->model1[mode*mod->vol] ;
	
	make_rot_quat(quaternion, rot) ;
	
	for (t = 0 ; t < det->num_pix ; ++t) {
		for (i = 0 ; i < 3 ; ++i) {
			rot_pix[i] = 0. ;
			for (j = 0 ; j < 3 ; ++j) 
				rot_pix[i] += rot[i][j] * det->qvals[t*3 + j] ;
			rot_pix[i] += mod->center ;
		}
		
		tx = rot_pix[0] ;
		ty = rot_pix[1] ;
		tz = rot_pix[2] ;
		
		x = floor(tx) ;
		y = floor(ty) ;
		z = floor(tz) ;
		
		if (x < 0 || x > size-2 || y < 0 || y > size-2 || z < 0 || z > size-2) {
			slice[t] = DBL_MIN ;
			continue ;
		}
		
		fx = tx - x ;
		fy = ty - y ;
		fz = tz - z ;
		cx = 1. - fx ;
		cy = 1. - fy ;
		cz = 1. - fz ;
		
		slice[t] = cx*cy*cz*model[x*size*size + y*size + z] +
		           cx*cy*fz*model[x*size*size + y*size + ((z+1)%size)] +
		           cx*fy*cz*model[x*size*size + ((y+1)%size)*size + z] +
		           cx*fy*fz*model[x*size*size + ((y+1)%size)*size + ((z+1)%size)] +
		           fx*cy*cz*model[((x+1)%size)*size*size + y*size + z] +
		           fx*cy*fz*model[((x+1)%size)*size*size + y*size + ((z+1)%size)] + 
		           fx*fy*cz*model[((x+1)%size)*size*size + ((y+1)%size)*size + z] + 
		           fx*fy*fz*model[((x+1)%size)*size*size + ((y+1)%size)*size + ((z+1)%size)] ;
		
		// Correct for solid angle and polarization
		slice[t] *= det->corr[t] ;
	}
}

/* Bi-linear interpolation:
 * Generates slice[t] from a stack of 2D models, model[x] using given angle 
 * angle[1] and angle[2] represent the x and y shifts respectively
 * The locations of the pixels in slice[t] are given by det->qvals[t]
 * The logartihm of the rescaled slice is generated unless rescale is set to 0.
 */
void slice_gen2d(double *angle_ptr, int mode, double *slice, struct detector *det, struct model *mod) {
	long t, i, j, x, y, size = mod->size ;
	double tx, ty, fx, fy, cx, cy ;
	double rot_pix[2], rot[2][2] = {{0}} ;
	double angle = angle_ptr[0] ;
	double *model = &mod->model1[mode*mod->vol] ;
	
	make_rot_angle(angle, rot) ;
	
	for (t = 0 ; t < det->num_pix ; ++t) {
		for (i = 0 ; i < 2 ; ++i) {
			rot_pix[i] = 0. ;
			for (j = 0 ; j < 2 ; ++j) 
				rot_pix[i] += rot[i][j] * (det->qvals[t*3 + j] - angle_ptr[j+1]) ;
			rot_pix[i] += mod->center ;
		}
		
		tx = rot_pix[0] ;
		ty = rot_pix[1] ;
		
		x = floor(tx) ;
		y = floor(ty) ;
		
		if (x < 0 || x > size-2 || y < 0 || y > size-2) {
			slice[t] = DBL_MIN ;
			continue ;
		}
		
		fx = tx - x ;
		fy = ty - y ;
		cx = 1. - fx ;
		cy = 1. - fy ;
		
		slice[t] = cx*cy*model[x*size + y] +
		           cx*fy*model[x*size + ((y+1)%size)] +
		           fx*cy*model[((x+1)%size)*size + y] +
		           fx*fy*model[((x+1)%size)*size + ((y+1)%size)] ;
		
		// Correct for solid angle and polarization
		slice[t] *= det->corr[t] ;
		
		if (slice[t] <= 0.)
			slice[t] = DBL_MIN ;
	}
}

/* RZ interpolation
 * Generates slice[t] from model[x] by interpolation at angle phi and beta
 * The locations of the pixels in slice[t] are given by detector[t]
*/
void slice_genrz(double *phibeta, int mode, double *slice, struct detector *det, struct model *mod) {
	int t, x, y, size = mod->size ;
	double tx, ty, fx, fy, cx, cy, fac ;
	double q_beta[3], q_0[3], rot_phi[2][2] = {{0}}, rot_beta[2][2] = {{0}} ;
	double *model = &mod->model1[mode*mod->vol] ;
	
	make_rot_angle(phibeta[0], rot_phi) ;
	make_rot_angle(phibeta[1], rot_beta) ;
	
	for (t = 0 ; t < det->num_pix ; ++t) {
		fac = det->detd * det->detd ;
		fac += pow(det->qvals[t*3+0], 2.) + pow(det->qvals[t*3+1], 2.) ;
		q_0[0] = (det->qvals[t*3+0]*rot_phi[0][0] + det->qvals[t*3+1]*rot_phi[0][1]) / fac ;
		q_0[1] = (det->qvals[t*3+0]*rot_phi[1][0] + det->qvals[t*3+1]*rot_phi[1][1]) / fac ;
		q_0[2] = det->detd/fac - 1. ;
		
		q_beta[0] = q_0[0] ;
		q_beta[1] = q_0[1]*rot_beta[0][0] + q_0[2]*rot_beta[0][1] ;
		q_beta[2] = q_0[1]*rot_beta[1][0] + q_0[2]*rot_beta[1][1] ;
		
		tx = mod->center * (1. + sqrt(q_beta[0]*q_beta[0] + q_beta[2]*q_beta[2])) ;
		ty = mod->center * (1. + q_beta[1]) ;
		
		x = floor(tx) ;
		y = floor(ty) ;
		
		if (x < 0 || x > size-2 || y < 0 || y > size-2) {
			slice[t] = DBL_MIN ;
			continue ;
		}
		
		fx = tx - x ;
		fy = ty - y ;
		cx = 1. - fx ;
		cy = 1. - fy ;
		
		slice[t] = cx*cy*model[x*size+ y] +
		           cx*fy*model[x*size+ ((y+1)%size)] +
		           fx*cy*model[((x+1)%size)*size + y] +
		           fx*fy*model[((x+1)%size)*size + ((y+1)%size)] ;
		
		// Correct for solid angle and polarization
		slice[t] *= det->corr[t] ;
		
		if (slice[t] <= 0.)
			slice[t] = DBL_MIN ;
	}
}

/* Tri-linear merging:
 * Merges slice[t] into model3d[x] using given quaternion
 * Also adds to weight[x] containing the interpolation weights
 * The locations of the pixels in slice[t] are given by det->qvals[t]
 * Only pixels with a mask value < 2 are merged
 */
void slice_merge3d(double *quaternion, int mode, double *slice, double *model2, double *inter_weight, long size, struct detector *det) {
	long t, i, j, x, y, z, center = size / 2, vol = size*size*size ;
	double tx, ty, tz, fx, fy, fz, cx, cy, cz, w, f ;
	double rot_pix[3], rot[3][3] = {{0}} ;
	double *model = &model2[mode*vol] ;
	double *weight = &inter_weight[mode*vol] ;
	
	make_rot_quat(quaternion, rot) ;
	
	for (t = 0 ; t < det->num_pix ; ++t) {
		if (det->raw_mask[t] > 1)
			continue ;
		
		for (i = 0 ; i < 3 ; ++i) {
			rot_pix[i] = 0. ;
			for (j = 0 ; j < 3 ; ++j)
				rot_pix[i] += rot[i][j] * det->qvals[t*3 + j] ;
			rot_pix[i] += center ;
		}
		
		tx = rot_pix[0] ;
		ty = rot_pix[1] ;
		tz = rot_pix[2] ;
		
		x = floor(tx) ;
		y = floor(ty) ;
		z = floor(tz) ;
		
		if (x < 0 || x > size-2 || y < 0 || y > size-2 || z < 0 || z > size-2)
			continue ;
		
		fx = tx - x ;
		fy = ty - y ;
		fz = tz - z ;
		cx = 1. - fx ;
		cy = 1. - fy ;
		cz = 1. - fz ;
		
		// Correct for solid angle and polarization
		w = slice[t] / det->corr[t] ;
		
		f = cx*cy*cz ;
		weight[x*size*size + y*size + z] += f ;
		model[x*size*size + y*size + z] += f * w ;
		
		f = cx*cy*fz ;
		weight[x*size*size + y*size + ((z+1)%size)] += f ;
		model[x*size*size + y*size + ((z+1)%size)] += f * w ;
		
		f = cx*fy*cz ;
		weight[x*size*size + ((y+1)%size)*size + z] += f ;
		model[x*size*size + ((y+1)%size)*size + z] += f * w ;
		
		f = cx*fy*fz ;
		weight[x*size*size + ((y+1)%size)*size + ((z+1)%size)] += f ;
		model[x*size*size + ((y+1)%size)*size + ((z+1)%size)] += f * w ;
		
		f = fx*cy*cz ;
		weight[((x+1)%size)*size*size + y*size + z] += f ;
		model[((x+1)%size)*size*size + y*size + z] += f * w ;
		
		f = fx*cy*fz ;
		weight[((x+1)%size)*size*size + y*size + ((z+1)%size)] += f ;
		model[((x+1)%size)*size*size + y*size + ((z+1)%size)] += f * w ;
		
		f = fx*fy*cz ;
		weight[((x+1)%size)*size*size + ((y+1)%size)*size + z] += f ;
		model[((x+1)%size)*size*size + ((y+1)%size)*size + z] += f * w ;
		
		f = fx*fy*fz ;
		weight[((x+1)%size)*size*size + ((y+1)%size)*size + ((z+1)%size)] += f ;
		model[((x+1)%size)*size*size + ((y+1)%size)*size + ((z+1)%size)] += f * w ;
	}
}

/* Bi-linear merging:
 * Merges slice[t] into a stack of 2D models, model[x] using given angle
 * angle[1] and angle[2] represent the x and y shift respectively
 * Also adds to weight[x] containing the interpolation weights
 * The locations of the pixels in slice[t] are given by det->qvals[t]
 * Only pixels with a mask value < 2 are merged
 */
void slice_merge2d(double *angle_ptr, int mode, double *slice, double *model2, double *inter_weight, long size, struct detector *det) {
	long t, i, j, x, y, center = size/2, vol = size*size ;
	double tx, ty, fx, fy, cx, cy, w, f ;
	double rot_pix[2], rot[2][2] = {{0}} ;
	double angle = angle_ptr[0] ;
	double *model = &model2[mode*vol] ;
	double *weight = &inter_weight[mode*vol] ;
	
	make_rot_angle(angle, rot) ;
	
	for (t = 0 ; t < det->num_pix ; ++t) {
		if (det->raw_mask[t] > 1)
			continue ;
		
		for (i = 0 ; i < 2 ; ++i) {
			rot_pix[i] = 0. ;
			for (j = 0 ; j < 2 ; ++j)
				rot_pix[i] += rot[i][j] * (det->qvals[t*3 + j] - angle_ptr[j+1]) ;
			rot_pix[i] += center ;
		}
		
		tx = rot_pix[0] ;
		ty = rot_pix[1] ;
		
		x = floor(tx) ;
		y = floor(ty) ;
		
		if (x < 0 || x > size-2 || y < 0 || y > size-2)
			continue ;
		
		fx = tx - x ;
		fy = ty - y ;
		cx = 1. - fx ;
		cy = 1. - fy ;
		
		// Correct for solid angle and polarization
		w = slice[t] / det->corr[t] ;
		
		f = cx*cy ;
		weight[x*size + y] += f ;
		model[x*size + y] += f * w ;
		
		f = cx*fy ;
		weight[x*size + ((y+1)%size)] += f ;
		model[x*size + ((y+1)%size)] += f * w ;
		
		f = fx*cy ;
		weight[((x+1)%size)*size + y] += f ;
		model[((x+1)%size)*size + y] += f * w ;
		
		f = fx*fy ;
		weight[((x+1)%size)*size + ((y+1)%size)] += f ;
		model[((x+1)%size)*size + ((y+1)%size)] += f * w ;
	}
}

/* RZ merging
 * Merges slice[t] into model[x] at the given phi, beta angles
 * Also adds to weight[x] containing the interpolation weights
 * The locations of the pixels in slice[t] are given by detector[t]
*/
void slice_mergerz(double *phibeta, int mode, double *slice, double *model2, double *inter_weight, long size, struct detector *det) {
	int t, x, y, center = size/2, vol = size*size ;
	double tx, ty, fx, fy, cx, cy, w, f, fac ;
	double q_beta[3], q_0[3], rot_phi[2][2] = {{0}}, rot_beta[2][2] = {{0}} ;
	double *model = &model2[mode*vol] ;
	double *weight = &inter_weight[mode*vol] ;
	
	make_rot_angle(phibeta[0], rot_phi) ;
	make_rot_angle(phibeta[1], rot_beta) ;
	
	for (t = 0 ; t < det->num_pix ; ++t) {
		if (det->raw_mask[t] > 1)
			continue ;
		
		fac = det->detd * det->detd ;
		fac += pow(det->qvals[t*3+0], 2.) + pow(det->qvals[t*3+1], 2.) ;
		q_0[0] = (det->qvals[t*3+0]*rot_phi[0][0] + det->qvals[t*3+1]*rot_phi[0][1]) / fac ;
		q_0[1] = (det->qvals[t*3+0]*rot_phi[1][0] + det->qvals[t*3+1]*rot_phi[1][1]) / fac ;
		q_0[2] = det->detd/fac - 1. ;
		
		q_beta[0] = q_0[0] ;
		q_beta[1] = q_0[1]*rot_beta[0][0] + q_0[2]*rot_beta[0][1] ;
		q_beta[2] = q_0[1]*rot_beta[1][0] + q_0[2]*rot_beta[1][1] ;
		
		tx = center * (1. + sqrt(q_beta[0]*q_beta[0] + q_beta[2]*q_beta[2])) ;
		ty = center * (1. + q_beta[1]) ;
		
		x = floor(tx) ;
		y = floor(ty) ;
		
		if (x < 0 || x > size-2 || y < 0 || y > size-2)
			continue ;
		
		fx = tx - x ;
		fy = ty - y ;
		cx = 1. - fx ;
		cy = 1. - fy ;
		
		// Correct for solid angle and polarization
		w = slice[t] / det->corr[t] ;
		
		f = cx*cy ;
		weight[x*size + y] += f ;
		model[x*size + y] += f * w ;
		
		f = cx*fy ;
		weight[x*size + ((y+1)%size)] += f ;
		model[x*size + ((y+1)%size)] += f * w ;
		
		f = fx*cy ;
		weight[((x+1)%size)*size + y] += f ;
		model[((x+1)%size)*size + y] += f * w ;
		
		f = fx*fy ;
		weight[((x+1)%size)*size + ((y+1)%size)] += f ;
		model[((x+1)%size)*size + ((y+1)%size)] += f * w ;
	}
}

/* Rotates cubic model according to given rotation matrix
 * 	Adds to rotated model. Does not zero output model.
 * 	Note that this function uses OpenMP so it should not be put in a parallel block
 * 	Arguments:
 * 		rot[3][3] - Rotation matrix
 * 		m - Pointer to model to rotate
 * 		s - Size of model. Center assumed to be at (s/2, s/2, s/2)
 * 		rotmodel - Pointer to rotated model
 */
void rotate_model(double rot[3][3], double *m, int s, int max_r, double *rotmodel) {
	int x, y, z, i, c = s/2, vx, vy, vz ;
	double fx, fy, fz, cx, cy, cz ;
	double rot_vox[3] ;
	if (max_r == 0)
		max_r = c ;
	
	for (vx = -max_r ; vx < max_r ; ++vx)
	for (vy = -max_r ; vy < max_r ; ++vy)
	for (vz = -max_r ; vz < max_r ; ++vz) {
		for (i = 0 ; i < 3 ; ++i) {
			rot_vox[i] = 0. ;
			rot_vox[i] += rot[i][0]*vx + rot[i][1]*vy + rot[i][2]*vz ;
			rot_vox[i] += c ;
		}
		
		if (rot_vox[0] < 0 || rot_vox[0] >= s - 1) continue ;
		if (rot_vox[1] < 0 || rot_vox[1] >= s - 1) continue ;
		if (rot_vox[2] < 0 || rot_vox[2] >= s - 1) continue ;
		
		x = (int) rot_vox[0] ;
		y = (int) rot_vox[1] ;
		z = (int) rot_vox[2] ;
		fx = rot_vox[0] - x ;
		fy = rot_vox[1] - y ;
		fz = rot_vox[2] - z ;
		cx = 1. - fx ;
		cy = 1. - fy ;
		cz = 1. - fz ;
		
		rotmodel[(vx+c)*s*s + (vy+c)*s + (vz+c)] +=
			cx*cy*cz*m[x*s*s + y*s + z] + 
			cx*cy*fz*m[x*s*s + y*s + ((z+1)%s)] + 
			cx*fy*cz*m[x*s*s + ((y+1)%s)*s + z] + 
			cx*fy*fz*m[x*s*s + ((y+1)%s)*s + ((z+1)%s)] + 
			fx*cy*cz*m[((x+1)%s)*s*s + y*s + z] +
			fx*cy*fz*m[((x+1)%s)*s*s + y*s + ((z+1)%s)] +
			fx*fy*cz*m[((x+1)%s)*s*s + ((y+1)%s)*s + z] +
			fx*fy*fz*m[((x+1)%s)*s*s + ((y+1)%s)*s + ((z+1)%s)] ;
	}
}

void rotate_model_openmp(double rot[3][3], double *m, int s, double *rotmodel) {
	#pragma omp parallel default(shared)
	{
		int x, y, z, i, c = s/2, vx, vy, vz ;
		double fx, fy, fz, cx, cy, cz ;
		double rot_vox[3] ;
		
		#pragma omp for schedule(static,1)
		for (vx = -c ; vx < s-c-1 ; ++vx)
		for (vy = -c ; vy < s-c-1 ; ++vy)
		for (vz = -c ; vz < s-c-1 ; ++vz) {
			for (i = 0 ; i < 3 ; ++i) {
				rot_vox[i] = 0. ;
				rot_vox[i] += rot[i][0]*vx + rot[i][1]*vy + rot[i][2]*vz ;
				rot_vox[i] += c ;
			}
			
			if (rot_vox[0] < 0 || rot_vox[0] >= s - 1) continue ;
			if (rot_vox[1] < 0 || rot_vox[1] >= s - 1) continue ;
			if (rot_vox[2] < 0 || rot_vox[2] >= s - 1) continue ;
			
			x = (int) rot_vox[0] ;
			y = (int) rot_vox[1] ;
			z = (int) rot_vox[2] ;
			fx = rot_vox[0] - x ;
			fy = rot_vox[1] - y ;
			fz = rot_vox[2] - z ;
			cx = 1. - fx ;
			cy = 1. - fy ;
			cz = 1. - fz ;
			
			rotmodel[(vx+c)*s*s + (vy+c)*s + (vz+c)] +=
				cx*cy*cz*m[x*s*s + y*s + z] + 
				cx*cy*fz*m[x*s*s + y*s + ((z+1)%s)] + 
				cx*fy*cz*m[x*s*s + ((y+1)%s)*s + z] + 
				cx*fy*fz*m[x*s*s + ((y+1)%s)*s + ((z+1)%s)] + 
				fx*cy*cz*m[((x+1)%s)*s*s + y*s + z] +
				fx*cy*fz*m[((x+1)%s)*s*s + y*s + ((z+1)%s)] +
				fx*fy*cz*m[((x+1)%s)*s*s + ((y+1)%s)*s + z] +
				fx*fy*fz*m[((x+1)%s)*s*s + ((y+1)%s)*s + ((z+1)%s)] ;
		}
	}
}

/* In-place Friedel symmetrization:
 * I(q) and I(-q) are replaced by their average
 */
void symmetrize_friedel(double *array, double *weights, int size) {
	int x, y, z, vox1, vox2, min = 0, center = size/2 ;
	double avg_intens, avg_weights ;
	if (size % 2 == 0)
		min = 1 ;
	
	for (x = min ; x < size ; ++x)
	for (y = min ; y < size ; ++y)
	for (z = min ; z <= center ; ++z) {
		vox1 = x*size*size + y*size + z ;
		vox2 = (2*center-x)*size*size + (2*center-y)*size + (2*center-z) ;
		avg_intens = 0.5 * (array[vox1]*weights[vox1] + array[vox2]*weights[vox2]) ;
		avg_weights = 0.5 * (weights[vox1] + weights[vox2]) ;
		if (avg_weights > 0)
			avg_intens /= avg_weights ;
		else
			avg_intens = 0. ;
		array[vox1] = avg_intens ;
		array[vox2] = avg_intens ;
	}
}

/* In-place Friedel symmetrization for 2D stack:
 * I(q) and I(-q) are replaced by their average for each 2D layer independently
 */
void symmetrize_friedel2d(double *array, double *weights, int num_layers, int size) {
	int x, y, z, vox1, vox2, min = 0, center = size/2 ;
	double avg_intens, avg_weights ;
	if (size % 2 == 0)
		min = 1 ;
	
	for (x = 0; x < num_layers; ++x)
	for (y = min ; y < size ; ++y)
	for (z = min ; z <= center ; ++z) {
		vox1 = x*size*size + y*size + z ;
		vox2 = x*size*size + (2*center-y)*size + (2*center-z) ;
		avg_intens = .5 * (array[vox1]*weights[vox1] + array[vox2]*weights[vox2]) ;
		avg_weights = .5 * (weights[vox1] + weights[vox2]) ;
		if (avg_weights > 0)
			avg_intens /= avg_weights ;
		else
			avg_intens = 0. ;
		array[vox1] = avg_intens ;
		array[vox2] = avg_intens ;
	}
}

/* Axial symmetrization
 * 	N-fold axial symmetrization about Z-axis
 * 	Arguments:
 * 		Pointer to model representing centered 3D volume
 * 		Pointer to weights
 * 		Size of model
 * 		Symmetry order
 * 	In main EMC code, Friedel symmetry is independently applied
 * 	No return value. Symmetrization performed in-place
 */
void symmetrize_axial(double *model, double *weights, int size, int order) {
	long i, vol = (long)size*size*size ;
	double *temp = malloc(vol*sizeof(double)) ;
	double angle, c, s, rot[3][3] = {{1,0,0},{0,1,0},{0,0,1}} ;

	// Calculate numerator: model <- SYM(model * weights)
	memcpy(temp, model, vol*sizeof(double)) ;
	for (i = 0 ; i < vol ; ++i)
		temp[i] *= weights[i] ;
	memset(model, 0, vol*sizeof(double)) ;
	for (i = 0 ; i < order ; ++i) {
		angle = 2. * M_PI * i / order ;
		c = cos(angle) ;
		s = sin(angle) ;
		rot[0][0] = c ;
		rot[0][1] = -s ;
		rot[1][1] = c ;
		rot[1][0] = s ;
		rotate_model_openmp(rot, temp, size, model) ;
	}

	// Calculate denominator: weights <- SYM(weights)
	memcpy(temp, weights, vol*sizeof(double)) ;
	memset(weights, 0, vol*sizeof(double)) ;
	for (i = 0 ; i < order ; ++i) {
		angle = 2. * M_PI * i / order ;
		c = cos(angle) ;
		s = sin(angle) ;
		rot[0][0] = c ;
		rot[0][1] = -s ;
		rot[1][1] = c ;
		rot[1][0] = s ;
		rotate_model_openmp(rot, temp, size, weights) ;
	}

	// Divide numerator by denominator
	for (i = 0 ; i < vol ; ++i)
	if (weights[i] > 0.)
		model[i] /= weights[i] ;
	else
		model[i] = 0. ;

	free(temp) ;
}


static double icos_list[60][3][3] = {
	{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}},
	{{-1., 0., 0.}, {0., -1., 0.}, {0., 0., 1.}},
	{{0., 0., 1.}, {1., 0., 0.}, {0., 1., 0.}},
	{{0., 1., 0.}, {0., 0., 1.}, {1., 0., 0.}},
	{{0., -1., 0.}, {0., 0., -1.}, {1., 0., 0.}},
	{{0., 0., 1.}, {-1., 0., 0.}, {0., -1., 0.}},
	{{0., 1., 0.}, {0., 0., -1.}, {-1., 0., 0.}},
	{{0., 0., -1.}, {-1., 0., 0.}, {0., 1., 0.}},
	{{1., 0., 0.}, {0., -1., 0.}, {0., 0., -1.}},
	{{0., 0., -1.}, {1., 0., 0.}, {0., -1., 0.}},
	{{0., -1., 0.}, {0., 0., 1.}, {-1., 0., 0.}},
	{{-1., 0., 0.}, {0., 1., 0.}, {0., 0., -1.}},
	{{0.5, -TAU2, IAU2}, {TAU2, IAU2, -0.5}, {IAU2, 0.5, TAU2}},
	{{-0.5, TAU2, IAU2}, {-TAU2, -IAU2, -0.5}, {-IAU2, -0.5, TAU2}},
	{{-TAU2, IAU2, 0.5}, {IAU2, -0.5, TAU2}, {0.5, TAU2, IAU2}},
	{{IAU2, 0.5, -TAU2}, {-0.5, TAU2, IAU2}, {TAU2, IAU2, 0.5}},
	{{IAU2, -0.5, TAU2}, {-0.5, -TAU2, -IAU2}, {TAU2, -IAU2, -0.5}},
	{{TAU2, -IAU2, 0.5}, {-IAU2, 0.5, TAU2}, {-0.5, -TAU2, IAU2}},
	{{-IAU2, 0.5, TAU2}, {0.5, TAU2, -IAU2}, {-TAU2, IAU2, -0.5}},
	{{TAU2, IAU2, -0.5}, {-IAU2, -0.5, -TAU2}, {-0.5, TAU2, -IAU2}},
	{{0.5, TAU2, -IAU2}, {TAU2, -IAU2, 0.5}, {IAU2, -0.5, -TAU2}},
	{{-TAU2, -IAU2, -0.5}, {IAU2, 0.5, -TAU2}, {0.5, -TAU2, -IAU2}},
	{{-IAU2, -0.5, -TAU2}, {0.5, -TAU2, IAU2}, {-TAU2, -IAU2, 0.5}},
	{{-0.5, -TAU2, -IAU2}, {-TAU2, IAU2, 0.5}, {-IAU2, 0.5, -TAU2}},
	{{-IAU2, -0.5, TAU2}, {0.5, -TAU2, -IAU2}, {TAU2, IAU2, 0.5}},
	{{IAU2, 0.5, TAU2}, {-0.5, TAU2, -IAU2}, {-TAU2, -IAU2, 0.5}},
	{{-0.5, TAU2, -IAU2}, {-TAU2, -IAU2, 0.5}, {IAU2, 0.5, TAU2}},
	{{TAU2, -IAU2, -0.5}, {-IAU2, 0.5, -TAU2}, {0.5, TAU2, IAU2}},
	{{TAU2, IAU2, 0.5}, {-IAU2, -0.5, TAU2}, {0.5, -TAU2, -IAU2}},
	{{0.5, -TAU2, -IAU2}, {TAU2, IAU2, 0.5}, {-IAU2, -0.5, TAU2}},
	{{-TAU2, -IAU2, 0.5}, {IAU2, 0.5, TAU2}, {-0.5, TAU2, -IAU2}},
	{{0.5, TAU2, IAU2}, {TAU2, -IAU2, -0.5}, {-IAU2, 0.5, -TAU2}},
	{{-IAU2, 0.5, -TAU2}, {0.5, TAU2, IAU2}, {TAU2, -IAU2, -0.5}},
	{{-0.5, -TAU2, IAU2}, {-TAU2, IAU2, -0.5}, {IAU2, -0.5, -TAU2}},
	{{-TAU2, IAU2, -0.5}, {IAU2, -0.5, -TAU2}, {-0.5, -TAU2, IAU2}},
	{{IAU2, -0.5, -TAU2}, {-0.5, -TAU2, IAU2}, {-TAU2, IAU2, -0.5}},
	{{-IAU2, 0.5, TAU2}, {-0.5, -TAU2, IAU2}, {TAU2, -IAU2, 0.5}},
	{{IAU2, -0.5, TAU2}, {0.5, TAU2, IAU2}, {-TAU2, IAU2, 0.5}},
	{{0.5, TAU2, -IAU2}, {-TAU2, IAU2, -0.5}, {-IAU2, 0.5, TAU2}},
	{{TAU2, -IAU2, 0.5}, {IAU2, -0.5, -TAU2}, {0.5, TAU2, -IAU2}},
	{{TAU2, IAU2, -0.5}, {IAU2, 0.5, TAU2}, {0.5, -TAU2, IAU2}},
	{{-0.5, -TAU2, -IAU2}, {TAU2, -IAU2, -0.5}, {IAU2, -0.5,   TAU2}},
	{{-TAU2, -IAU2, -0.5}, {-IAU2, -0.5,   TAU2}, {-0.5, TAU2, IAU2}},
	{{-0.5, TAU2,   IAU2}, {TAU2, IAU2, 0.5}, {IAU2,   0.5, -TAU2}},
	{{-IAU2, -0.5, -TAU2}, {-0.5,   TAU2, -IAU2}, {TAU2, IAU2, -0.5}},
	{{0.5, -TAU2,   IAU2}, {-TAU2, -IAU2,   0.5}, {-IAU2, -0.5, -TAU2}},
	{{-TAU2, IAU2,   0.5}, {-IAU2,   0.5, -TAU2}, {-0.5, -TAU2, -IAU2}},
	{{IAU2,   0.5, -TAU2}, {0.5, -TAU2, -IAU2}, {-TAU2, -IAU2, -0.5}},
	{{0.5, TAU2, IAU2}, {-TAU2, IAU2,   0.5}, {IAU2, -0.5, TAU2}},
	{{-0.5, -TAU2,   IAU2}, {TAU2, -IAU2, 0.5}, {-IAU2, 0.5,   TAU2}},
	{{TAU2, IAU2, 0.5}, {IAU2,   0.5, -TAU2}, {-0.5, TAU2, IAU2}},
	{{IAU2, 0.5,   TAU2}, {0.5, -TAU2, IAU2}, {TAU2,   IAU2, -0.5}},
	{{IAU2, -0.5, -TAU2}, {0.5,   TAU2, -IAU2}, {TAU2, -IAU2,   0.5}},
	{{-TAU2, -IAU2, 0.5}, {-IAU2, -0.5, -TAU2}, {0.5, -TAU2, IAU2}},
	{{-IAU2, 0.5, -TAU2}, {-0.5, -TAU2, -IAU2}, {-TAU2, IAU2, 0.5}},
	{{-TAU2, IAU2, -0.5}, {-IAU2, 0.5, TAU2}, {0.5, TAU2, -IAU2}},
	{{0.5, -TAU2, -IAU2}, {-TAU2, -IAU2, -0.5}, {IAU2, 0.5, -TAU2}},
	{{TAU2, -IAU2, -0.5}, {IAU2, -0.5, TAU2}, {-0.5, -TAU2, -IAU2}},
	{{-IAU2, -0.5, TAU2}, {-0.5, TAU2, IAU2}, {-TAU2, -IAU2, -0.5}},
	{{-0.5, TAU2, -IAU2}, {TAU2, IAU2, -0.5}, {-IAU2, -0.5, -TAU2}}
} ;

static double cube_list[48][3][3] = {
	{{ 1,  0,  0}, { 0,  1,  0}, { 0,  0,  1}},
	{{ 1,  0,  0}, { 0,  1,  0}, { 0,  0, -1}},
	{{ 1,  0,  0}, { 0, -1,  0}, { 0,  0,  1}},
	{{ 1,  0,  0}, { 0, -1,  0}, { 0,  0, -1}},
	{{-1,  0,  0}, { 0,  1,  0}, { 0,  0,  1}},
	{{-1,  0,  0}, { 0,  1,  0}, { 0,  0, -1}},
	{{-1,  0,  0}, { 0, -1,  0}, { 0,  0,  1}},
	{{-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1}},
	{{ 0,  1,  0}, { 1,  0,  0}, { 0,  0,  1}},
	{{ 0,  1,  0}, { 1,  0,  0}, { 0,  0, -1}},
	{{ 0,  1,  0}, {-1,  0,  0}, { 0,  0,  1}},
	{{ 0,  1,  0}, {-1,  0,  0}, { 0,  0, -1}},
	{{ 0, -1,  0}, { 1,  0,  0}, { 0,  0,  1}},
	{{ 0, -1,  0}, { 1,  0,  0}, { 0,  0, -1}},
	{{ 0, -1,  0}, {-1,  0,  0}, { 0,  0,  1}},
	{{ 0, -1,  0}, {-1,  0,  0}, { 0,  0, -1}},
	{{ 0,  0,  1}, { 0,  1,  0}, { 1,  0,  0}},
	{{ 0,  0,  1}, { 0,  1,  0}, {-1,  0,  0}},
	{{ 0,  0,  1}, { 0, -1,  0}, { 1,  0,  0}},
	{{ 0,  0,  1}, { 0, -1,  0}, {-1,  0,  0}},
	{{ 0,  0, -1}, { 0,  1,  0}, { 1,  0,  0}},
	{{ 0,  0, -1}, { 0,  1,  0}, {-1,  0,  0}},
	{{ 0,  0, -1}, { 0, -1,  0}, { 1,  0,  0}},
	{{ 0,  0, -1}, { 0, -1,  0}, {-1,  0,  0}},
	{{ 1,  0,  0}, { 0,  0,  1}, { 0,  1,  0}},
	{{ 1,  0,  0}, { 0,  0,  1}, { 0, -1,  0}},
	{{ 1,  0,  0}, { 0,  0, -1}, { 0,  1,  0}},
	{{ 1,  0,  0}, { 0,  0, -1}, { 0, -1,  0}},
	{{-1,  0,  0}, { 0,  0,  1}, { 0,  1,  0}},
	{{-1,  0,  0}, { 0,  0,  1}, { 0, -1,  0}},
	{{-1,  0,  0}, { 0,  0, -1}, { 0,  1,  0}},
	{{-1,  0,  0}, { 0,  0, -1}, { 0, -1,  0}},
	{{ 0,  1,  0}, { 0,  0,  1}, { 1,  0,  0}},
	{{ 0,  1,  0}, { 0,  0,  1}, {-1,  0,  0}},
	{{ 0,  1,  0}, { 0,  0, -1}, { 1,  0,  0}},
	{{ 0,  1,  0}, { 0,  0, -1}, {-1,  0,  0}},
	{{ 0, -1,  0}, { 0,  0,  1}, { 1,  0,  0}},
	{{ 0, -1,  0}, { 0,  0,  1}, {-1,  0,  0}},
	{{ 0, -1,  0}, { 0,  0, -1}, { 1,  0,  0}},
	{{ 0, -1,  0}, { 0,  0, -1}, {-1,  0,  0}},
	{{ 0,  0,  1}, { 1,  0,  0}, { 0,  1,  0}},
	{{ 0,  0,  1}, { 1,  0,  0}, { 0, -1,  0}},
	{{ 0,  0,  1}, {-1,  0,  0}, { 0,  1,  0}},
	{{ 0,  0,  1}, {-1,  0,  0}, { 0, -1,  0}},
	{{ 0,  0, -1}, { 1,  0,  0}, { 0,  1,  0}},
	{{ 0,  0, -1}, { 1,  0,  0}, { 0, -1,  0}},
	{{ 0,  0, -1}, {-1,  0,  0}, { 0,  1,  0}},
	{{ 0,  0, -1}, {-1,  0,  0}, { 0, -1,  0}},
} ;

/* Octahedral symmetrization
 * 	Assumes vertices are at permutations of (+-1, +-1, +-1)
 * 	Arguments:
 * 		Pointer to model representing centered 3D volume
 * 		Pointer to weights
 * 		Size of model
 * 	No return value. Symmetrization performed in-place
 */
void symmetrize_octahedral(double *model, double *weights, int size) {
	long i, vol = (long)size*size*size ;
	double *temp = malloc(vol*sizeof(double)) ;
	
	// Calculate numerator: model <- SYM(model * weights)
	memcpy(temp, model, vol*sizeof(double)) ;
	for (i = 0 ; i < vol ; ++i)
		temp[i] *= weights[i] ;
	memset(model, 0, vol*sizeof(double)) ;
	for (i = 0 ; i < 48 ; ++i)
		rotate_model_openmp(cube_list[i], temp, size, model) ;
	
	// Calculate denominator: weights <- SYM(weights)
	memcpy(temp, weights, vol*sizeof(double)) ;
	memset(weights, 0, vol*sizeof(double)) ;
	for (i = 0 ; i < 48 ; ++i)
		rotate_model_openmp(cube_list[i], temp, size, weights) ;
	
	// Divide numerator by denominator
	for (i = 0 ; i < vol ; ++i)
	if (weights[i] > 0.)
		model[i] /= weights[i] ;
	else
		model[i] = 0. ;
	
	free(temp) ;
}

/* Icosahedral symmetrization
 * 	Assumes vertices are at permutations of (0, +-1, +-tau)
 * 	Arguments:
 * 		Pointer to model representing centered 3D volume
 * 		Pointer to weights
 * 		Size of model
 * 	No return value. Symmetrization performed in-place
 */
void symmetrize_icosahedral(double *model, double *weights, int size) {
	long i, vol = (long)size*size*size ;
	double *temp = malloc(vol*sizeof(double)) ;
	
	// Calculate numerator: model <- SYM(model * weights)
	memcpy(temp, model, vol*sizeof(double)) ;
	for (i = 0 ; i < vol ; ++i)
		temp[i] *= weights[i] ;
	memset(model, 0, vol*sizeof(double)) ;
	for (i = 0 ; i < 60 ; ++i)
		rotate_model_openmp(icos_list[i], temp, size, model) ;
	
	// Calculate denominator: weights <- SYM(weights)
	memcpy(temp, weights, vol*sizeof(double)) ;
	memset(weights, 0, vol*sizeof(double)) ;
	for (i = 0 ; i < 60 ; ++i)
		rotate_model_openmp(icos_list[i], temp, size, weights) ;
	
	// Divide numerator by denominator
	for (i = 0 ; i < vol ; ++i)
	if (weights[i] > 0.)
		model[i] /= weights[i] ;
	else
		model[i] = 0. ;
	
	free(temp) ;
}

void make_rot_quat(double *quaternion, double rot[3][3]) {
	double q0, q1, q2, q3, q01, q02, q03, q11, q12, q13, q22, q23, q33 ;
	
	q0 = quaternion[0] ;
	q1 = quaternion[1] ;
	q2 = quaternion[2] ;
	q3 = quaternion[3] ;
	
	q01 = q0*q1 ;
	q02 = q0*q2 ;
	q03 = q0*q3 ;
	q11 = q1*q1 ;
	q12 = q1*q2 ;
	q13 = q1*q3 ;
	q22 = q2*q2 ;
	q23 = q2*q3 ;
	q33 = q3*q3 ;
	
	rot[0][0] = (1. - 2.*(q22 + q33)) ;
	rot[0][1] = 2.*(q12 + q03) ;
	rot[0][2] = 2.*(q13 - q02) ;
	rot[1][0] = 2.*(q12 - q03) ;
	rot[1][1] = (1. - 2.*(q11 + q33)) ;
	rot[1][2] = 2.*(q01 + q23) ;
	rot[2][0] = 2.*(q02 + q13) ;
	rot[2][1] = 2.*(q23 - q01) ;
	rot[2][2] = (1. - 2.*(q11 + q22)) ;
}

void make_rot_angle(double angle, double rot[2][2]) {
	double c = cos(angle), s = sin(angle) ;
	
	rot[0][0] = c ;
	rot[0][1] = -s ;
	rot[1][0] = s ;
	rot[1][1] = c ;
}

