#ifndef _DATA_HANDLER_H	
	#define _DATA_HANDLER_H

	#include <hdf5.h>
  #include "mat_utils.h"

	mat_t* arrange_data(double*, hsize_t, hsize_t);
	mat_t* get_data(char*, char*);
	void free_1d_data(double*);
	void free_mat_t(mat_t*);
	void create_new_dset(hid_t, hid_t, char*, void*, char);
	void save_to_h5(double*, int*, int*, long int, char*);
  void save_vec_to_h5(double*, int, char*);
#endif
