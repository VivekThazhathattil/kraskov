#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include "data_handler.h"

/* function to get the h5 data */
/*-------------------------------------------------------------*/
mat_t* get_data(char* file_path, char* h5_dpath)
/*-------------------------------------------------------------*/
{
	hid_t file_id, dset_id, dataspace_id;	
	herr_t status;
	hsize_t dims[2]; // in this case our dspace has dim Npts x Nobs

	file_id = H5Fopen(file_path, H5F_ACC_RDONLY, H5P_DEFAULT);
	dset_id = H5Dopen2(file_id, h5_dpath, H5P_DEFAULT);
	dataspace_id = H5Dget_space(dset_id);
	H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

	double* dset = (double*) malloc(dims[0]*dims[1]*sizeof(double));
	status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		H5P_DEFAULT, dset);
	
	mat_t* res = arrange_data(dset, dims[1], dims[0]);		

	H5Dclose(dset_id);
	H5Sclose(dataspace_id);
	H5Fclose(file_id);

	free_1d_data(dset);

	return res;
}

/*-------------------------------------------------------------*/
mat_t* arrange_data(double* dat, hsize_t m, hsize_t n)
/*-------------------------------------------------------------*/
{
	int i, j;
	mat_t* mat = (mat_t*) malloc(sizeof(mat_t));
	mat->nr = m;
	mat->nc = n;
	mat->m = (double**) malloc(sizeof(double*)*m);
	for(i = 0; i < m; ++i){
		mat->m[i] = (double*) malloc(sizeof(double)*n);
		for(j = 0; j < n; ++j){
			mat->m[i][j] = dat[m*j + i];
		}
	}
	return mat;
}

/*-------------------------------------------------------------*/
void free_1d_data(double* data)
/*-------------------------------------------------------------*/
{
	free(data);
}

/*-------------------------------------------------------------*/
void free_mat_t(mat_t* mat)
/*-------------------------------------------------------------*/
{
	int i;
	for(i = 0; i < mat->nr; ++i)
		free(mat->m[i]);
	free(mat);	
}

/*-------------------------------------------------------------*/
void show_mat(mat_t* mat)
/*-------------------------------------------------------------*/
{
	int i, j;
	for(i = 0; i < mat->nr; ++i){
		for(j = 0; j < mat->nc; ++j)
			printf("%f ", mat->m[i][j]);
		printf("\n");
	}
}

/*-------------------------------------------------------------*/
void create_new_dset(hid_t file_id, hid_t dspace_id, 
	char* dset_name, void* data, char dtype)
/*-------------------------------------------------------------*/
{
	hid_t dset_id, type_id;

	if(dtype == 'd'){
		type_id = H5T_NATIVE_DOUBLE;
		data = (double*) data;
	}

	else if(dtype == 'i'){
		type_id = H5T_NATIVE_INT;
		data = (int*) data;
	}

	dset_id = H5Dcreate2(
							file_id, 			/* location identifier for file */
							dset_name,		/* name of dataset to create */
						  type_id,	 		/* datatype identifier */
						  dspace_id,	 	/* dataspace identifier */
						  H5P_DEFAULT,	/* link creation property list identifer */
						  H5P_DEFAULT,	/* dataset creation property list identifer */
						  H5P_DEFAULT		/* dataset access property list identifer */
					);	
	H5Dwrite(dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	H5Dclose(dset_id);
}

/*-------------------------------------------------------------*/
void save_to_h5(double* val, int* row, int* col, 
	long int num_elms, char* op_file_name)
/*-------------------------------------------------------------*/
{
	hid_t file_id, dspace_id;
	file_id = H5Fcreate(op_file_name, H5F_ACC_TRUNC, H5P_DEFAULT,
		H5P_DEFAULT);
	// H5F_ACC_TRUNC: overwrite existing files
	hsize_t dims[1] = {num_elms};
	dspace_id = H5Screate_simple(1, dims, NULL);
	create_new_dset(file_id, dspace_id, "/row", row, 'i');
	create_new_dset(file_id, dspace_id, "/col", col, 'i');
	create_new_dset(file_id, dspace_id, "/val", val, 'd');
	H5Sclose(dspace_id);
	H5Fclose(file_id);
}


/*-------------------------------------------------------------*/
void save_vec_to_h5(double* vals, int nelms, char* op_file_name)
/*-------------------------------------------------------------*/
{
	hid_t file_id, dspace_id;
	file_id = H5Fcreate(op_file_name, H5F_ACC_TRUNC, H5P_DEFAULT,
		H5P_DEFAULT);
	hsize_t dims[1] = {nelms};
	dspace_id = H5Screate_simple(1, dims, NULL);
	create_new_dset(file_id, dspace_id, "/res", vals, 'd');
	H5Sclose(dspace_id);
	H5Fclose(file_id);
  return;
}

/*-------------------------------------------------------------*/
