#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include "BrainDist.h"

static PyObject*  BrainDist_wm_dist(PyObject* self,  PyObject* args) {
	char* mesh_fn, *matrix_fn, *surface_mask_fn, *example_fn, *density_fn;
	long int label;
    int	wm_search_depth, max_threads, write_vertex, subsample, subsample_factor, VERBOSE;
	PyObject* img_vol_obj, *starts_obj, *steps_obj, *maxima_obj;

	PyArg_ParseTuple(args, "OOOOsssssiiiiiii", &img_vol_obj, &steps_obj, &maxima_obj, &starts_obj,  &mesh_fn,  &matrix_fn, &surface_mask_fn, &example_fn, &density_fn, &label, &wm_search_depth, &max_threads, &write_vertex, &subsample, &subsample_factor, &VERBOSE );
	PyObject *img_vol_array = PyArray_FROM_OTF(img_vol_obj, NPY_INT64, NPY_IN_ARRAY);
	PyObject *starts_array = PyArray_FROM_OTF(starts_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *steps_array = PyArray_FROM_OTF(steps_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *maxima_array = PyArray_FROM_OTF(maxima_obj, NPY_INT64, NPY_IN_ARRAY);
	
	long int *img_vol = (long int*) PyArray_DATA(img_vol_array);
	double *starts = (double*) PyArray_DATA(starts_array);
	double *steps = (double*) PyArray_DATA(steps_array);
	long int *maxima = (int*) PyArray_DATA(maxima_array);

	wm_dist( img_vol,  mesh_fn,  matrix_fn, surface_mask_fn, example_fn, density_fn, label, wm_search_depth, max_threads,  write_vertex,  subsample,  subsample_factor, steps, maxima, starts, VERBOSE);
	Py_RETURN_NONE;
}

/*static PyObject*  BrainDist_surf_dist(PyObject* self,  PyObject* args) {
	char *mask_fn, *output_fn;
	double dt;
    long int	inside_label, subsample, subsample_factor, n, VERBOSE;
	PyObject* coords_array, *nngh_array, *ngh_array;

	//opts.inside_label, subsample,subsample_factor, np.int64(n), np.int64(opts.verbose), opts.mask_fn, opts.output_fn,  coords, nngh, ngh
	PyArg_ParseTuple(args, "iiiiissOOO",&inside_label, &subsample, &subsample_factor, &n,&VERBOSE, &mask_fn, &output_fn,   &coords_array, &nngh_array, &ngh_array  );

	double *coords_0 = (double*) PyArray_DATA(coords_array), **coords;
	long int *ngh_0 = (int**) PyArray_DATA(ngh_array), **ngh;
	long int *nngh = (long int*) PyArray_DATA(nngh_array);
	

	printf("label %d\n", inside_label);
	printf("subsample: %d\n", subsample);
	printf("subsample_factor: %d\n", subsample_factor);
	printf("Mask: %s\n", mask_fn);
	printf("Output: %s\n", output_fn);
	printf("n: %d\n", n);
	printf("verbose %d\n", VERBOSE);

	int index=0;

	ngh = malloc(sizeof(*ngh) * n );
	coords = malloc(sizeof(*coords) * n );
	for(int i=0; i<n; i++){
	  	//prin("%d (%d) ", i, nngh[i]); 
		ngh[i] = malloc( sizeof(**ngh) * nngh[i] );
		coords[i] = malloc( sizeof(**coords) * 3 );
		for(int j=0; j < nngh[i]; j++){
			ngh[i][j] = ngh_0[index];
			index++;
		}
	   coords[i][0] = coords_0[ i * 3 + 0  ];	
	   coords[i][1] = coords_0[ i * 3 + 1  ];	
	   coords[i][2] = coords_0[ i * 3 + 2  ];	
	}

	//printf("ngh_0\t%d %d %d\n", ngh_0[0], ngh_0[1], ngh_0[2]); 
	//printf("ngh\t%d %d %d\n", ngh[0][0], ngh[0][1], ngh[0][2]); 
	//printf("coords_0\t%f %f %f\n", coords_0[0], coords_0[1], coords_0[2]);
	//printf("coords\t%f %f %f\n", coords[0][0], coords[0][1], coords[0][2]);

	surf_dist(inside_label, subsample, subsample_factor, n,  mask_fn, output_fn, VERBOSE, coords, nngh, ngh);

	Py_RETURN_NONE;
}*/

static char wm_dist_docs[] = "wm_dist( ): Calculate minimum distances through WM\n";
//static char surf_dist_docs[] = "surf_dist( ): Calculate minimum geodesic distance across surface\n";

static PyMethodDef BrainDist_funcs[] = {
  {"wm_dist", (PyCFunction) BrainDist_wm_dist, METH_VARARGS, wm_dist_docs},
//   {"surf_dist", (PyCFunction) BrainDist_surf_dist, METH_VARARGS, surf_dist_docs},
   {NULL}
};



//#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
       PyModuleDef_HEAD_INIT,
        "BrainDist",     // m_name 
        "This is the wm_dist module",  // m_doc
        -1,                  // m_size 
        BrainDist_funcs,    // m_methods
        NULL,                // m_reload 
        NULL,                // m_traverse 
        NULL,                // m_clear 
        NULL,                // m_free 
    };

PyMODINIT_FUNC
PyInit_BrainDist(void){
	PyObject* m;
	m = PyModule_Create(&moduledef);
	import_array();
	return m;
}

