#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include "wm_dist.h"

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

static char wm_dist_docs[] =
   "wm_dist( ): Any message you want to put here!!\n";

static PyMethodDef BrainDist_funcs[] = {
   {"wm_dist", (PyCFunction) BrainDist_wm_dist, METH_VARARGS, wm_dist_docs},
      {NULL}
};



//#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
       PyModuleDef_HEAD_INIT,
        "BrainDist",     /* m_name */
        "This is the wm_dist module",  /* m_doc */
        -1,                  /* m_size */
        BrainDist_funcs,    /* m_methods */
        NULL,                /* m_reload */
        NULL,                /* m_traverse */
        NULL,                /* m_clear */
        NULL,                /* m_free */
    };

PyMODINIT_FUNC
PyInit_BrainDist(void){
	PyObject* m;
	m = PyModule_Create(&moduledef);
	import_array();
	return m;
}


//#else
  
//  	m = PyInit_BrainDist("BrainDist", BrainDist_funcs, "This is a module");
//#endif
