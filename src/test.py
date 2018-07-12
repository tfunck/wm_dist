import BrainDist
import h5py
import numpy as np


def safe_h5py_open(filename, mode):
    '''open hdf5 file, exit elegantly on failure'''
    try :
        f = h5py.File(filename, mode)
        return f
    except OSError :
        print('Error: Could not open', filename)
        exit(1)


img_fn = "CTRL_C01_pve_classify.mnc"

img_vol_h5py = safe_h5py_open(img_fn, 'r')

img_vol = np.array(img_vol_h5py['minc-2.0/']['image']['0']['image'])



mesh_fn = "CTRL_C01_white_surface_left_81920.obj"

label=3

matrix_fn = "matrix.csv"

surface_mask_fn=""

example_fn=""

density_fn=""

wm_search_depth=3

max_threads=8

write_vertex=0

subsample=0

subsample_factor=0

VERBOSE=1

zstep=ystep=xstep=1

zmax= img_vol.shape[0]
ymax= img_vol.shape[1]
xmax= img_vol.shape[2]
img_vol = img_vol.astype(np.int).flatten()

starts = np.array([-72., -126., -90.]).astype(np.double)

BrainDist.wm_dist( img_vol, mesh_fn, matrix_fn, surface_mask_fn, example_fn, density_fn,label, wm_search_depth, max_threads, write_vertex, subsample, subsample_factor, zstep, ystep, xstep, zmax, ymax, xmax, starts,  VERBOSE  )
