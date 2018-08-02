
import c_wm_dist
from optparse import OptionParser, OptionGroup
import nibabel as nib
import numpy as np

def calculate(classified_fn, mesh_fn, matrix_fn, example_fn="", density_fn="", surface_mask_fn="", label=1, wm_search_depth=3, max_threads=1, example_vertex=1, subsample_list="0,0", verbose=0 ):
    img = nib.load(opts.classified_fn)  
    affine = img.get_affine()

    steps = np.array(img.header.get_zooms()).astype(np.double)
    maxima  = np.array(img.header.get_data_shape()).astype(np.int64)

    #dim_order = np.array(nib.aff2axcodes(img.affine))
    #zdim = np.argmax(dim_order == 'S' )
    #ydim = np.argmax(dim_order == 'A' )
    #xdim = np.argmax(dim_order == 'R' )

    starts = affine[[2,1,0], 3].astype(np.double)
    
    img_vol = img.get_data().astype(np.int64)
    img_vol = img_vol.flatten()
    #img_vol = [ img_vol[z][y][x] for z in range(maxima[0]) for y in range(maxima[1]) for x in range(maxima[2])  ]
    #print(np.sum(img_vol == img_vol_0), maxima[0]*maxima[1]*maxima[2])
    #img_vol = img_vol.flatten()
    subsample, subsample_factor = [ int(i) for i in  subsample_list.split(',') ]
    
    c_wm_dist.calculate(img_vol, steps, maxima, starts, opts.mesh_fn, opts.matrix_fn, opts.surface_mask_fn, opts.example_fn, opts.density_fn, np.int64(opts.label), opts.wm_search_depth, opts.max_threads, opts.example_vertex, subsample, subsample_factor,  opts.verbose  )

    return 0

if __name__ == "__main__":
    usage = "usage: "
    version=1
    parser = OptionParser(usage=usage,version=version)
    group= OptionGroup(parser,"Mandatory")
    group.add_option("-c","--classified",  dest="classified_fn", type="string", default="",  help="Classified image with integer values coding anatomic regions.")
    group.add_option("-s","--surface",  dest="mesh_fn", type="string", default="",  help="Cortical surface mesh file (.obj).")
    group.add_option("-m","--matrix",  dest="matrix_fn", type="string", default="",  help="Output matrix file with distances between vertex points.")
    parser.add_option_group(group)
    group.add_option("-l","--label", dest="label", type='int',  help="Integer label to identify anantomic region in classified image.")

    group= OptionGroup(parser,"Optional")
    group.add_option("-S","--surface-mask",  dest="surface_mask_fn", type="string", default="",  help="Surface mask for cortical surface mask (will only use non-zero labeled vertices) ")
    group.add_option("-e","--example", dest="example_fn", type='string',  default="", help="Save an exemplar volumetric image with distances from a vertex.")
    group.add_option("-d","--density", dest="density_fn", type='string', default="",  help="Save volumetric image with density of minimum path distances at each voxel.")
    group.add_option("-v","--example-vertex", dest="example_vertex", type='int', default=1, help="Vertex value (integer) for which to save an example volume")
    
    group.add_option("-w","--wm-search-depth", dest="wm_search_depth", type='int', default=3, help="Number of voxels over which to search from surface vertex point to labeled region.")
    group.add_option("-t","--max-threads", dest="max_threads", type='int', default=-1, help="Number of threads to use for parallel processing. (Default=use all available threads)")
    group.add_option("-u","--subsample", dest="subsample_list", type='string',default="0,4", help="Comma-separated list with number of times to subsample surface mesh by given factor: <subsample times>,<subsample factor> (Default=0,0)")
    group.add_option("--verbose", dest="verbose", type='int', default=1, help="Level of verbosity (0=silent, 1=verbose) (Default=1)")
   
    print("\nWarning: Dimension order of image must be Z, Y, X!\n");

    (opts, args) = parser.parse_args()

    calculate( opts.classified_fn, opts.mesh_fn, opts.matrix_fn, opts.example_fn, opts.density_fn, opts.surface_mask_fn, opts.label, opts.wm_search_depth, opts.max_threads, opts.example_vertex, opts.subsample_list, opts.verbose)




