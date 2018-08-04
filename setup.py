import setuptools
from distutils.core import setup, Extension
import numpy as np
#, "wm_dist-4.0.c"

long_description="wm_dist calculates distances from points on a mesh through a volume. Distances are calculated by solving the Eikonal equation.\nThe primary use case for wm_dist is brain imaging for calculating distances through the cortical white matter from vertex points defined on a mesh representation of the cortical gray matter. However, in principal, any distances can be calculated through any volume from any points defined on an surface mesh."
setup(name='wm_dist', version='1.0.2',   ext_modules=[Extension('c_wm_dist', ["wm_dist.c"], include_dirs=[np.get_include()])], 
        author="Thomas Funck",
        author_email="tffunck@gmail.com",
        description="Calculates distances through volume with Eikonal equation",
        long_description=long_description,
        #long_description_content_type="text/markdown",
        url="https://github.com/tfunck/wm_dist",
        packages=setuptools.find_packages(),
        classifiers=(
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            #"Operating System :: POSIX :: Linux",
        ),
        
        
        )


