from distutils.core import setup, Extension
import numpy as np

# the c++ extension module
#extension_mod = Extension("c_wm_dist", ["c_wm_dist.c", "c_dist-4.0.c"])

setup(name='c_wm_dist', version='1.0',   ext_modules=[Extension('c_wm_dist', ["c_wm_dist.c", "wm_dist-4.0.c"], include_dirs=[np.get_include()])] 
        author="Thomas Funck",
        author_email="tffunck@gmail.com",
        description="Calculates distances through volume with Eikonal equation",
        #long_description=long_description,
        #long_description_content_type="text/markdown",
        url="https://github.com/tfunck/wm_dist",
        packages=setuptools.find_packages(),
        classifiers=(
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ),
        
        
        )


