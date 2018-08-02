from distutils.core import setup, Extension
from distutils.core import setup, Extension

# the c++ extension module
extension_mod = Extension("BrainDist", ["BrainDist.c", "wm_dist-4.0.c"])

setup(name='BrainDist', version='1.0',   ext_modules=[Extension('BrainDist', ["BrainDist.c", "wm_dist-4.0.c"])])


