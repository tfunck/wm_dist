from distutils.core import setup, Extension

# the c++ extension module
extension_mod = Extension("c_wm_dist", ["c_wm_dist.c", "c_dist-4.0.c"])

setup(name='c_wm_dist', version='1.0',   ext_modules=[Extension('c_wm_dist', ["c_wm_dist.c", "wm_dist-4.0.c"])])


