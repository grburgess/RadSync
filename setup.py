
from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install_headers import install_headers
from distutils.command.build_clib import build_clib
from Cython.Distutils import build_ext
from Cython.Build import cythonize



import os

import numpy



os.system("g++ -shared -o Synchronator/libradsync.so -fPIC Synchronator/radsync.cxx -I/usr/local/include  -L/usr/local/lib -lgsl -lgslcblas")
os.system("cp Synchronator/libradsync.so /usr/local/lib/")
os.system("cp Synchronator/radsync.h /usr/local/include/")

ext_modules = [Extension("Synchronator/radsync_glue",["Synchronator/radsync_glue.pyx"],
            library_dirs=['/usr/local/lib'],language='c++',
            libraries=["radsync"],include_dirs = [numpy.get_include()])]


setup(
    
    name="radsync",

    
    packages = ['Synchronator'],

    include_dirs = [numpy.get_include()],                
    
    version = 'v0.0.1',
        
    author = 'J. Michael Burgess',
    
    author_email = 'jmichael.grb@gmail.com',
    
    
    ext_modules=cythonize(ext_modules),
        
      
      )


