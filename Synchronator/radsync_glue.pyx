cimport cython
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
#from pygsl.sf import synchrotron_1
#import pygsl.errors

cdef extern from "/usr/local/include/radsync.h":
    cdef vector[double] cooling_equation(double, double,double, double, double, double, int, int)
    cdef vector[double] synchrotron(double *, double, double *, double *, double *, int, int)
    cdef vector[double] coolAndSynch(double *, int,double, double, double, double, double, double,double, int, int)

def electrons(ne, A,  gamma_min,  gamma_max,  index,  tcool, N ,steps):
    return cooling_equation(ne, A,  gamma_min,  gamma_max,  index,  tcool, N, steps )

def emission(np.ndarray[double, ndim=1, mode="c"] ene not None, double ne, A, bulkGamma, gamma_min, gamma_max, index,DT, N, steps):
    return coolAndSynch(<double*> np.PyArray_DATA(ene),
                        ene.shape[0],
                        ne,
                        A,
                        bulkGamma,
                        gamma_min,
                        gamma_max,
                        index,
                        DT,
                        N,
                        steps
                        )




@cython.cdivision(True)
def synchrotronPy(np.ndarray[double, ndim=1, mode="c"] ene not None, A,np.ndarray[double, ndim=1, mode="c"] gamma not None, np.ndarray[double, ndim=1, mode="c"] fgamma not None, N ):
    Nene = ene.shape[0]
    gamma2 = np.power(gamma,2)
    
    val = synchrotron(<double*> np.PyArray_DATA(ene),
                      A,
                      <double*> np.PyArray_DATA(gamma),
                       <double*> np.PyArray_DATA(gamma2),
                        <double*> np.PyArray_DATA(fgamma),
                        N,
                        Nene)
    
    return np.array(val)
    
#@cython.cdivision(True)
# def synchrotron_old(ene, double A, gamma,fgamma, int N):
#     ene=np.array(ene,ndmin=1,copy=False)
#     val = np.zeros_like(ene)
#     cdef double h
#     cdef double s
#     for j,e in enumerate(ene):
#         y=np.zeros(N)
#         for i in range(N):
#             try:
#                 y[i] = fgamma[i]*synchrotron_1(e/(A*gamma[i]*gamma[i]))[0]
#             except pygsl.errors.gsl_Error, err:
#                 pass
    
#         s = y[0] + 2.0 * y[1:N-1].sum() + y[N-1]
#         h = float(gamma[-1] - gamma[0]) / N
            
#         val[j] = s * h / (2.0*e)
#     return A*val
               
