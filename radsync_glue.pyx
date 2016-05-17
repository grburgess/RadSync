cimport cython
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from pygsl.sf import synchrotron_1
import pygsl.errors

cdef extern from "/usr/local/include/radsync.h":
    cdef vector[double] cooling_equation(double, double, double, double, double, int)



def electrons(A,  gamma_min,  gamma_max,  index,  tcool, N ):
    return cooling_equation(A,  gamma_min,  gamma_max,  index,  tcool, N )

    
@cython.cdivision(True)
def synchrotron(ene, double A, gamma,fgamma, int N):
    ene=np.array(ene,ndmin=1,copy=False)
    val = np.zeros_like(ene)
    cdef double h
    cdef double s
    for j,e in enumerate(ene):
        y=np.zeros(N)
        for i in range(N):
            try:
                y[i] = fgamma[i]*synchrotron_1(e/(A*gamma[i]*gamma[i]))[0]
            except pygsl.errors.gsl_Error, err:
                pass
    
        s = y[0] + 2.0 * y[1:N-1].sum() + y[N-1]
        h = float(gamma[-1] - gamma[0]) / N
            
        val[j] = s * h / (2.0*e)
    return val
               
