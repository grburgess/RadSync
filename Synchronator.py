import numpy as np
import matplotlib.pyplot as plt
from radsync_glue import electrons, synchrotron


class Synchronator(object):


    def __init__(self, Ngrid=500):
        self._Ngrid = Ngrid
        self._constFactor = 23532.535

    def SetParameters(A,  gamma_min,  gamma_max,  index,  tCoolFrac):

        self._A = A
        self._gamma_min = gamma_min
        self._gamma_max = gamma_max
        self._index = index

        self._tCoolFrac = tCoolFrac
        self._syncCool = 1./(A*A*self._factor)  
        self._tCool = self._syncCool*tCoolFrac

        print "Cooling electrons...\n"
        print "Synchrotron cooling time is %.2f s"%self._syncCool
        print "Cooling electrons for %.2f s"%self._tCool
        
        self._fgammma = electrons(A,  gamma_min,  gamma_max,  index,  self._tCool, self._Ngrid )
        step = np.exp(1./N*log(gamma_max*1.1))
        self._gamma = np.power(step,arange(self._Ngrid))
        self._gammaFlag = True

    def GetEmisson(self,ene):
        if self._gammaFlag:
            return synchrotron(ene, self._A, self._gamma,self._fgammma,self._Ngrid)
        else:
            print "No Parameters Set!"
            return
    def GetGammaGrid(self):
        return self._gamma

    def GetGammaDist(self):
        return self._fgammma
   
