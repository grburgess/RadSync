#include "radsync.h"



double Source(double x, double gamma_min,double gamma_max, double index)
{
 
	if((x<gamma_min) || (x>gamma_max)){
		return 0.;
	}
	else{
		return pow(x,-index);
	}
}

std::vector<double> cooling_equation(double A, double gamma_min, double gamma_max, double index, double tcool, int N)
{
  int i,j;
  
  


  std::vector<double> gamma(N), fgamma(N), fgammatp1(N), G(N+1);
  double gdotp,gdotm;
  double V2,V3;
  double deltaGamma, step;
  double cool;

  cool = 25352.525*A*A;
  step = exp(1./N*log(gamma_max*1.1));
  
  for(i=0;i<N;i++){
    
    gamma[i]=pow(step,i);
    if(i<N-1){
      G[i]=0.5*(gamma[i]+gamma[i]*step);
    }
    else{
      G[N-1]=0.5*(gamma[i]+gamma[i]*step);
    }
    fgamma[i]=pow(gamma[i],-index);
    if((gamma[i]<gamma_min) || (gamma[i]>gamma_max)){
      fgamma[i]=0.;
    }
    
  }

  deltaGamma = G[N-1]-G[N-2];
  fgammatp1[N-1]=fgamma[N-1]/(1. + (tcool*cool*gamma[N-1]*gamma[N-1] )/ deltaGamma);

  for(j = N-2; j>=1; j--){

    deltaGamma = 0.5*(G[j]-G[j-1]); //Half steps are at j+.5 and j-.5


    //Set the coeffs. At some point this should be virtualized
    gdotp=cool*gamma[j+1]*gamma[j+1];		// Initialisation: no SSC
    gdotm=cool*gamma[j]*gamma[j];
    V3 = (tcool*gdotp )/ deltaGamma;
    V2 = 1. + (tcool*gdotm )/ deltaGamma;
    fgammatp1[j] = (fgamma[j]+Source(gamma[j],gamma_min,gamma_max,index)*tcool + V3*fgammatp1[j+1])/V2;
    
  }
  deltaGamma = G[1]-G[0];
  fgammatp1[0]=fgamma[0]+(tcool*cool*gamma[1]*gamma[1]* fgammatp1[1])/ deltaGamma;
      	

  return fgammatp1;


}
