#include "radsync.h"



double Source(double x, double gamma_min,double gamma_max, double p)
{
 
	if((x<gamma_min) || (x>gamma_max)){
		return 0.;
	}
	else{
	  return pow(x,-p) * (1-p)* 1./ (pow(gamma_max,1-p)-pow(gamma_min,1-p));
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
  
  for(i=0;i<N;i++)
    {
    
      gamma[i]=pow(step,i);
      if(i<N-1)
	{
	  G[i]=0.5*(gamma[i]+gamma[i]*step);
	}
      else
	{
	  G[N-1]=0.5*(gamma[i]+gamma[i]*step);
	}
      fgamma[i]=Source(gamma[i],gamma_min,gamma_max,index);
    
    
    }

  deltaGamma = G[N-1]-G[N-2];
  fgammatp1[N-1]=fgamma[N-1]/(1. + (tcool*cool*gamma[N-1]*gamma[N-1] )/ deltaGamma);

  for(j = N-2; j>=1; j--)
    {
      
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

std::vector<double> synchrotron(double * energy, double A, double * gamma, double * gamma2  ,double * fgamma, int N, int Nene)
{
  std::vector<double> out_val(Nene);
  double h,s, ec, y[N];
  double syncArg, sum;
  int i,j,k, status;
  gsl_sf_result result;
  double Bcritical = 4.14E13; // Gauss
  
    ec = 1.5*A/Bcritical; //Put proper units in!

  

  for(i=0;i<Nene;i++)
    {
      for(j=0;j<N;j++)
	{
	  syncArg = energy[i]/(ec * gamma2[j] );
	  y[j] = fgamma[j] * gsl_sf_synchrotron_1(syncArg);
	}

      sum=0.;
      for(k=1;k<N-1;k++)
	{
	  sum+=y[k];
	}
      
      
      s = y[0] + 2.*sum + y[N-1];
      
      h = (double)(gamma[N-1] - gamma[0])/(double)N;
      out_val[i] = s*h/(2.* energy[i] );
      
    }

  return out_val;

}


