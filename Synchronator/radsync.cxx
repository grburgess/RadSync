#include "radsync.h"



double Source(double x, double ne ,double gamma_min,double gamma_max, double p)
{
 
	if((x<gamma_min) || (x>gamma_max)){
		return 0.;
	}
	else{
	  return ne * pow(x,-p) * (1-p)* 1./ (pow(gamma_max,1-p)-pow(gamma_min,1-p));
	}
}

std::vector<double> cooling_equation(double ne, double A, double gamma_min, double gamma_max, double index, double DT, int N,int steps)
{
  int i,j;
  
  
  //N=100;

  std::vector<double> gamma(N), fgamma(N), fgammatp1(N), G(N+1), total(N);
  double gdotp,gdotm;
  double V2,V3;
  double deltaGamma, step;
  double cool;
  //A=.0001;
  cool = 25352.525*A*A;
  step = exp(1./N*log(gamma_max*1.1));

  //DT = 1E-5;
  //cool= 0.005;

  for(j=0;j<N;j++)
    {
    
      gamma[j]=pow(step,j);
      if(j<N-1)
	{
	  G[j]=0.5*(gamma[j]+gamma[j]*step);
	}
      else
	{
	  G[N-1]=0.5*(gamma[j]+gamma[j]*step);
	}
      fgamma[j]=0.;//Source(gamma[j],gamma_min,gamma_max,index);
       
        total[j]=fgamma[j];
    }
  


  for(i=0;i<steps;i++)
    {

      // We will walk backwards so we must set the end points
      
      deltaGamma = G[N-1]-G[N-2];
      fgammatp1[N-1]=fgamma[N-1]/(1. + (DT*cool*gamma[N-1]*gamma[N-1] )/ deltaGamma);

      // Try to fix the [j=0] point first
      //fgammatp1[0]=;

      for(j = N-2; j>=1; j--)
	{


	  
	  deltaGamma = 0.5*(G[j]-G[j-1]); //Half steps are at j+.5 and j-.5
	  
	  
	  
	  gdotp=cool*gamma[j+1]*gamma[j+1]; // Forward  step cooling 
	  gdotm=cool*gamma[j]*gamma[j];     // Backward step cooling

	  V3 = (DT*gdotp )/ deltaGamma;     // Tridiagonal coeff.
	  V2 = 1. + (DT*gdotm )/ deltaGamma;// Tridiagonal coeff.

	 
	   


	  // Solve for forward electron distribution
	  fgammatp1[j] = (fgamma[j]+Source(gamma[j], ne, gamma_min,gamma_max,index)*DT + V3*fgammatp1[j+1])/V2;
	  //fgammatp1[j] = (fgamma[j] + V3*fgammatp1[j+1])/V2;
	 
	 
	}

      // Set the end point 
      deltaGamma = G[1]-G[0];
      fgammatp1[0]=fgamma[0]+(DT*cool*gamma[1]*gamma[1]* fgammatp1[1])/ deltaGamma;
      // fgammatp1[0]=fgamma[0]+(DT*cool*gamma[1]*gamma[1]* fgammatp1[1])/ deltaGamma;



      // Set fgamma to the soultion of the energy
      // equation from this run.
      

      for(j=0;j<N;j++)
	{
	  
	  fgamma[j]=fgammatp1[j];
	  if(fgamma[j]<1e-30)
	    {
	      fgamma[j]=0.;
	    }
	  //printf("%e\n",fgamma[j]);
	  total[j]+=fgamma[j];
	  //printf("%e\n\n",total[j]);
	 
	}
      
      // Now fgamma will be cooled in the round
	
    }
  return fgamma;
  //return total;
}










std::vector<double> synchrotron(double * energy, double A, double * gamma, double * gamma2  ,double * fgamma, int N, int Nene)
{
  std::vector<double> out_val(Nene);
  double h,s, ec, y[N];
  double syncArg, sum;
  int i,j,k, status;
  gsl_sf_result result;
  double Bcritical = 4.14E13; // Gauss
   gsl_set_error_handler_off();
  
    ec = 1.5*A/Bcritical; //Put proper units in!

    

  for(i=0;i<Nene;i++)
    {
      for(j=0;j<N;j++)
	{
	  syncArg = energy[i]/(ec * gamma2[j] );
	  status = gsl_sf_synchrotron_1_e(syncArg, &result);
	    
	  //if (status)
	  //  {
	      if (status == GSL_ERANGE)
		{
		  y[j]=0.;
		  
		}
	      else
		{
		  //printf("%e\n",result.val);
		  y[j] = fgamma[j] * result.val;	
		}
	      //  }
	  
	  //y[j] = fgamma[j] * gsl_sf_synchrotron_1(syncArg);
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





////// synch + electrons


std::vector<double> coolAndSynch(double * energy, int Nene, double ne, double A, double bulkGamma, double gamma_min, double gamma_max, double index, double DT, int N,int steps)
{
  int i,j,k;
  
  


  std::vector<double> emission(Nene), emissiontmp(Nene);
  // std::vector<double> gamma(N), fgamma(N), fgammatp1(N), G(N+1), total(N);

  double gamma[N], gamma2[N], fgamma[N], fgammatp1[N], G[N+1];
  double gdotp,gdotm;
  double V2,V3;
  double deltaGamma, step;
  double cool;
  //A=.1;
  cool = 25352.525*A*A;
  step = exp(1./N*log(gamma_max*1.1));

  //DT = 0.0005;
  //cool= 0.005;

  for(j=0;j<N;j++)
    {
    
      gamma[j]=pow(step,j);
      gamma2[j]=gamma[j]*gamma[j];
      
      if(j<N-1)
	{
	  G[j]=0.5*(gamma[j]+gamma[j]*step);
	}
      else
	{
	  G[N-1]=0.5*(gamma[j]+gamma[j]*step);
	}
      fgamma[j]=0;//Source(gamma[j],gamma_min,gamma_max,index);
     
    }
  


  for(i=0;i<steps;i++)
    {

      // We will walk backwards so we must set the end points
      
      deltaGamma = G[N-1]-G[N-2];
      fgammatp1[N-1]=fgamma[N-1]/(1. + (DT*cool*gamma[N-1]*gamma[N-1] )/ deltaGamma);

      // Try to fix the [j=0] point first
      //fgammatp1[0]=;

      for(j = N-2; j>=1; j--)
	{


	  
	  deltaGamma = 0.5*(G[j]-G[j-1]); //Half steps are at j+.5 and j-.5
	  
	  
	  
	  gdotp=cool*gamma[j+1]*gamma[j+1]; // Forward  step cooling 
	  gdotm=cool*gamma[j]*gamma[j];     // Backward step cooling

	  V3 = (DT*gdotp )/ deltaGamma;     // Tridiagonal coeff.
	  V2 = 1. + (DT*gdotm )/ deltaGamma;// Tridiagonal coeff.

	 
	   


	  // Solve for forward electron distribution
	  fgammatp1[j] = (fgamma[j]+Source(gamma[j],ne,gamma_min,gamma_max,index)*DT + V3*fgammatp1[j+1])/V2;
	  //fgammatp1[j] = (fgamma[j] + V3*fgammatp1[j+1])/V2;
	 
	 
	}

      // Set the end point 
      deltaGamma = G[1]-G[0];
      fgammatp1[0]=fgamma[0]+(DT*cool*gamma[1]*gamma[1]* fgammatp1[1])/ deltaGamma;
     

      // Set fgamma to the soultion of the energy
      // equation from this run.
      

      for(j=0;j<N;j++)
	{
	  
	  fgamma[j]=fgammatp1[j];
	  if(fgamma[j]<1e-30)
	    {
	      fgamma[j]=0.;
	    }

	 
	}
      
      emissiontmp = synchrotron(energy,A*bulkGamma,gamma,gamma2,fgamma,N,Nene);
      for(k=0;k<Nene;k++)
	{
	  emission[k]+= emissiontmp[k];
	}
      
      // Now fgamma will be cooled in the round
	
    }
  return emission;
  //return total;
}
