#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_result.h>
#include <stdlib.h>
#include <vector>
#include <stdio.h>
double Source(double x, double gamma_min,double gamma_max, double index);

std::vector<double> cooling_equation(double ne,double A, double gamma_min, double gamma_max, double index, double tcool, int N,int steps);
std::vector<double> synchrotron(double * energy, double A, double * gamma, double * gamma2  ,double * fgamma, int N, int Nene);
std::vector<double> coolAndSynch(double * energy, int Nene, double ne, double A, double bulkGamma, double gamma_min, double gamma_max, double index, double DT, int N,int steps);
