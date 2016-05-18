#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>
#include <stdlib.h>
#include <vector>

double Source(double x, double gamma_min,double gamma_max, double index);

std::vector<double> cooling_equation(double A, double gamma_min, double gamma_max, double index, double tcool, int N);
