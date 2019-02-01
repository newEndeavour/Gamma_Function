/*
  File:         Gamma_Function.cpp
  Version:      0.0.2
  Date:         27-Jan-2019
  Revision:     30-Jan-2019
  Author:       Jerome Drouin (jerome.p.drouin@gmail.com)

  Editions:	Please go to Gamma_Function.h for Edition Notes.

  Credits:	https://en.wikipedia.org/wiki/Gamma_function
		http://www.mymathlib.com/functions/gamma_beta.html

  Gamma_Function.cpp - Library for 'duino
  https://github.com/newEndeavour/Gamma_Function

  Gamma_Function implements the Gamma function. 

  Copyright (c) 2018-2019 Jerome Drouin  All rights reserved.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along double with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

// Includes
#include <Arduino.h>   // required for Serial.print()
#include <math.h>      // required for powl(), sinl(), fabsl() and ldexpl().
#include <float.h>     // required for DBL_MAX and LDBL_MAX
#include <limits.h>    // required for LONG_MAX
#include <float.h>     // required for DBL_EPSILON
#include <stdlib.h>    // required for malloc
#include <Gamma_Function.h>
#include <Factorial_Function.h>


// Defines
#define Nterms 20

// Library definitions
double 			Gamma_Function(double x);
long double 		xGamma_Function(long double x);
static long double 	xGamma(long double x);
static long double 	Duplication_Formula(long double two_x);

double 			Gamma_Function_Max_Arg(void); 

double 			Ln_Gamma_Function(double x);
long double 		xLn_Gamma_Function(long double x);
static long double 	xLnGamma_Asymptotic_Expansion(long double x); 

double 			DiGamma_Function(double x);
long double 		xDiGamma_Function(long double x);

static long double 	xDiGamma(long double x);
static long double 	xDiGamma_Asymptotic_Expansion(long double x);


double 			Lower_Incomplete_Gamma_Function(double x, double nu);
long double 		xLower_Incomplete_Gamma_Function(long double x, long double nu);
double 			Entire_Incomplete_Gamma_Function(double x, double nu); 
long double 		xEntire_Incomplete_Gamma_Function(long double x, long double nu);
static long double 	xSmall_x(long double x, long double nu);
static long double 	xMedium_x(long double x, long double nu);
static long double 	xLarge_x(long double x, long double nu);

// Constants
static long double const 	e 	=  2.71828182845904523536028747L;
static long double const 	pi 	=  3.14159265358979323846264338L;
static long double const 	g 	=  9.65657815377331589457187L;
static const long double 	log_sqrt_2pi 		= 9.18938533204672741780329736e-1L;
static long double const 	exp_g_o_sqrt_2pi 	= +6.23316569877722552586386e+3L;
static double 			max_double_arg 		= 171.0;
static long double 		max_long_double_arg 	= 1755.5L;
static double 			cutoff 			= 171.0;

static long double const a[] = {
                                 +1.14400529453851095667309e+4L,
                                 -3.23988020152318335053598e+4L,
                                 +3.50514523505571666566083e+4L,
                                 -1.81641309541260702610647e+4L,
                                 +4.63232990536666818409138e+3L,
                                 -5.36976777703356780555748e+2L,
                                 +2.28754473395181007645155e+1L,
                                 -2.17925748738865115560082e-1L,
                                 +1.08314836272589368860689e-4L
                              };


// Bernoulli numbers B(2),B(4),B(6),...,B(20).  Only B(2),...,B(6) currently used.
static const long double B[] = {   1.0L / (long double)(6 * 2 * 1),
                                  -1.0L / (long double)(30 * 4 * 3),
                                   1.0L / (long double)(42 * 6 * 5),
                                  -1.0L / (long double)(30 * 8 * 7),
                                   5.0L / (long double)(66 * 10 * 9),
                                -691.0L / (long double)(2730 * 12 * 11),
                                   7.0L / (long double)(6 * 14 * 13),
                               -3617.0L / (long double)(510 * 16 * 15),
                               43867.0L / (long double)(796 * 18 * 17),
                             -174611.0L / (long double)(330 * 20 * 19) 
                           };

static const int NB = sizeof(B) / sizeof(long double);



// Active /////////////////////////////////////////////////////////////////
// Gamma_Function
// This function uses Lanczos' expression to calculate Gamma(x) for real
// x, where -(max_double_arg - 1) < x < max_double_arg.                 
// Note the Gamma function is meromorphic in the complex plane and has
// poles at the nonpositive integers. 
// Tests for x a positive integer or a half positive integer give a
// maximum absolute relative error of about 1.9e-16.              
// If x > max_double_arg, then one should either use xGamma_Function(x)
// or calculate lnGamma(x). Note that for x < 0, ln (Gamma(x)) may be a complex number.

double Gamma_Function(double x)
{
long double g;

	if (x > max_double_arg) 
		return DBL_MAX;
   
	g = xGamma_Function((long double) x);

	if (fabsl(g) < DBL_MAX) 
		return (double) g;
	
	return (g < 0.0L) ? -DBL_MAX : DBL_MAX;

}

// xGamma_Function
long double xGamma_Function(long double x)
{
long double sin_x;
long double rg;
long int ix;

	// For a positive argument (x > 0)
        // if x <= max_long_double_arg return xGamma(x)
        // otherwise return LDBL_MAX.
	if (x > 0.0L) {
		if (x <= max_long_double_arg)
			return xGamma(x);
      		else 	
			return LDBL_MAX;
	}

        // For a nonpositive argument (x <= 0)
        // if x is a pole return LDBL_MAX
   	if (x > -(long double)LONG_MAX) {
		ix = (long int) x;
		if (x == (long double)ix) 
			return LDBL_MAX;
   	}
   	
	sin_x = sinl(pi * x);
	if (sin_x == 0.0L) 
		return LDBL_MAX;

        // if x is not a pole and x < -(max_long_double_arg - 1)
	// then return 0.0L
	if (x < -(max_long_double_arg - 1.0L)) 
		return 0.0L;

	// if x is not a pole and x >= -(max_long_double - 1)
	// then return Gamma(x) //
	rg = xGamma(1.0L - x) * sin_x / pi;
	if (rg != 0.0L) 
		return (1.0L / rg);
	
	return LDBL_MAX;

}


// xGamma
// xGamma uses Lanczos' expression to calculate Gamma(x) for real
// x, where 0 < x <= 900. For 900 < x < 1755.5, the duplication formula
// is used.                                                            
// The major source of relative error is in the use of the c library   
// function powl().  The results have a relative error of about 10^-16.
// except near x = 0.                                                  
static long double xGamma(long double x)
{
long double xx = (x < 1.0L) ? x + 1.0L : x;
long double temp;
int const n = sizeof(a) / sizeof(long double);
int i;

	if (x > 1755.5L) 
		return LDBL_MAX;

   	if (x > 900.0L) 
		return Duplication_Formula(x);

   	temp = 0.0L;
   	
	for (i = n-1; i >= 0; i--) {
      		temp += ( a[i] / (xx + (long double) i) );
   	}
   
	temp += 1.0L;
   	temp *= ( powl((g + xx - 0.5L) / e, xx - 0.5L) / exp_g_o_sqrt_2pi );
   	
	/*
	//DEBUG
	Serial.print("I was here...\n");
	Serial.print("x:");
	Serial.print((double)x);
	Serial.print("return:");
	long double re = (x < 1.0L) ?  temp / x : temp;
	Serial.print((double)re);
	Serial.print("\n");
	*/

	return (x < 1.0L) ?  temp / x : temp;
}


// Duplication_Formula
// Duplication_Formula returns the Gamma(two_x) using the duplication formula
// Gamma(2x) = (2^(2x-1) / sqrt(pi)) Gamma(x) Gamma(x+1/2).
static long double Duplication_Formula(long double two_x)
{
long double x = 0.5L * two_x;
long double g;
double two_n = 1.0;
int n = (int) two_x - 1;

	g  = powl(2.0L, two_x - 1.0L - (long double) n);
	g  = ldexpl(g,n);	
	g /= sqrt(pi);
	g *= xGamma_Function(x);	
	g *= xGamma_Function(x + 0.5L);
	
	return g;

}


// Gamma_Function_Max_Arg returns the maximum argument of Gamma_Function for which
// a number < DBL_MAX is returned, for arguments greater than 1.         
double Gamma_Function_Max_Arg(void) 
{ 
	return max_double_arg; 
}


// This function returns the maximum argument of xGamma_Function for which
// a number < LDBL_MAX is returned, for arguments greater than 1.
long double xGamma_Function_Max_Arg(void) 
{ 
	return max_long_double_arg; 
}


// Lower_Incomplete_Gamma_Function
// The Lower incomplete gamma function is defined as the integral from 0 to x
// of:
//		x
//		S t^(nu-1) exp(-t) dt.  
//		0
//
// The parameter nu is sometimes referred to as the shape parameter.
double Lower_Incomplete_Gamma_Function(double x, double Nu) 
{
    return (double) xLower_Incomplete_Gamma_Function((long double)x, (long double)Nu);
}


// xIncomplete_Gamma_Function
long double xLower_Incomplete_Gamma_Function(long double x, long double nu)
{

	if (x == 0.0L) 
		return 0.0L;

	if (nu<=Gamma_Function_Max_Arg()) {
		
		return xEntire_Incomplete_Gamma_Function(x,nu) * xGamma_Function(nu);
	}
	else
      		return expl(logl(xEntire_Incomplete_Gamma_Function(x,nu))
                                                  + xLn_Gamma_Function(nu));
}



// Ln_Gamma_Function
// Ln_Gamma_Function calculates the natural log of Gamma(x) for positive real x.
// Assuming that Gamma_Function_Max_Arg() = 171, then
// If 0 < x <= 171, then ln(gamma(x)) is calculated by taking the natural
// log of the result from Gamma_Function(x).  
// If x > 171, then ln(gamma(x)) is calculated using the asymptotic expansion
// 	ln(gamma(x)) ~ ln(2sqrt(2pi)) - x + (x - 1/2) ln x +              
// 			Sum B[2j] / [ 2j * (2j-1) * x^(2j-1) ], 
// 	summed over j from 1 to 3 and where B[2j] is the 2j-th Bernoulli number.
//
double Ln_Gamma_Function(double x)
{

	// For a positive argument, 0 < x <= Gamma_Function_Max_Arg()
       	// then  return log Gamma(x).
   	if (x <= Gamma_Function_Max_Arg()) 
		return log(Gamma_Function(x));

    	// otherwise return result from asymptotic expansion of ln Gamma(x). //
   	return (double) xLnGamma_Asymptotic_Expansion( (long double) x );
}


// xLn_Gamma_Function
long double xLn_Gamma_Function(long double x)
{
	// For a positive argument, 0 < x <= Gamma_Function_Max_Arg()
       	// then return log Gamma(x).
   	if (x <= Gamma_Function_Max_Arg()) 
		return logl(xGamma_Function(x));

    	// otherwise return result from asymptotic expansion of ln Gamma(x).
   	return xLnGamma_Asymptotic_Expansion( x );
}


// xLnGamma_Asymptotic_Expansion
// xLnGamma_Asymptotic_Expansion estimates log(gamma(x)) by evaluating the asymptotic
// expression:
// 	ln(Gamma(x)) ~ ln(2sqrt(2pi)) - x + (x - 1/2) ln x +
//                        Sum B[2j] / [ 2j * (2j-1) * x^(2j-1) ], 
//	summed over j from 1 to 3 and where B[2j] is the 2j-th Bernoulli number.
static long double xLnGamma_Asymptotic_Expansion(long double x) 
{
const int  m = 3;
long double term[3];
long double sum = 0.0L;
long double xx = x * x;
long double xj = x;
long double lngamma = log_sqrt_2pi - xj + (xj - 0.5L) * logl(xj);
int i;

	for (i = 0; i < m; i++) {
		term[i] = B[i] / xj; 
		xj *= xx; 
	}
   	
	for (i = m - 1; i >= 0; i--) 
		sum += term[i]; 

 	return lngamma + sum;
}


//-----------------------------


// DiGamma_Function                                      
// This function uses the derivative of the log of Lanczos's expression   
// for the Gamma function to calculate the DiGamma function for x > 1 and 
// x <= cutoff = 171.  An asymptotic expression for the DiGamma function  
// for x > cutoff.  The reflection formula                                
//                   DiGamma(x) = DiGamma(1+x) - 1/x                       
// for 0 < x < 1. and the reflection formula                              
//            DiGamma(x) = DiGamma(1-x) - pi * cot(pi*x)                  
// for x < 0.                                                             
// The DiGamma function has singularities at the nonpositive integers.    
// At a singularity, DBL_MAX is returned.                                 
double DiGamma_Function(double x)
{
long double psi = xDiGamma_Function((long double) x);

	if (fabsl(psi) < DBL_MAX) 
		return (double) psi;

   	return (psi < 0.0L) ? -DBL_MAX : DBL_MAX;
}


// xDiGamma_Function
long double xDiGamma_Function(long double x)
{
long double sin_x, cos_x;
long int ix;

	// For a positive argument (x > 0)                
        // if x <= cutoff return Lanzcos approximation
        // otherwise return Asymptotic approximation

   	if ( x > 0.0L )
      		if (x <= cutoff)
         		if (x >= 1.0L) 
				return xDiGamma(x);
         		else 
				return xDiGamma( x + 1.0L ) - (1.0L / x);
      		else 
			return xDiGamma_Asymptotic_Expansion(x);

        // For a nonpositive argument (x <= 0)
        // If x is a singularity then return LDBL_MAX.
	if ( x > -(long double)LONG_MAX) {
      		ix = (long int) x;
      		if ( x == (long double) ix) return LDBL_MAX;
   	}

	sin_x = sinl(pi * x);
   	if (sin_x == 0.0L) 
		return LDBL_MAX;
   
	cos_x = cosl(pi * x);
   	if (fabsl(cos_x) == 1.0L) 
		return LDBL_MAX;

        // If x is not a singularity then return DiGamma(x)
   	return xDiGamma(1.0L - x) - pi * cos_x / sin_x;
}


// xDiGamma
// This function uses the derivative of the log of Lanczos's expression   
// for the Gamma function to calculate the DiGamma function for x > 1 and 
// x <= cutoff = 171.                                                     
static long double xDiGamma(long double x) 
{
long double lnarg = (g + x - 0.5L);
long double temp;
long double term;
long double numerator = 0.0L;
long double denominator = 0.0L;
int const n = sizeof(a) / sizeof(long double);
int i;

	for (i = n-1; i >= 0; i--) {
      		temp 		= x + (long double) i;
      		term 		= a[i] / temp;
      		denominator    += term;
      		numerator      += term / temp;
   	} 
   
	denominator += 1.0L;
   
   	return logl(lnarg) - (g / lnarg) - numerator / denominator;

}


// xDiGamma_Asymptotic_Expansion
//     This function estimates DiGamma(x) by evaluating the asymptotic        //
//     expression:                                                            //
//         DiGamma(x) ~ ln(x) - (1/2) x +                                     //
//                        Sum B[2j] / [ 2j * x^(2j) ], summed over            //
//     j from 1 to 3 and where B[2j] is the 2j-th Bernoulli number.           //
static long double xDiGamma_Asymptotic_Expansion(long double x ) {
const int  m = 3;
long double term[3];
long double sum = 0.0L;
long double xx = x * x;
long double xj = x;
long double digamma = logl(xj) - 1.0L / (xj + xj);
int i;

	xj = xx;
   	for (i = 0; i < m; i++) { 
		term[i] = B[i] / xj; 
		xj *= xx; 
	}
   
	for (i = m - 1; i >= 0; i--) 
		sum += term[i]; 
   
	return digamma - sum;
}


//-----------------------------


// Entire_Incomplete_Gamma_Function
// The entire incomplete gamma function, also called the regularized
// incomplete gamma function, is defined as the integral from 0 to x of
// the integrand t^(nu-1) exp(-t) / gamma(nu) dt.  The parameter nu is
// sometimes referred to as the shape parameter.
double Entire_Incomplete_Gamma_Function(double x, double nu)
{
	return (double) xEntire_Incomplete_Gamma_Function((long double)x,(long double)nu);
}


// xEntire_Incomplete_Gamma_Function
long double xEntire_Incomplete_Gamma_Function(long double x, long double nu)
{

	if (x == 0.0L) 
		return 0.0L;

	if (fabsl(x) <= 1.0L) 
		return xSmall_x(x, nu);

	if (fabsl(x) < (nu + 1.0L) ) 
		return xMedium_x(x, nu);

	return xLarge_x(x, nu);
}


// xSmall_x
// This function approximates the entire incomplete gamma function for
// x, where -1 <= x <= 1.
static long double xSmall_x(long double x, long double nu)
{
long double terms[Nterms];
long double x_term = -x;
long double x_power = 1.0L;
long double sum;
int i;

	for (i = 0; i < Nterms; i++) {
      		terms[i] = (x_power / xFactorial(i)) / (i + nu);
      		x_power *= x_term;
   	}
	sum = terms[Nterms-1];
   
	for (i = Nterms-2; i >= 0; i--) 
		sum += terms[i];
   
	if ( nu <= Gamma_Function_Max_Arg() )
      		return powl(x,nu) * sum / xGamma_Function(nu);
   	else 
		return expl(nu * logl(x) + logl(sum) - xLn_Gamma_Function(nu));

}


// xMedium_x
// This function approximates the entire incomplete gamma function for
// x, where 1 < x < nu + 1.
// If nu + 1 < x, then one should use xLarge_x(x,nu).
static long double xMedium_x(long double x, long double nu)
{
long double coef;
long double term = 1.0L / nu;
long double corrected_term = term;
long double temp_sum = term;
long double correction = -temp_sum + corrected_term;
long double sum1 = temp_sum;
long double sum2;
long double epsilon = 0.0L;
long double EpsCorrection = 1E-9;
int i;


	if (nu > Gamma_Function_Max_Arg()) {
      		coef = expl( nu * logl(x) - x - xLn_Gamma_Function(nu) );
      		if (coef > 0.0L) 
			epsilon = DBL_EPSILON/coef;
   	} else {
		coef = powl(x, nu) * expl(-x) / xGamma_Function(nu);
      		epsilon = DBL_EPSILON/coef;
   	}
   
	if (epsilon <= 0.0L) 
		epsilon = (long double) DBL_EPSILON;

   	for (i = 1; term > epsilon * sum1; i++) {
      		term *= x / (nu + i);
      		corrected_term = term + correction;
      		temp_sum = sum1 + corrected_term;
      		correction = (sum1 - temp_sum) + corrected_term;
      		sum1 = temp_sum;
   	}

   	sum2 = sum1;
   	sum1 *= coef;
   	correction += sum2 - sum1 / coef;
   	term *= x / (nu + i);
   	sum2 = term + correction;
   
	//Condition too stringent for not much improvement
	//for (i++; (term + correction) > (epsilon * sum2); i++) {
	for (i++; ((term + correction) - (epsilon * sum2)) > EpsCorrection; i++) {
      		term 		*= x / (nu + i);
      		corrected_term   = term + correction;
      		temp_sum 	 = sum2 + corrected_term;
      		correction 	 = (sum2 - temp_sum) + corrected_term;
      		sum2 		 = temp_sum;
   	}
   
   	sum2 += correction;
   	sum2 *= coef;
   
	return sum1 + sum2;
}


//
// This function approximates the entire incomplete gamma function 
// for x, where nu + 1 <= x.
// If 0 <= x < nu + 1, then one should use xSmall_x(x,nu).
static long double xLarge_x(long double x, long double nu)
{
long double temp = 1.0L / nu;
long double sum = temp;
long double coef;
long double medi;
int i = 0;
int n;

	n = (int)(x - nu - 1.0L) + 1;
   
	for (i = 1; i < n; i++) {
      		temp *= x / (nu + i);
      		sum += temp;
   	}
		
	if (nu <= Gamma_Function_Max_Arg() ) {
      		coef = powl(x, nu) * expl(-x) / xGamma_Function(nu);
		medi = xMedium_x(x, nu + n);
		return medi + coef * sum;
   	} else {
      		return expl(logl(sum) + nu * logl(x) - x - xLn_Gamma_Function(nu)) +
                                                        xMedium_x(x, nu + n); 
   	}   

}



// Retired /////////////////////////////////////////////////////////////////
// Lanczos approximation to the gamma function. 
double Retired_Gamma(double x) 
{
    double ret = (1.000000000190015 + 
                 76.18009172947146 / (x + 1) +  
                 -86.50532032941677 / (x + 2) + 
                 24.01409824083091 / (x + 3) +  
                 -1.231739572450155 / (x + 4) + 
                 1.208650973866179e-3 / (x + 5) + 
                 -5.395239384953e-6 / (x + 6));
    
    return ret * sqrt(2*M_PI)/x * pow(x + 5.5, x+.5) * exp(-x-5.5);
}


// Natural log of the gamma function Gammaln
double Retired_Gammaln(double xx)
{
    double x, tmp, ser;

    const static double cof[6]={76.18009172947146,    -86.50532032941677,
                                24.01409824083091,    -1.231739572450155,
                                0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    x	 = xx-1.0;
    tmp	 = x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser	 = 1.000000000190015;

    for (j=0;j<=5;j++) {
        x   += 1.0;
        ser += cof[j] / x;
    }

    return -tmp + log(2.5066282746310005 * ser);
}


