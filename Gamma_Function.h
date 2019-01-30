/*
  File:         Gamma_Function.h
  Version:      0.0.2
  Date:         27-Jan-2019
  Revision:     30-Jan-2019
  Author:       Jerome Drouin (jerome.p.drouin@gmail.com)

  Credits:	https://en.wikipedia.org/wiki/Gamma_function
		http://www.mymathlib.com/functions/gamma_beta.html

  Gamma_Function.h - Library for 'duino
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
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Editions:	
	0.0.1	: first version
	0.0.2	: Correction of a few bugs in .cpp (long-> long double and other casting 
		  niceties)
	
*/


#ifndef GAMMA_FUNCTION_H
#define GAMMA_FUNCTION_H

double Gamma_Function(double x);
double Lower_Incomplete_Gamma_Function(double x, double Nu);


#endif
