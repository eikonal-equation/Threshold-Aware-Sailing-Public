/*=============================================================================
 * Copyright (C) 2023 MingYi Wang
 *
 * This program is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see http://www.gnu.org/licenses/.
 *============================================================================*/


 /*==============================================================================
  * File: SimpleFunctions.h
  *
  * Author: MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu
  *
  * Description: This file contains some simple helper functions and definitions.
  *
  *============================================================================*/


#pragma once
#ifndef SIMPLEFUNCTIONS_H
#define SIMPLEFUNCTIONS_H
#include <iostream>
#include <assert.h>
#include <cmath>
#include <array>
#include<numeric>

#define PI 3.141592653589793 // constant for pi
#define root_Pi  sqrt(PI)
#define inv_GR_square ((3 - sqrt(5))/2) // 1 / golden ratio^2


__inline int find_index(const double aXloc, const double aDeltaX, const double aX0, const double aXmax)
{
	//assert(aXloc <= aXmax);
	//assert(aXloc >= aX0);
	int right = int(std::ceil((aXloc + 1e-15 - aX0) / aDeltaX));
	return right;
}

__inline double mod_2pi(const double aValue) {
	return aValue - 2 * PI * std::floor(aValue / (2 * PI));
}


#endif // !SIMPLEFUNCTIONS_H

