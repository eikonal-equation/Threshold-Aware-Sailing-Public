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
  * File: CGrid_2D.h
  *
  * Author: MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu
  *
  * Description: This file contains the declarations of the class "CGrid_2D"
  * that sets up a unifrom grid in 2D as well as functions for
  * bilinear interpolation and linear interpolation.
  *
  * Details of all of these functions are found in CGrid_2D.cpp.
  *
  *============================================================================*/
#pragma once
#ifndef GRID_2D_H
#define CGRID_2D_H

  //---------------------------Libraries-----------------------------------------
#include <iostream>
#include <cmath>
#include <array>
#include <boost/multi_array.hpp>

//---------------------------Definitions---------------------------------------
using namespace std;
typedef boost::multi_array<double, 3> array_3D;
typedef boost::multi_array<double, 1> array_1D;


class CGrid_2D // Class declaration begins here
{
public:

	// Constructors
	CGrid_2D() = default;
	CGrid_2D(double  aDr, double  armin, int arDim, double  aDtheta, double  athetamin, int athetaDim)
	{
		fDr = aDr; // Delta r
		fr0 = armin; // starting point at r-axis
		fr = arDim; // number of points on spatial grid = fr + 1
		fDtheta = aDtheta; // Delta theta
		ftheta0 = athetamin; // starting point at theta-axis
		ftheta = athetaDim; // number of theta points = ftheta + 1

	}

	// Construct an array with length 4 for Bilinear Interpolation
	// Returns array of length 4 containing the value of the four gridpoints for bilinear interpolation
	array<double, 4> Stencil_for_Bilinear_Interp(const array_3D& aValMat, const int aQ, const int aIr, const int aItheta);

	// Bilinear Interpolation in 2D
	// aValArray is the output array from Stencil_for_Bilinear_Interp
	double  Bilinear_Interp(const array<double, 4>& aValArray, const int aIr, const int aItheta,
		const double  aRloc, const double  aThetaloc);

	// Construct an array with length 2 for Linear Interpolation with respect to theta
	// Returns array of length 2 containing the value of the two gridpoints for linear interpolation
	array<double, 2> Stencil_for_Linear_Interp_Theta(const array_3D& aValMat, const int aQ,
		const int aIr, const int aItheta);
	// Linear Interpolation in 1D in the Theta variable
	// aValArray is the output array from Stencil_for_Linear_Interp_Theta
	double  Linear_Interp_Theta(const array<double, 2>& aValArray, const int aIndex, const double  aThetaloc);

	// Construct an array with length 2 for Linear Interpolation
	// Returns array of length 2 containing the value of the two gridpoints for linear interpolation
	array<double, 2> Stencil_for_Linear_Interp(const array_1D& aArray, const int aI);

	// Linear Interpolation in 1D
	// aValStencil is the output array from Stencil_for_Linear_Interp
	double  Linear_Interp(const array<double, 2>& aValStencil, const array<double, 2>& aLocStencil, const double  aXloc);

//member variables of the class
protected:
	double  fDr; // Delta r
	double  fr0; // starting r position
	int fr; // numeber of points of the r array  = fr + 1
	double  fDtheta; // Delta theta
	double  ftheta0; // starting theta position
	int ftheta; // numeber of points of the theta array  = ftheta + 1
};
#endif // !CGRID_1D_H
