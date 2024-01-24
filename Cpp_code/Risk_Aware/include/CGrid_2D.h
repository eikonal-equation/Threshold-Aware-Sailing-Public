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
#include <cmath>
#include <array>
#include <boost/multi_array.hpp>

//---------------------------Definitions---------------------------------------
using namespace std;
typedef boost::multi_array<double, 4> array_4D;
typedef boost::multi_array<short, 4> array_4D_short;
typedef boost::multi_array<double, 2> array_2D;
typedef boost::multi_array<double, 1> array_1D;
enum class Where_to_interp{marching, switching};

class CGrid_2D // Class declaration begins here
{
public:

	// Constructors
	CGrid_2D() = default;
	CGrid_2D(double aDr, double aRmin, int aRDim, double aDtheta, double aThetamin, int aThetaDim, double aDs)
	{
		fDr = aDr; // Delta r
		fr0 = aRmin; // starting point at r-axis
		fr = aRDim; // number of points on spatial grid = fr + 1
		fDtheta = aDtheta; // Delta theta
		ftheta0 = aThetamin; // starting point at theta-axis
		ftheta = aThetaDim; // number of theta points = ftheta + 1
		fDs = aDs; // s (threshold variable) increment
	}

	// Construct an array with length 4 for Bilinear Interpolation
	// Returns array of length 4 containing the value of the four gridpoints for bilinear interpolation
	tuple<array<double,4>, double, double> Stencil_for_Bilinear_Interp(const array_4D& aValMat, const int aMode, const int aI,
		const int aJ, const int aK, const Where_to_interp aPlace4Interp);


	// Bilinear Interpolation in 2D
	// aValArray is the output array from Stencil_for_Bilinear_Interp
	double Bilinear_Interp(const array<double, 4>& aValStencil, const double aXplus, const double aYplus,
		const double aXloc, const double aYloc, const double aDeltaX, const double aDeltaY);

	// Construct an array with length 2 for Linear Interpolation with respect to theta
   // Returns array of length 2 containing the value of the two gridpoints for linear interpolation
	array<double, 2> Stencil_for_Linear_Interp_Theta(const array_4D& aValMat, const int aMode, const int aI,
		const int aJ, const int aK);


	// Construct an array with length 2 for Linear Interpolation
	// Returns array of length 2 containing the value of the two gridpoints for linear interpolation
	array<double, 2> Stencil_for_Linear_Interp(const array_1D& aArray, const int aI);


	//Linear Interpolation in 1D
	// aValStencil is the output array from Stencil_for_Linear_Interp
	double Linear_Interp(const array<double, 2>& aValStencil, const array<double, 2>& aLocStencil, const double aXloc);


//member variables of the class
protected:
	double fDr; // Delta r
	double fr0; // starting r position
	int fr; // numeber of points of the r array  = fr + 1
	double fDtheta; // Delta theta
	double ftheta0; // starting theta position
	int ftheta; // numeber of points of the theta array  = ftheta + 1
	double fDs; //Delta s (threshold variable)
};

inline double CGrid_2D::Bilinear_Interp(const array<double, 4>& aValStencil, const double aXplus, const double aYplus,
	const double aXloc, const double aYloc, const double aDeltaX, const double aDeltaY)
{
	double DxDy = aDeltaX * aDeltaY; //dx*dy
	double pos_diff_y = aYplus - aYloc; // y(+) - yloc
	double neg_diff_y = aYloc - aYplus + aDeltaY; // yloc - y(-)
	double pos_diff_x = aXplus - aXloc; // x(+) -  xloc
	double neg_diff_x = aXloc - aXplus + aDeltaX; // xloc - x(-)

	return (pos_diff_y * (pos_diff_x * aValStencil[3] + neg_diff_x * aValStencil[2])
		+ neg_diff_y * (pos_diff_x * aValStencil[1] + neg_diff_x * aValStencil[0])) / DxDy;
}

inline double CGrid_2D::Linear_Interp(const array<double, 2>& aValStencil, const array<double, 2>& aLocStencil, const double aXloc) {
	double inc = aLocStencil[1] - aLocStencil[0];
	return aValStencil[0] + ((aValStencil[1] - aValStencil[0]) / inc) * (aXloc - aLocStencil[0]);
}
#endif // ! CGRID_1D_H
