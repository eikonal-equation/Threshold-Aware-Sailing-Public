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
  * File: CGrid_2D.cpp
  *
  * Author: MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu
  *
  * Description: This file contains the acutal implementation of member functions
  * that sets up a unifrom grid in 2D as well as functions for
  * bilinear interpolation and linear interpolation.
  *
  *============================================================================*/

//----------------------Project specific header files---------------------------
#include "CGrid_2D.h"

//---------------------------Definitions---------------------------------------
using namespace std;

// This function returns an array of 4 double that stores the value at
// topright,topleft, bottomright, and bottomleft gridpoints.
// Order of Dimensions for aValMat: mode, slice, radius, theta
array<double, 4> CGrid_2D::Stencil_for_Bilinear_Interp(const array_3D& aValMat, const int aQ,
	const int aIr, const int aItheta) {
	array<double, 4> theValArray;

	//checking if the bilinear interpolation is between ftheta-1 and ftheta,
	// if so, we will replace the value at ftheta to be the value at theta=0
	if (aItheta == ftheta) {
		theValArray[0] = aValMat[aQ - 1][aIr][0]; // value at theta = 0 (instead of theta = 2pi)
		theValArray[1] = aValMat[aQ - 1][aIr][aItheta - 1];//topleft
		theValArray[2] = aValMat[aQ - 1][aIr - 1][0]; //value at theta = 0 (instead of theta = 2pi)
		theValArray[3] = aValMat[aQ - 1][aIr - 1][aItheta - 1];//bottomleft
	}
	else if (aItheta == 0) {
		theValArray[0] = aValMat[aQ - 1][aIr][aItheta];
		theValArray[1] = aValMat[aQ - 1][aIr][ftheta - 1];//topleft
		theValArray[2] = aValMat[aQ - 1][aIr - 1][aItheta];
		theValArray[3] = aValMat[aQ - 1][aIr - 1][ftheta - 1];//bottomleft
	}
	else {
		theValArray[0] = aValMat[aQ - 1][aIr][aItheta];// topright
		theValArray[1] = aValMat[aQ - 1][aIr][aItheta - 1];//topleft
		theValArray[2] = aValMat[aQ - 1][aIr - 1][aItheta];//bottomright
		theValArray[3] = aValMat[aQ - 1][aIr - 1][aItheta - 1];//bottomleft
	}

	return theValArray;
}


// This function takes the array outputted by Stencil_for_Bilinear_Interp as the stencil
// and performs bi-linear interpolation at the specified (r, theta) location
double  CGrid_2D::Bilinear_Interp(const array<double, 4>& aValArray, const int aIr, const int aItheta,
	const double  aRloc, const double  aThetaloc)
{
	double  r_pos = aIr * fDr + fr0; // radial position
	double  theta_pos = aItheta * fDtheta + ftheta0; // theta position
	double  return_value = aValArray[0] * ((aRloc - (r_pos - fDr)) / fDr) * ((aThetaloc - (theta_pos - fDtheta)) / fDtheta) // top right
		+ aValArray[1] * ((aRloc - (r_pos - fDr)) / fDr) * ((theta_pos - aThetaloc) / fDtheta) // top left
		+ aValArray[2] * ((r_pos - aRloc) / fDr) * ((aThetaloc - (theta_pos - fDtheta)) / fDtheta) // bottom right
		+ aValArray[3] * ((r_pos - aRloc) / fDr) * ((theta_pos - aThetaloc) / fDtheta); // top left
	return return_value;
}


// This function returns an array of 2 double  that stores the value at
// left and right gridpoints.
// Order of Dimensions for aValMat: mode, slice, radius, theta
array<double, 2> CGrid_2D::Stencil_for_Linear_Interp_Theta(const array_3D& aValMat, const int aQ,
	const int aIr, const int aItheta) {
	array<double, 2> theValArray;
	// checking if the interpolation is between ftheta-1 and ftheta,
	// if so, we will replace the value at ftheta to be the value at theta=0
	if (aItheta == ftheta) {
		theValArray[0] = aValMat[aQ - 1][aIr][aItheta - 1]; // left
		theValArray[1] = aValMat[aQ - 1][aIr][0]; //value at theta = 0 (instead of theta = 2pi)
	}
	else if (aItheta == 0) {
		theValArray[0] = aValMat[aQ - 1][aIr][ftheta - 1]; // left
		theValArray[1] = aValMat[aQ - 1][aIr][aItheta]; // right
	}
	else {
		theValArray[0] = aValMat[aQ - 1][aIr][aItheta - 1]; // left
		theValArray[1] = aValMat[aQ - 1][aIr][aItheta]; // right
	}

	return theValArray;
}


// This function takes the array outputted by Stencil_for_Linear_Interp_Theta as the stencil
// and performs interpolation at the specified theta location
double  CGrid_2D::Linear_Interp_Theta(const array<double, 2>& aValArray, const int aIndex, const double  aThetaloc) {
	return aValArray[0] + ((aValArray[1] - aValArray[0]) / fDtheta) * (aThetaloc - (ftheta0 + (double(aIndex) - 1.0) * fDtheta));
}


// This function constructs an array with length 2 for Linear Interpolation
// and returns array of length 2 containing the value of the two gridpoints for linear interpolation
array<double, 2> CGrid_2D::Stencil_for_Linear_Interp(const array_1D& aArray, const int aI) {
	array<double, 2> theValArray;
	if (aI == 0) {
		theValArray[0] = aArray[aI]; // left
		theValArray[1] = aArray[aI]; // right
	}
	else {
		theValArray[0] = aArray[aI - 1]; // left
		theValArray[1] = aArray[aI]; // right
	}

	return theValArray;
}

// Linear Interpolation in 1D
// aValStencil is the output array from Stencil_for_Linear_Interp
double  CGrid_2D::Linear_Interp(const array<double, 2>& aValStencil, const array<double, 2>& aLocStencil, const double  aXloc) {
	//return ((fx0 + ((double )aIndex)*fDx - aXState) * aValArray[0] + (aXState - (fx0 + ((double )aIndex - 1) * fDx)) * aValArray[1])/fDx;
	double  inc = aLocStencil[1] - aLocStencil[0];
	return aValStencil[0] + ((aValStencil[1] - aValStencil[0]) / inc) * (aXloc - aLocStencil[0]);
}