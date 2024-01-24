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
  * Description: This file contains the actual implementations of the member functions
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
// Stencil order: (-) === loc === (+)
//
tuple<array<double, 4>, double, double> CGrid_2D::Stencil_for_Bilinear_Interp(const array_4D& aValMat, const int aMode, const int aI,
	const int aJ, const int aK, const Where_to_interp aPlace4Interp) {
	array<double, 4> theValStencil;
	double theYplus, theXplus;

	if (aPlace4Interp == Where_to_interp::marching) {
		//checking if the bilinear interpolation is between ftheta-1 and ftheta,
		// if so, we will replace the value at ftheta to be the value at theta=0
		if (aK == ftheta) {
			theValStencil[0] = aValMat[aMode - 1][aI - 1][aJ][0]; // value at theta = 0 (instead of theta = 2pi)
			theValStencil[1] = aValMat[aMode - 1][aI - 1][aJ][aK - 1];//topleft
			theValStencil[2] = aValMat[aMode - 1][aI - 1][aJ - 1][0]; //value at theta = 0 (instead of theta = 2pi)
			theValStencil[3] = aValMat[aMode - 1][aI - 1][aJ - 1][aK - 1];//bottomleft
		}
		else if (aK == 0) { // will this ever happen?
			theValStencil[0] = aValMat[aMode - 1][aI - 1][aJ][aK];
			theValStencil[1] = aValMat[aMode - 1][aI - 1][aJ][ftheta - 1];//topleft
			theValStencil[2] = aValMat[aMode - 1][aI - 1][aJ - 1][aK];
			theValStencil[3] = aValMat[aMode - 1][aI - 1][aJ - 1][ftheta - 1];//bottomleft
		}
		else {
			theValStencil[0] = aValMat[aMode - 1][aI - 1][aJ][aK];// topright
			theValStencil[1] = aValMat[aMode - 1][aI - 1][aJ][aK - 1];//topleft
			theValStencil[2] = aValMat[aMode - 1][aI - 1][aJ - 1][aK];//bottomright
			theValStencil[3] = aValMat[aMode - 1][aI - 1][aJ - 1][aK - 1];//bottomleft
		}
		theXplus = double(aK) * fDtheta + ftheta0; // theta position (+)
		theYplus = double(aJ) * fDr + fr0; // radial position (+)
	}
	else if (aPlace4Interp == Where_to_interp::switching)
	{
		if (aK == ftheta) {
			theValStencil[0] = aValMat[aMode - 1][aI][aJ][0]; // value at theta = 0 (instead of theta = 2pi)
			theValStencil[1] = aValMat[aMode - 1][aI][aJ][aK - 1];//topleft
			theValStencil[2] = aValMat[aMode - 1][aI - 1][aJ][0]; //value at theta = 0 (instead of theta = 2pi)
			theValStencil[3] = aValMat[aMode - 1][aI - 1][aJ][aK - 1];//bottomleft
		}
		else if (aK == 0) {
			theValStencil[0] = aValMat[aMode - 1][aI][aJ][aK];
			theValStencil[1] = aValMat[aMode - 1][aI][aJ][ftheta - 1];//topleft
			theValStencil[2] = aValMat[aMode - 1][aI - 1][aJ][aK];
			theValStencil[3] = aValMat[aMode - 1][aI - 1][aJ][ftheta - 1];//bottomleft
		}
		else {
			theValStencil[0] = aValMat[aMode - 1][aI][aJ][aK];// topright
			theValStencil[1] = aValMat[aMode - 1][aI][aJ][aK - 1];//topleft
			theValStencil[2] = aValMat[aMode - 1][aI - 1][aJ][aK];//bottomright
			theValStencil[3] = aValMat[aMode - 1][aI - 1][aJ][aK - 1];//bottomleft
		}
		theXplus = double(aK) * fDtheta + ftheta0; // theta position (+)
		theYplus = double(aI) * fDs; // threshold position (+)

	}
	return make_tuple(theValStencil,theXplus,theYplus);
}


// This function returns an array of 2 double that stores the value at
// left and right gridpoints.
// Order of Dimensions for aValMat: mode, slice, radius, theta
array<double, 2> CGrid_2D::Stencil_for_Linear_Interp_Theta(const array_4D& aValMat, const int aMode, const int aI,
	const int aJ, const int aK) {
	array<double, 2> theValArray;

	//checking if the interpolation is between ftheta-1 and ftheta,
	// if so, we will replace the value at ftheta to be the value at theta=0
	if (aK == ftheta) {
		theValArray[0] = aValMat[aMode - 1][aI][aJ][aK - 1]; // left
		theValArray[1] = aValMat[aMode - 1][aI][aJ][0]; //value at theta = 0 (instead of theta = 2pi)
	}
	else if (aK == 0) {
		theValArray[0] = aValMat[aMode - 1][aI][aJ][ftheta - 1]; // left
		theValArray[1] = aValMat[aMode - 1][aI][aJ][aK]; // right
	}
	else {
		theValArray[0] = aValMat[aMode - 1][aI][aJ][aK - 1]; // left
		theValArray[1] = aValMat[aMode - 1][aI][aJ][aK]; // right
	}

	return theValArray;
}

// This function construct an array with length 2 for Linear Interpolation
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