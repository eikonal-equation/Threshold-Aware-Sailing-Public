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
  * File: CENO3.h
  *
  * Author: MingYi Wang
  *
  * Description: This file contains the declarations of the class "CENO" and
  * its member functions that compute
  * Newton 2nd and 3rd divided difference; interpolation in Newton form;
  * 1D quadratic/cubic ENO interpolation with and without a periodic boundary
  * condition; and 2D bi-quadratic/cubic ENO interpolation with a periodic
  * boundary condition in the "theta" variable.

  *
  * Details of all of these functions are found in CENO3.cpp.
  *
  *============================================================================*/

#pragma once
#ifndef CENO_H
#define CENO_H

//----------------------Project specific header files---------------------------
#include "CGrid_2D.h"

class CENO:public CGrid_2D // Class declaration begins here
{
public:
	// constructors
	CENO() = default;
	CENO(double aDr, double aRmin, int aRDim, double aDtheta, double aThetamin, int aThetaDim, double aDs) : CGrid_2D(aDr,aRmin,aRDim,aDtheta,aThetamin,aThetaDim,aDs) {};
	~CENO() {};

	//---------------------------------ENO-Interpolation--------------------------------------------------------
	// Main ENO quadratic interploation function in 1D (non-periodic normal case in R)
	double ENO_quad_interp_Rdim(const array<double, 4>& aValueArray, const int aRIndex, const double aRloc, const double aDr, const double ar0, const int aRDim);

	// Main ENO quadratic interploation function in 1D (perodic in Theta)
	double ENO_quad_interp_Thetadim(const array<double, 4>& aValueArray, const int aThetaIndex, const double aThetaloc, const double aDtheta, const double atheta0, const int aThetaDim);

	// Main ENO bi-quadratic interpolation function in 2D (perodic in Theta)
	double ENO_quad_interp_RTheta(const array_2D& aValueMatrix, const int aRIndex, const int aThetaIndex, const double aRloc, const double aThetaloc);

	// Main ENO cubic interploation function in 1D (non-periodic normal case in R)
	double ENO_cubic_interp_Rdim(const array<double, 6>& aValueArray, const int aRIndex, const double aRloc, const double aDr, const double ar0, const int aRDim);

	// Main ENO cubic interploation function in 1D (perodic in Theta)
	double ENO_cubic_interp_Thetadim(const array<double, 6>& aValueArray, const int aThetaIndex, const double aThetaloc, const double aDtheta, const double atheta0, const int aThetaDim);

	// Main ENO bi-cubic interpolation function in 2D (perodic in Theta)
	double ENO_cubic_interp_RTheta(const array_2D& aValueMatrix, const int aRIndex, const int aThetaIndex, const double aRloc, const double aThetaloc);

	// This function computes the second divided difference
	template <typename ValType, size_t Dim>
	array<ValType,Dim> Second_Divided_Difference(const double D2, const double D12, const array<ValType, Dim>& aValueArray, const array<ValType, Dim>& aStencilArray);

	// This function computes the third divided difference
	template <typename ValType, size_t Dim>
	array<ValType, Dim> Third_Divided_Difference(const array<ValType,Dim> aDDarr, const array<ValType, Dim>& aValueArray, const array<ValType, Dim>& aStencilArray);

	// This function computes interpolation in Newton form
	template <typename ValType, size_t Dim>
	double NewtonInterp(const array<ValType, Dim>& aStencilArray, const array<ValType, Dim>& aCoefficientArray, const double aXloc);


	//====================================================== SUBFUNCTIONS FOR 1D ===========================================================================================
	//Construct the very first ENO stencil in 1D
	template <typename ValType, size_t Dim>
	void Construct_ENO_stencil(const int aIndex, const double ax0, const double aDx, const int aOrderOfInterp, const array<ValType, 2*(Dim - 1)>& aValueArray,
		array<ValType, Dim>& aLeft_stencil, array<ValType, Dim>& aRight_stencil, array<ValType, Dim>& aLeft_value_arr, array<ValType, Dim>& aRight_value_arr);

	//First computations (1st and 2nd Divided Differences) needed for ENO
	template <typename ValType, size_t Dim>
	void FirstStep_for_ENO(const array<ValType, Dim>& aLeft_stencil, const array<ValType, Dim>& aRight_stencil, const array<ValType, Dim>& aLeft_value, const array<ValType, Dim>& aRight_value, const double aDx,
		array<ValType, Dim>& aDD_diag, array<ValType, Dim>& aDD_left_arr, array<ValType, Dim>& aDD_right_arr);

	//Choosing stencil
	template <typename ValType, size_t Dim>
	void Choosing_stencil(const double DDleft, const double DDright, bool& LeftFlag, array<ValType, Dim>& aDD_diag, const int aCurrent_order, int& aLeft_count, int& aRight_count,
		array<ValType, Dim>& aLeft_stencil, array<ValType, Dim>& aRight_stencil, array<ValType, Dim>& aLeft_value_arr, array<ValType, Dim>& aRight_value_arr);

	//Adding left/right points to the current stencil (for order higher than 2)
	template <typename ValType, size_t Dim>
	void Adding_points_to_ENOstencil(const bool LeftFlag, const int aLeft_count, const int aRight_count, const int aCurrent_order, const int aInterp_order, const int aIndex, const double ax0, const double aDx,
		const array<ValType, 2 * (Dim - 1)>& aValueArray, const array<ValType, Dim>& aDD_left_arr, const array<ValType, Dim>& aDD_right_arr,
		array<ValType, Dim>& aDD_temp_arr, array<ValType, Dim>& aLeft_stencil, array<ValType, Dim>& aRight_stencil, array<ValType, Dim>& aLeft_value_arr, array<ValType, Dim>& aRight_value_arr);


	//====================================================== SUBFUNCTIONS FOR 2D ===========================================================================================
	//Construct the very first ENO stencil in 2D
	template <typename ValType, size_t Dim>
	void Construct_ENO_stencil_2D_quad(const int aRIndex, const int aThetaIndex, const double aThetaloc, const int aOrderOfInterp, const array_2D& aValueMatrix,
		array<ValType, 2 * (Dim-1)>& aLeft_col, array<ValType, 2 * (Dim-1)>& aRight_col,
		array<ValType, Dim>& aLeft_stencil, array<ValType, Dim>& aRight_stencil, array<ValType, Dim>& aLeft_value_arr, array<ValType, Dim>& aRight_value_arr, double& aTemp_val_left, double& aTemp_val_right);


	//Construct the very first ENO stencil in 2D
	template <typename ValType, size_t Dim>
	void Construct_ENO_stencil_2D_cubic(const int aRIndex, const int aThetaIndex, const double aThetaloc, const int aOrderOfInterp, const array_2D& aValueMatrix,
		array<ValType, 2 * (Dim - 1)>& aLeft_col, array<ValType, 2 * (Dim - 1)>& aRight_col,
		array<ValType, Dim>& aLeft_stencil, array<ValType, Dim>& aRight_stencil, array<ValType, Dim>& aLeft_value_arr, array<ValType, Dim>& aRight_value_arr, double& aTemp_val_left, double& aTemp_val_right);


	//Adding left/right points to the current stencil (for order higher than 2)
	template <typename ValType, size_t Dim>
	void Adding_points_to_ENOstencil_2D(const bool LeftFlag, const int aLeft_count, const int aRight_count, const int aCurrent_order, const int aInterp_order, const int aRIndex, const int aThetaIndex, const double aThetaloc,
		const array_2D& aValueMatrix, const array<ValType, Dim>& aDD_left_arr, const array<ValType, Dim>& aDD_right_arr, array<ValType, Dim>& aDD_temp_arr,
		array<ValType, 2*(Dim-1)>& aLeft_col, array<ValType, 2 * (Dim - 1)>& aRight_col, array<ValType, Dim>& aLeft_stencil, array<ValType, Dim>& aRight_stencil, array<ValType, Dim>& aLeft_value_arr, array<ValType, Dim>& aRight_value_arr,
		const double aTemp_val_left, const double aTemp_val_right);

	// This function constructs the matrix needed for ENO interpolation in 2D
	array_2D Matrix_for_ENO_Interp_RTheta(const array_4D& aValMat, const int aInterpOrder, const int aMode, const int aSliceNum, const int aRIndex, const int aThetaIndex);

	//Subfunction to build the matrix for bi-ENO interpolation
	void Matrix_for_ENO_Interp_RTheta_subcases(const int aSpeicalIndex_for_R, const int aEndindex_R, const int aInterpOrder, const array_4D& aValmat, array_2D& aValueMatrix,
		const int aMode, const int aSliceNum, const int aRIndex, const int aThetaIndex);

	//This function constructs the array needed for ENO interpolation in 1D (in the theta direction only; i.e., it is periodic)
	template <int Dim>
inline	array<double, Dim> Array_for_ENO_Interp_Theta(const array_4D& aValMat, const int aInterpOrder, const int aMode, const int aSliceNum, const int aThetaIndex, const int aJ)
{
	array<double, Dim> aValueArray;

	//Taking care of the 2*pi periodicity here
	if (aThetaIndex <= aInterpOrder - 1)
	{ //left boundary of [0, 2*pi]
		int aGhostNum = aInterpOrder - aThetaIndex;

		for (int kk = aGhostNum; kk > 0; kk--) {
			aValueArray[aGhostNum - kk] = aValMat[aMode - 1][aSliceNum][aJ][ftheta - kk];
		}

		for (int k = aGhostNum; k < 2 * aInterpOrder; k++) {
			aValueArray[k] = aValMat[aMode - 1][aSliceNum][aJ][k - aGhostNum];
		}
	}
	else if (aThetaIndex >= ftheta - aInterpOrder + 2)
	{//right boundary of [0, 2*pi]
		int aGhostNum = aInterpOrder - 1 - (ftheta - aThetaIndex);

		for (int k = 0; k < 2 * aInterpOrder - aGhostNum - 1; k++) {
			aValueArray[k] = aValMat[aMode - 1][aSliceNum][aJ][ftheta - 2 * aInterpOrder + 1 + aGhostNum + k];
		}

		for (int k = 2 * aInterpOrder - aGhostNum - 1; k < 2 * aInterpOrder; k++) {
			aValueArray[k] = aValMat[aMode - 1][aSliceNum][aJ][k + aGhostNum + 1 - 2 * aInterpOrder];
		}
	}
	else {//interior
		for (int k = 0; k < 2 * aInterpOrder; k++) {
			aValueArray[k] = aValMat[aMode - 1][aSliceNum][aJ][k + aThetaIndex - aInterpOrder];
		}
	}
	return aValueArray;
}

private:
};
#endif // !CENO_H

