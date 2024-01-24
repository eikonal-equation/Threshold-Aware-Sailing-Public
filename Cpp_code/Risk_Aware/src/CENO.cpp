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
  * File: CENO3.cpp
  *
  * Author: MingYi Wang
  *
  * Description: This file contains acutal implementations of the member functions
  * that compute Newton 2nd and 3rd divided difference; interpolation in Newton form;
  * 1D quadratic/cubic ENO interpolation with and without a periodic boundary
  * condition; and 2D bi-quadratic/cubic ENO interpolation with a periodic
  * boundary condition in the "theta" variable.
  *
  *============================================================================*/

//----------------------Project specific header files---------------------------
#include "CENO.h"


// This function computes the second divided difference
// D2(input) : the second 0-th order divided difference
// D12(input) : the first 1st order divided difference
// aValueArray(input) : 3x1 vector of corresponding values of sample points
//
// D3(output) : the third 0-th order divided difference
// D23(output) : the second 1st order divided difference
// D123(output) : the first 2nd order divided difference
//
template <typename ValType, size_t Dim>
array<ValType, Dim>
CENO::Second_Divided_Difference(const double D2, const double D12, const array<ValType, Dim>& aValueArray, const array<ValType, Dim>& aStencilArray)
{
	// D123 : = f[x1, x2, x3]
	//        = (f[x2, x3] - f[x1, x2]) / (x3 - x1).
	array<ValType, Dim> Return_array{};
	Return_array[0] = aValueArray[2]; // D3
	Return_array[1] = (Return_array[0] - D2) / (aStencilArray[2] - aStencilArray[1]); // D23
	Return_array[2] = (Return_array[1] - D12) / (aStencilArray[2] - aStencilArray[0]); // D123

	return Return_array;
}


// This function computes the third divided difference
// D3(input) : the third 0-th order divided difference
// D23(input) : the second 1st order divided difference
// D123(input) : the first 2nd order divided difference
// aValueArray(input) : 4x1 vector of corresponding values of sample points
//
// D1234(output) : the first 3rd order divided difference
//
template <typename ValType, size_t Dim>
array<ValType, Dim> CENO::Third_Divided_Difference(const array<ValType, Dim> aDDarr, const array<ValType, Dim>& aValueArray, const array<ValType, Dim>& aStencilArray)
{
	// D1234 : = f[x1, x2, x3, x4]
	//         = (f[x2, x3, x4] - f[x1, x2, x3]) / (x4 - x1).
	array<ValType, Dim> Return_array{};
	Return_array[0] = aValueArray[3]; // D4
	Return_array[1] = (Return_array[0] - aDDarr[0]) / (aStencilArray[3] - aStencilArray[2]); // D34
	Return_array[2] = (Return_array[1] - aDDarr[1]) / (aStencilArray[3] - aStencilArray[1]); // D234
	Return_array[3] = (Return_array[2] - aDDarr[2]) / (aStencilArray[3] - aStencilArray[0]); // D1234
	return Return_array;
}

// This function computes cubic interpolation in Newton form
// aStencilArray(input) : 4x1 vector of sample points
// aCoefficientArray(input) : 4x1 coefficient vector of x computed from Newton divided difference.
// xloc(input) : coordinate of the query point
// val(output) : interpolated value
template <typename ValType, size_t Dim>
double CENO::NewtonInterp(const array<ValType, Dim>& aStencilArray, const array<ValType, Dim>& aCoefficientArray, const double aXloc)
{
	const int n = aCoefficientArray.size() - 1; //degree of the interpolating polynomial
	double val = aCoefficientArray[n];
	for (int i = 0; i < n; i++)
	{
		val = aCoefficientArray[n - i - 1] + (aXloc - aStencilArray[n - i - 1]) * val;
	}
	return val;
}


//Adding left/right points to a given stencil
template <typename ValType, size_t Dim>
void CENO::Construct_ENO_stencil(const int aIndex, const double ax0, const double aDx, const int aOrderOfInterp, const array<ValType, 2 * (Dim - 1)>& aValueArray,
	array<ValType, Dim>& aLeft_stencil, array<ValType, Dim>& aRight_stencil, array<ValType, Dim>& aLeft_value_arr, array<ValType, Dim>& aRight_value_arr)
{
	// for given aIndex, set the original two-point stencil to be [x(aIndex - 1), x(aIndex)]
	//  since the grid is uniform, and fx0 = x(0), then x(k) = fx0 + kIndex * fDx
	//  Here as we passed the whole 4-point stencil in, the value function corresponding to
	//  the first two-point stencil is fixed to be [aValueArray(1), aValueArray(2)]

	aLeft_stencil[0] = ax0 + (double(aIndex) - 1) * aDx;
	aLeft_stencil[1] = ax0 + double(aIndex) * aDx;
	aRight_stencil = aLeft_stencil;

	aLeft_value_arr[0] = aValueArray[aOrderOfInterp - 1];
	aLeft_value_arr[1] = aValueArray[aOrderOfInterp];
	aRight_value_arr = aLeft_value_arr;

	// add one point either from the left or from the right to form left and
	// right stencil respectively. Note here we just concatenate them because
	//	the order does not matter for divided differences
	aLeft_stencil[2] = ax0 + (double(aIndex) - 2) * aDx;
	aRight_stencil[2] = ax0 + (double(aIndex) + 1) * aDx;
	aLeft_value_arr[2] = aValueArray[aOrderOfInterp - 2];
	aRight_value_arr[2] = aValueArray[aOrderOfInterp + 1];

}


//First computations (1st and 2nd Divided Differences) needed for ENO
template <typename ValType, size_t Dim>
void CENO::FirstStep_for_ENO(const array<ValType, Dim>& aLeft_stencil, const array<ValType, Dim>& aRight_stencil, const array<ValType, Dim>& aLeft_value_arr, const array<ValType, Dim>& aRight_value_arr, const double aDx,
	array<ValType, Dim>& aDD_diag, array<ValType, Dim>& aDD_left_arr, array<ValType, Dim>& aDD_right_arr)
{
	// compute first divided difference
	double D1 = aLeft_value_arr[0];
	double D2 = aLeft_value_arr[1];
	double D12 = (D2 - D1) / aDx;
	aDD_diag[0] = D1;
	aDD_diag[1] = D12;

	// compute the second divided differences
	aDD_left_arr = Second_Divided_Difference(D2, D12, aLeft_value_arr, aLeft_stencil);
	aDD_right_arr = Second_Divided_Difference(D2, D12, aRight_value_arr, aRight_stencil);
}

//Choosing stencil
template <typename ValType, size_t Dim>
void CENO::Choosing_stencil(const double DDleft, const double DDright, bool& LeftFlag, array<ValType, Dim>& aDD_diag, const int aCurrent_order, int& aLeft_count, int& aRight_count,
	array<ValType, Dim>& aLeft_stencil, array<ValType, Dim>& aRight_stencil, array<ValType, Dim>& aLeft_value_arr, array<ValType, Dim>& aRight_value_arr)
{
	if (abs(DDleft) <= abs(DDright))
	{
		// if the seond divided difference of the left stencil is less than the one of the right, then choose the left stencil.
		aRight_stencil = aLeft_stencil;
		aRight_value_arr = aLeft_value_arr;
		LeftFlag = true;
		aDD_diag[aCurrent_order] = DDleft;
		aLeft_count++;
	}
	else // Otherwise, choose the right stencil
	{
		aLeft_stencil = aRight_stencil;
		aLeft_value_arr = aRight_value_arr;
		LeftFlag = false;
		aDD_diag[aCurrent_order] = DDright;
		aRight_count++;
	}

}


//Adding left/right points to the current stencil (or order higher than 2)
template <typename ValType, size_t Dim>
void CENO::Adding_points_to_ENOstencil(const bool LeftFlag, const int aLeft_count, const int aRight_count, const int aCurrent_order, const int aInterp_order, const int aIndex, const double ax0, const double aDx,
	const array<ValType, 2 * (Dim - 1)>& aValueArray, const array<ValType, Dim>& aDD_left_arr, const array<ValType, Dim>& aDD_right_arr,
	array<ValType, Dim>& aDD_temp_arr, array<ValType, Dim>& aLeft_stencil, array<ValType, Dim>& aRight_stencil, array<ValType, Dim>& aLeft_value_arr, array<ValType, Dim>& aRight_value_arr)
{
	aLeft_stencil[aCurrent_order] = ax0 + (double(aIndex) - double(aLeft_count) - 2) * aDx;
	aRight_stencil[aCurrent_order] = ax0 + (double(aIndex) + double(aRight_count) + 1) * aDx;
	aLeft_value_arr[aCurrent_order] = aValueArray[aInterp_order - aLeft_count - 2];
	aRight_value_arr[aCurrent_order] = aValueArray[aInterp_order + aRight_count + 1];

	if (LeftFlag)
	{
		aDD_temp_arr = aDD_left_arr;
	}
	else
	{
		aDD_temp_arr = aDD_right_arr;
	}
}


// This function computes a more efficient version of ENO quadratic interpolation on a uniform grid
// aValueArray(input) : a 4x1 vector of values of the 6-point stencil to be used for ENO3 interpolation
// xloc(input) : coordinate of the query point
// kIndex(input) : the first index such that xloc <= x, where x is the vector of sample points
//
// val(output) : interpolated value
//
double CENO::ENO_quad_interp_Rdim(const array<double, 4>& aValueArray, const int aRIndex, const double aRloc, const double aDr, const double ar0, const int aRDim)
{
	double val = 0; // initialize the interpolated value to be 0
	// initialize arrays with length 4 as storage
	array<double, 3> left_stencil;
	array<double, 3> right_stencil;
	// values at left stencil and right stencil respectively
	array<double, 3> yleft;
	array<double, 3> yright;
	array<double, 3> DD_diag;
	array <double, 3> DDleft, DDright;

	if (aRIndex == 1)
	{
		array<double, 3> stencil{ ar0, ar0 + aDr, ar0 + 2 * aDr };
		array<double, 3> temp_value_arrary{ aValueArray[0],aValueArray[1],aValueArray[2] };
		val = NewtonInterp(stencil, temp_value_arrary, aRloc);
	}
	else
	{	//Construct the very first 4-pt ENO stencil to be used (Interp Order = 2)
		Construct_ENO_stencil(aRIndex, ar0, aDr, 2, aValueArray, left_stencil, right_stencil, yleft, yright);
		//Compute the 1st and 2nd DD
		FirstStep_for_ENO(left_stencil, right_stencil, yleft, yright, aDr, DD_diag, DDleft, DDright);

		if (abs(DDleft[2]) <= abs(DDright[2])) // choosing stencil for the quadractic order and compute the interpolated value
		{
			DD_diag[2] = DDleft[2];
			val = NewtonInterp(left_stencil, DD_diag, aRloc);

		}
		else // Otherwise, choose the right stencil
		{
			DD_diag[2] = DDright[2];
			val = NewtonInterp(right_stencil, DD_diag, aRloc);
		}
	}
	return val;
}

// Main ENO qudaratic interploation function in 1D (perodic in Theta)
double CENO::ENO_quad_interp_Thetadim(const array<double, 4>& aValueArray, const int aThetaIndex, const double aThetaloc, const double aDtheta, const double atheta0, const int aThetaDim)
{
	double val = 0; // initialize the interpolated value to be 0
	// initialize arrays with length 4 as storage
	array<double, 3> left_stencil;
	array<double, 3> right_stencil;
	// values at left stencil and right stencil respectively
	array<double, 3> yleft;
	array<double, 3> yright;
	array<double, 3> DD_diag;
	array <double, 3> DDleft, DDright;


	//Construct the very first 4-pt ENO stencil to be used (Interp Order = 2)
	Construct_ENO_stencil(aThetaIndex, atheta0, aDtheta, 2, aValueArray, left_stencil, right_stencil, yleft, yright);
	//Compute the 1st and 2nd DD
	FirstStep_for_ENO(left_stencil, right_stencil, yleft, yright, aDtheta, DD_diag, DDleft, DDright);

	if (abs(DDleft[2]) <= abs(DDright[2])) // choosing stencil for the quadractic order and compute the interpolated value
	{
		DD_diag[2] = DDleft[2];
		val = NewtonInterp(left_stencil, DD_diag, aThetaloc);

	}
	else // Otherwise, choose the right stencil
	{
		DD_diag[2] = DDright[2];
		val = NewtonInterp(right_stencil, DD_diag, aThetaloc);
	}
	return val;

}



// This function computes a more efficient version of ENO cubic interpolation on a uniform grid
// aValueArray(input) : a 6x1 vector of values of the 6-point stencil to be used for ENO3 interpolation
// xloc(input) : coordinate of the query point
// kIndex(input) : the first index such that xloc <= x, where x is the vector of sample points
//
// val(output) : interpolated value
//
double CENO::ENO_cubic_interp_Rdim(const array<double, 6>& aValueArray, const int aRIndex, const double aRloc, const double aDr, const double ar0, const int aRDim)
{
	double val = 0; // initialize the interpolated value to be 0
	bool Leftflag = true;

	// initialize arrays with length 4 as storage
	array<double, 4> left_stencil;
	array<double, 4> right_stencil;
	// values at left stencil and right stencil respectively
	array<double, 4> yleft;
	array<double, 4> yright;
	array<double, 4> DD_diag;
	array <double, 4> DDleft, DDright, DD_temp;
	int left_count = 0;
	int right_count = 0;

	if (aRIndex <= 2)
	{
		// if kIndex = 1 or kIndex = 2, we simply choose the first 4-point stencil and use pchip interpolation
		/*vector<double> stencil{ fx0,fx0 + fDx,fx0 + f2h,fx0 + f3h };
		vector<double> temp_value_arrary{ aValueArray[0],aValueArray[1],aValueArray[2],aValueArray[3] };*/
		/*pchip<vector<double>> interp(std::move(stencil), std::move(temp_value_arrary));
		val = interp(aXloc);*/
		array<double, 4> stencil{ ar0,ar0 + aDr,ar0 + 2 * aDr, ar0 + 3 * aDr };
		array<double, 4> temp_value_arrary{ aValueArray[0],aValueArray[1],aValueArray[2],aValueArray[3] };
		val = NewtonInterp(stencil, temp_value_arrary, aRloc);
	}
	else
	{   //Construct the very first 4-pt ENO stencil to be used (Interp Order = 3)
		Construct_ENO_stencil(aRIndex, ar0, aDr, 3, aValueArray, left_stencil, right_stencil, yleft, yright);
		//Compute the 1st and 2nd DD
		FirstStep_for_ENO(left_stencil, right_stencil, yleft, yright, aDr, DD_diag, DDleft, DDright);
		//Choose the left/right stencil based on 2nd DD (current order = 2)
		Choosing_stencil(DDleft[2], DDright[2], Leftflag, DD_diag, 2, left_count, right_count, left_stencil, right_stencil, yleft, yright);
		//Add a point to the right and a point to the left (current order = 3)
		Adding_points_to_ENOstencil(Leftflag, left_count, right_count, 3, 3, aRIndex, ar0, aDr, aValueArray, DDleft, DDright, DD_temp,
			left_stencil, right_stencil, yleft, yright);

		// compute the third divivded differences
		DDleft = Third_Divided_Difference(DD_temp, yleft, left_stencil);
		DDright = Third_Divided_Difference(DD_temp, yright, right_stencil);

		if (abs(DDleft[3]) <= abs(DDright[3]))// choosing stencil for the cubic order and compute the interpolated value
		{
			DD_diag[3] = DDleft[3];
			val = NewtonInterp(left_stencil, DD_diag, aRloc); // compute the interpolated value
		}
		else
		{
			DD_diag[3] = DDright[3];
			val = NewtonInterp(right_stencil, DD_diag, aRloc); // compute the interpolated value
		}
	}
	return val;
}

// Main ENO cubic interploation function in 1D (perodic in Theta)
double CENO::ENO_cubic_interp_Thetadim(const array<double, 6>& aValueArray, const int aThetaIndex, const double aThetaloc, const double aDtheta, const double atheta0, const int aThetaDim)
{
	double val = 0; // initialize the interpolated value to be 0
	bool Leftflag = true;

	// initialize arrays with length 4 as storage
	array<double, 4> left_stencil;
	array<double, 4> right_stencil;
	// values at left stencil and right stencil respectively
	array<double, 4> yleft;
	array<double, 4> yright;
	array<double, 4> DD_diag;
	array <double, 4> DDleft, DDright, DD_temp;
	int left_count = 0;
	int right_count = 0;


	//Construct the very first 4-pt ENO stencil to be used (Interp Order = 3)
	Construct_ENO_stencil(aThetaIndex, atheta0, aDtheta, 3, aValueArray, left_stencil, right_stencil, yleft, yright);
	//Compute the 1st and 2nd DD
	FirstStep_for_ENO(left_stencil, right_stencil, yleft, yright, aDtheta, DD_diag, DDleft, DDright);
	//Choose the left/right stencil based on 2nd DD (current order = 2)
	Choosing_stencil(DDleft[2], DDright[2], Leftflag, DD_diag, 2, left_count, right_count, left_stencil, right_stencil, yleft, yright);
	//Add a point to the right and a point to the left (current order = 3)
	Adding_points_to_ENOstencil(Leftflag, left_count, right_count, 3, 3, aThetaIndex, atheta0, aDtheta, aValueArray, DDleft, DDright, DD_temp,
		left_stencil, right_stencil, yleft, yright);

	// compute the third divivded differences
	DDleft = Third_Divided_Difference(DD_temp, yleft, left_stencil);
	DDright = Third_Divided_Difference(DD_temp, yright, right_stencil);

	if (abs(DDleft[3]) <= abs(DDright[3]))// choosing stencil for the cubic order and compute the interpolated value
	{
		DD_diag[3] = DDleft[3];
		val = NewtonInterp(left_stencil, DD_diag, aThetaloc); // compute the interpolated value
	}
	else
	{
		DD_diag[3] = DDright[3];
		val = NewtonInterp(right_stencil, DD_diag, aThetaloc); // compute the interpolated value
	}

	return val;

}

//Construct the very first ENO stencil in 2D
template <typename ValType, size_t Dim>
void CENO::Construct_ENO_stencil_2D_quad(const int aRIndex, const int aThetaIndex, const double aThetaloc, const int aOrderOfInterp, const array_2D& aValueMatrix,
	array<ValType, 2 * (Dim - 1)>& aLeft_col, array<ValType, 2 * (Dim - 1)>& aRight_col,
	array<ValType, Dim>& aLeft_stencil, array<ValType, Dim>& aRight_stencil, array<ValType, Dim>& aLeft_value_arr, array<ValType, Dim>& aRight_value_arr, double& aTemp_val_left, double& aTemp_val_right)
{
	aLeft_stencil[0] = fr0 + (double(aRIndex) - 1) * fDr;
	aLeft_stencil[1] = fr0 + double(aRIndex) * fDr;
	aRight_stencil = aLeft_stencil;

	for (int k = 0; k < 2 * aOrderOfInterp; k++) {
		aLeft_col[k] = aValueMatrix[aOrderOfInterp - 1][k];
		aRight_col[k] = aValueMatrix[aOrderOfInterp][k];
	}

	aLeft_value_arr[0] = ENO_quad_interp_Thetadim(aLeft_col, aThetaIndex, aThetaloc, fDtheta, ftheta0, ftheta);
	aLeft_value_arr[1] = ENO_quad_interp_Thetadim(aRight_col, aThetaIndex, aThetaloc, fDtheta, ftheta0, ftheta);

	aRight_value_arr = aLeft_value_arr;


	aLeft_stencil[2] = fr0 + (double(aRIndex) - 2) * fDr;
	aRight_stencil[2] = fr0 + (double(aRIndex) + 1) * fDr;

	for (int k = 0; k < 2 * aOrderOfInterp; k++) {
		aLeft_col[k] = aValueMatrix[aOrderOfInterp - 2][k];
		aRight_col[k] = aValueMatrix[aOrderOfInterp + 1][k];
	}


	aLeft_value_arr[2] = ENO_quad_interp_Thetadim(aLeft_col, aThetaIndex, aThetaloc, fDtheta, ftheta0, ftheta);
	aRight_value_arr[2] = ENO_quad_interp_Thetadim(aRight_col, aThetaIndex, aThetaloc, fDtheta, ftheta0, ftheta);

	aTemp_val_left = aLeft_value_arr[2];
	aTemp_val_right = aRight_value_arr[2];

}



//Construct the very first ENO stencil in 2D
template <typename ValType, size_t Dim>
void CENO::Construct_ENO_stencil_2D_cubic(const int aRIndex, const int aThetaIndex, const double aThetaloc, const int aOrderOfInterp, const array_2D& aValueMatrix,
	array<ValType, 2 * (Dim - 1)>& aLeft_col, array<ValType, 2 * (Dim - 1)>& aRight_col,
	array<ValType, Dim>& aLeft_stencil, array<ValType, Dim>& aRight_stencil, array<ValType, Dim>& aLeft_value_arr, array<ValType, Dim>& aRight_value_arr, double& aTemp_val_left, double& aTemp_val_right)
{
	aLeft_stencil[0] = fr0 + (double(aRIndex) - 1) * fDr;
	aLeft_stencil[1] = fr0 + double(aRIndex) * fDr;
	aRight_stencil = aLeft_stencil;

	for (int k = 0; k < 2 * aOrderOfInterp; k++) {
		aLeft_col[k] = aValueMatrix[aOrderOfInterp - 1][k];
		aRight_col[k] = aValueMatrix[aOrderOfInterp][k];
	}



	aLeft_value_arr[0] = ENO_cubic_interp_Thetadim(aLeft_col, aThetaIndex, aThetaloc, fDtheta, ftheta0, ftheta);
	aLeft_value_arr[1] = ENO_cubic_interp_Thetadim(aRight_col, aThetaIndex, aThetaloc, fDtheta, ftheta0, ftheta);

	aRight_value_arr = aLeft_value_arr;


	aLeft_stencil[2] = fr0 + (double(aRIndex) - 2) * fDr;
	aRight_stencil[2] = fr0 + (double(aRIndex) + 1) * fDr;

	for (int k = 0; k < 2 * aOrderOfInterp; k++) {
		aLeft_col[k] = aValueMatrix[aOrderOfInterp - 2][k];
		aRight_col[k] = aValueMatrix[aOrderOfInterp + 1][k];
	}


	aLeft_value_arr[2] = ENO_cubic_interp_Thetadim(aLeft_col, aThetaIndex, aThetaloc, fDtheta, ftheta0, ftheta);
	aRight_value_arr[2] = ENO_cubic_interp_Thetadim(aRight_col, aThetaIndex, aThetaloc, fDtheta, ftheta0, ftheta);

	aTemp_val_left = aLeft_value_arr[2];
	aTemp_val_right = aRight_value_arr[2];

}



//Adding left/right points to the current stencil (for order higher than 2)
template <typename ValType, size_t Dim>
void CENO::Adding_points_to_ENOstencil_2D(const bool LeftFlag, const int aLeft_count, const int aRight_count, const int aCurrent_order, const int aInterp_order, const int aRIndex, const int aThetaIndex, const double aThetaloc,
	const array_2D& aValueMatrix, const array<ValType, Dim>& aDD_left_arr, const array<ValType, Dim>& aDD_right_arr, array<ValType, Dim>& aDD_temp_arr,
	array<ValType, 2 * (Dim - 1)>& aLeft_col, array<ValType, 2 * (Dim - 1)>& aRight_col, array<ValType, Dim>& aLeft_stencil, array<ValType, Dim>& aRight_stencil, array<ValType, Dim>& aLeft_value_arr, array<ValType, Dim>& aRight_value_arr,
	const double aTemp_val_left, const double aTemp_val_right)
{
	aLeft_stencil[aCurrent_order] = fr0 + (double(aRIndex) - double(aLeft_count) - 2) * fDr;
	aRight_stencil[aCurrent_order] = fr0 + (double(aRIndex) + double(aRight_count) + 1) * fDr;

	if (LeftFlag)
	{
		aDD_temp_arr = aDD_left_arr;
		// copy the column at (aRIndex - 3) for cubic interpolation
		for (int k = 0; k < 6; k++)
		{
			aLeft_col[k] = aValueMatrix[aInterp_order - aLeft_count - 2][k];
		}
		aLeft_value_arr[aCurrent_order] = ENO_cubic_interp_Thetadim(aLeft_col, aThetaIndex, aThetaloc, fDtheta, ftheta0, ftheta);
		aRight_value_arr[aCurrent_order] = aTemp_val_right;
	}
	else
	{
		aDD_temp_arr = aDD_right_arr;
		// copy the column at (aRIndex + 2) for cubic interpolation
		for (int k = 0; k < 6; k++)
		{
			aRight_col[k] = aValueMatrix[aInterp_order + aRight_count + 1][k];
		}
		aLeft_value_arr[aCurrent_order] = aTemp_val_left;
		aRight_value_arr[aCurrent_order] = ENO_cubic_interp_Thetadim(aRight_col, aThetaIndex, aThetaloc, fDtheta, ftheta0, ftheta);
	}
}

// Main ENO bi-quadractic interpolation function in 2D (perodic in Theta, Theta direction first then R direction)
double CENO::ENO_quad_interp_RTheta(const array_2D& aValueMatrix, const int aRIndex, const int aThetaIndex, const double aRloc, const double aThetaloc)
{
	double val = 0; // initialize the interpolated value to be 0
	// initialize arrays with length 4 as storage
	array<double, 3> left_stencil;
	array<double, 3> right_stencil;
	// values at left stencil and right stencil respectively
	array<double, 3> yleft;
	array<double, 3> yright;
	array<double, 3> DD_diag;
	array <double, 3> DDleft, DDright;
	array<double, 4> left_col, right_col;
	double yl_temp, yr_temp;

	if (aRIndex == 1) {
		array<double, 3> left_stencil{ fr0, fr0 + fDr, fr0 + 2 * fDr };
		for (int ii = 0; ii < 3; ii++) {
			for (int k = 0; k < 4; k++) {
				left_col[k] = aValueMatrix[ii][k];
			}
			yleft[ii] = ENO_quad_interp_Thetadim(left_col, aThetaIndex, aThetaloc, fDtheta, ftheta0, ftheta);
		}
		val = NewtonInterp(left_stencil, yleft, aRloc);
	}
	else //interior
	{
		//Construct the very first 4-pt ENO stencil (in the R-direction) to be used (Interp Order = 2)
		Construct_ENO_stencil_2D_quad(aRIndex, aThetaIndex, aThetaloc, 2, aValueMatrix, left_col, right_col, left_stencil, right_stencil, yleft, yright, yl_temp, yr_temp);
		//Compute the 1st and 2nd DD
		FirstStep_for_ENO(left_stencil, right_stencil, yleft, yright, fDr, DD_diag, DDleft, DDright);

		if (abs(DDleft[2]) <= abs(DDright[2])) // choosing stencil for the quadractic order and compute the interpolated value
		{
			DD_diag[2] = DDleft[2];
			val = NewtonInterp(left_stencil, DD_diag, aRloc);

		}
		else // Otherwise, choose the right stencil
		{
			DD_diag[2] = DDright[2];
			val = NewtonInterp(right_stencil, DD_diag, aRloc);
		}
	}
	return val;
}


// Main ENO bi-cubic interpolation function in 2D
double CENO::ENO_cubic_interp_RTheta(const array_2D& aValueMatrix, const int aRIndex, const int aThetaIndex, const double aRloc, const double aThetaloc)
{
	double val = 0; // initialize the interpolated value to be 0
	bool Leftflag = true;

	// initialize arrays with length 4 as storage
	array<double, 4> left_stencil;
	array<double, 4> right_stencil;
	// values at left stencil and right stencil respectively
	array<double, 4> yleft;
	array<double, 4> yright;
	array<double, 4> DD_diag;
	array <double, 4> DDleft, DDright, DD_temp;
	int left_count = 0;
	int right_count = 0;
	array<double, 6> left_col, right_col;
	double yl_temp, yr_temp;

	if (aRIndex <= 2)
	{
		// if kIndex = 1 or kIndex = 2, we simply choose the first 4-point stencil and use pchip interpolation
		array<double, 4> left_stencil{ fr0,fr0 + fDr,fr0 + 2 * fDr, fr0 + 3 * fDr };
		for (int ii = 0; ii < 4; ii++) {
			for (int k = 0; k < 6; k++) {
				left_col[k] = aValueMatrix[ii][k];
			}
			yleft[ii] = ENO_cubic_interp_Thetadim(left_col, aThetaIndex, aThetaloc, fDtheta, ftheta0, ftheta);
		}
		val = NewtonInterp(left_stencil, yleft, aRloc);
	}
	else
	{   //Construct the very first 4-pt ENO stencil to be used (Interp Order = 3)
		Construct_ENO_stencil_2D_cubic(aRIndex, aThetaIndex, aThetaloc, 3, aValueMatrix, left_col, right_col, left_stencil, right_stencil, yleft, yright, yl_temp, yr_temp);
		//Compute the 1st and 2nd DD
		FirstStep_for_ENO(left_stencil, right_stencil, yleft, yright, fDr, DD_diag, DDleft, DDright);
		//Choose the left/right stencil based on 2nd DD (current order = 2)
		Choosing_stencil(DDleft[2], DDright[2], Leftflag, DD_diag, 2, left_count, right_count, left_stencil, right_stencil, yleft, yright);
		//Add a point to the right and a point to the left (current order = 3)
		Adding_points_to_ENOstencil_2D(Leftflag, left_count, right_count, 3, 3, aRIndex, aThetaIndex, aThetaloc, aValueMatrix, DDleft, DDright, DD_temp, left_col, right_col,
			left_stencil, right_stencil, yleft, yright, yl_temp, yr_temp);

		// compute the third divivded differences
		DDleft = Third_Divided_Difference(DD_temp, yleft, left_stencil);
		DDright = Third_Divided_Difference(DD_temp, yright, right_stencil);

		if (abs(DDleft[3]) <= abs(DDright[3]))// choosing stencil for the cubic order and compute the interpolated value
		{
			DD_diag[3] = DDleft[3];
			val = NewtonInterp(left_stencil, DD_diag, aRloc); // compute the interpolated value
		}
		else
		{
			DD_diag[3] = DDright[3];
			val = NewtonInterp(right_stencil, DD_diag, aRloc); // compute the interpolated value
		}
	}
	return val;
}


// This function constructs the matrix needed for ENO interpolation in 2D
array_2D CENO::Matrix_for_ENO_Interp_RTheta(const array_4D& aValMat, const int aInterpOrder, const int aMode, const int aSliceNum, const int aRIndex, const int aThetaIndex)
{
	array_2D aValueMatrix(boost::extents[2 * aInterpOrder][2 * aInterpOrder]);
	int End_indx_R = aInterpOrder + 1;
	if (aRIndex <= aInterpOrder - 1) {
		//Inner boundary of R (could take care of it as the next case, but since our target radius = 0.1 by default, it is fine as of now.)
		Matrix_for_ENO_Interp_RTheta_subcases(0, End_indx_R, aInterpOrder, aValMat, aValueMatrix, aMode, aSliceNum, aRIndex, aThetaIndex);
	}
	else if (aRIndex >= fr - aInterpOrder + 2) {
		//Outerboundary of R: Need to extend boundary conditions
		int aDiff_count = fr - aRIndex;
		End_indx_R += aDiff_count;
		Matrix_for_ENO_Interp_RTheta_subcases(fr - aInterpOrder - aDiff_count, End_indx_R, aInterpOrder, aValMat,
			aValueMatrix, aMode, aSliceNum, aRIndex, aThetaIndex);
		//Adding ghost points outside Rmax
		for (int j = End_indx_R; j < 2 * aInterpOrder; j++) {
			for (int k = 0; k < 2 * aInterpOrder; k++) {
				aValueMatrix[j][k] = 0.0;
			}
		}
	}
	else { //Interior
		End_indx_R = 2 * aInterpOrder;
		Matrix_for_ENO_Interp_RTheta_subcases(aRIndex - aInterpOrder, End_indx_R, aInterpOrder, aValMat,
			aValueMatrix, aMode, aSliceNum, aRIndex, aThetaIndex);
	}
	return aValueMatrix;
}

//Subfunction to build the matrix for bi-ENO interpolation
void CENO::Matrix_for_ENO_Interp_RTheta_subcases(const int aSpeicalIndex_for_R, const int aEndindex_R, const int aInterpOrder, const array_4D& aValMat, array_2D& aValueMatrix,
	const int aMode, const int aSliceNum, const int aRIndex, const int aThetaIndex)
{
	if (aThetaIndex <= aInterpOrder - 1) {
		int aGhostNum = aInterpOrder - aThetaIndex;
		for (int j = 0; j < aEndindex_R; j++) {
			for (int kk = aGhostNum; kk > 0; kk--) {
				aValueMatrix[j][aGhostNum - kk] = aValMat[aMode - 1][aSliceNum - 1][j + aSpeicalIndex_for_R][ftheta - kk];
			}
		}
		for (int j = 0; j < aEndindex_R; j++) {
			for (int k = aGhostNum; k < 2 * aInterpOrder; k++) {
				aValueMatrix[j][k] = aValMat[aMode - 1][aSliceNum - 1][j + aSpeicalIndex_for_R][k - aGhostNum];
			}
		}
	}
	else if (aThetaIndex >= ftheta - aInterpOrder + 2) {
		int aGhostNum = aInterpOrder - 1 - (ftheta - aThetaIndex);
		for (int j = 0; j < aEndindex_R; j++) {
			for (int k = 0; k < 2 * aInterpOrder - aGhostNum - 1; k++) {
				aValueMatrix[j][k] = aValMat[aMode - 1][aSliceNum - 1][j + aSpeicalIndex_for_R][ftheta - 2 * aInterpOrder + 1 + aGhostNum + k];
			}
		}
		for (int j = 0; j < aEndindex_R; j++) {
			for (int k = 2 * aInterpOrder - aGhostNum - 1; k < 2 * aInterpOrder; k++) {
				aValueMatrix[j][k] = aValMat[aMode - 1][aSliceNum - 1][j + aSpeicalIndex_for_R][k + aGhostNum + 1 - 2 * aInterpOrder];
			}
		}
	}
	else {
		for (int j = 0; j < aEndindex_R; j++) {
			for (int k = 0; k < 2 * aInterpOrder; k++) {
				aValueMatrix[j][k] = aValMat[aMode - 1][aSliceNum - 1][j + aSpeicalIndex_for_R][k + aThetaIndex - aInterpOrder];
			}
		}
	}
}
