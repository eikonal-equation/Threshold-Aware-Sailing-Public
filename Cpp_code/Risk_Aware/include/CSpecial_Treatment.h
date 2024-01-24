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
  * File: CSpecial_Treatment.h
  *
  * Author: MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu
  *
  * Description: This file contains the declarations of the class "CSpecial_Treatment"
  * and its member functions that implement the "Deadline-upgrade" policy if the value function
  * is sufficiently close to 1, and assign the risk-neutral policy to our risk-aware policy
  * if the value function is sufficiently close to 0.
  *
  *
  * Details of all of these functions are found in CSpecial_Treatment.cpp.
  *
  *============================================================================*/

#pragma once
#ifndef SPECIAL_TREATMENT_H
#define SPECIAL_TREATMENT_H

//----------------------Project specific header files---------------------------
#include "CGrid_2D.h"
//---------------------------Libraries-----------------------------------------
#include <cmath>
#include <array>
#include <boost/multi_array.hpp>
#include <boost/bimap.hpp>

//---------------------------Definitions---------------------------------------
typedef boost::bimap<int, double> boostMap;
typedef boostMap::value_type mapEntry;
typedef boost::multi_array<double, 3> array_3D;
typedef boost::multi_array<double, 1> array_1D;


class CSpecial_Treatment {// Class declaration begins here
public:
	CSpecial_Treatment() = default;
	CSpecial_Treatment(int aNum_slice_in_mem, const int aR, const int aTheta, double aTol_one, double aTol_zero) {
		fNum_slices = aNum_slice_in_mem; //Number of s-slices kept in the RAM
		fR = aR; // number of points on spatial grid = fR + 1
		fTheta = aTheta;  // number of theta points = fTheta + 1
		fTol_1 = aTol_one; // how close is close to 1?
		fTol_zero = aTol_zero; // how close is close to 0?
	}

	//This function assigns the risk-neutral policy to our risk-aware policy if the value function is exactly equal to 0
	void Assign_if_Exactly_0(const int aMode, const int aI, const int aJ, const int aK,
		const array_3D& aStat_policy, const array_3D& aStat_switch,
		array_4D& aVal_mat, array_4D_short& aSwitch_mat, array_4D& aPolicy_mat, const boostMap aEnum_map_for_steer_angle,
		array_4D_short& aVal_mat_out, array_4D_short& aPolicy_mat_out);

	//This function assigns the risk-neutral policy to our risk-aware policy if the value function is sufficiently close to 0
	void Assign_if_Near_0(const int aI, const int aJ, const int aK, const int aMode, const double aValfunc_new,
		array_4D& aVal_mat, array_4D& aPolicy_mat, array_4D_short& aSwitch_mat,
		const array_3D& aStat_policy, const array_3D& aStat_switch, const boostMap aEnum_map_for_steer_angle,
		array_4D_short& aVal_mat_out, array_4D_short& aPolicy_mat_out);

	//This function checks if the value function at slice i is close to 1. If so, assign "w = 1" to all slices > i
	void Assign_if_Near_1(const int aI, const int aActual_SliceNum, const int aJ, const int aK, const int aMode,
		const double aValfunc_new, const double aCtrl_new,
		const bool aSwitch_flag, array_3D& aCheck_mat, array_4D& aVal_mat, array_4D& aPolicy_mat, array_4D_short& aSwitch_mat,
		array_4D_short& aVal_mat_out, array_4D_short& aPolicy_mat_out);

	//This function implements the "Deadline-upgrade" policy when the value function computed is sufficiently close to 1.
	//That is, we just need to assign the value function / optimal policy from the previous s-slice to the current s-slice.
	void Deadline_upgrade_policy(bool aIf_reached_top, const int aStart, const int aEnd,
		const int aI, const int aJ, const int aK, array_4D& aVal_mat, array_4D& aPolicy_mat, array_4D_short& aSwitch_mat,
		array_4D_short& aVal_mat_out, array_4D_short& aPolicy_mat_out);

	//============================Methods for if we are forced to switch===============================
	//This function checks if the value function for the case if we are forced to switch initially at slice i is close to 1. If so, assign "prob = 1" to all slices > i
	void Assign_if_Near_1_force_to_switch(const int aI, const int aActual_SliceNum, const int aJ, const int aK, const int aMode, const double aW_switch,
	 array_3D& aForce_to_switch_Check_mat, array_4D& aWswitch_mat, array_4D_short& aWswitch_mat_out);

	//This function assigns the risk-neutral policy to our risk-aware policy if the value function for the case if we are forced to switch initially is sufficiently close to 0
	void Assign_if_Near_0_force_to_switch(const int aI, const int aJ, const int aK, const int aMode, const double aW_switch, array_4D& aWswitch_mat, array_4D_short& aWswitch_mat_out);

//member variables of the class
private:
	int fNum_slices; //Number of s-slices kept in the RAM
	int fR; // number of points on spatial grid = fR + 1
	int fTheta; // number of theta points = fTheta + 1
	double fTol_1; // how close is close to 1?
	double fTol_zero;// how close is close to 0?
};

#endif // !SPECIAL_TREATMENT_H
