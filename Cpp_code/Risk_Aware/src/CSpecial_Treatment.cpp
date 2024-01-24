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
  * File: CSpecial_Treatment.cpp
  *
  * Author: MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu
  *
  * Description: This file contains the implementations of the member functions
  * that carry out the "Deadline-upgrade" policy if the value function
  * is sufficiently close to 1, and assign the risk-neutral policy to our risk-aware policy
  * if the value function is sufficiently close to 0.
  *
  *============================================================================*/

//---------------------------Libraries-----------------------------------------
#include "CSpecial_Treatment.h"
#include <fstream>
#include <iostream>


//Checks if the value function at slice i is close to 1. If so, assign "w = 1" to all slices > i
void CSpecial_Treatment::Assign_if_Near_1(const int aI, const int aActual_SliceNum, const int aJ, const int aK, const int aMode,
	const double aValfunc_new, const double aCtrl_new,
	const bool aSwitch_flag, array_3D& aCheck_mat, array_4D& aVal_mat, array_4D& aPolicy_mat, array_4D_short& aSwitch_mat,
	array_4D_short& aVal_mat_out, array_4D_short& aPolicy_mat_out)
{
	if (abs(aValfunc_new - 1.0) < 1.0*fTol_1) {
		//if it is the first time of reaching prob 1
		aCheck_mat[aMode][aJ][aK] = aActual_SliceNum;
		aVal_mat[aMode][aI][aJ][aK] = 1.0;

		//for storage
		aVal_mat_out[aMode][aI][aJ][aK] = short(1.0 * 1e4);

		if (aSwitch_flag) {
			aSwitch_mat[aMode][aI][aJ][aK] = 1; //switching mode
		}
		else {
			aSwitch_mat[aMode][aI][aJ][aK] = 0; //staying on the current mode
		}
        aPolicy_mat[aMode][aI][aJ][aK] = aCtrl_new; //record the "best" steering angle no matter switching or not
		//for storage
		aPolicy_mat_out[aMode][aI][aJ][aK] = (short)round(aCtrl_new*1e4);
	}
}

//Assign the risk-neutral policy when the prob is near 0
void CSpecial_Treatment::Assign_if_Near_0(const int aI, const int aJ, const int aK, const int aMode, const double aValfunc_new,
	array_4D& aVal_mat, array_4D& aPolicy_mat, array_4D_short& aSwitch_mat,
	const array_3D& aStat_policy, const array_3D& aStat_switch, const boostMap aEnum_map_for_steer_angle,
	array_4D_short& aVal_mat_out, array_4D_short& aPolicy_mat_out)
{
	if (abs(aValfunc_new - 0.0) < fTol_zero) {

		aVal_mat[aMode][aI][aJ][aK] = 0.0;

		//for storage
		aVal_mat_out[aMode][aI][aJ][aK] = 0;

		if (aStat_switch[aMode][aJ][aK] == 1) {
			aPolicy_mat[aMode][aI][aJ][aK] = -1000; //record "-1000" for the steering angle
			aSwitch_mat[aMode][aI][aJ][aK] = 2; //switching mode (risk-neutral case)
			//for storage
			aPolicy_mat_out[aMode][aI][aJ][aK] = -1000;
		}
		else {
			aPolicy_mat[aMode][aI][aJ][aK] = aStat_policy[aMode][aJ][aK]; //assign the best "risk-neutral" steering angle
			aSwitch_mat[aMode][aI][aJ][aK] = 3; //staying on the current mode (risk-neutral case)
			//for storage
			aPolicy_mat_out[aMode][aI][aJ][aK] = (short)round(aStat_policy[aMode][aJ][aK] * 1e4);
		}
	}
}


//Assign the risk-neutral policy when prob is exactly 0
void CSpecial_Treatment::Assign_if_Exactly_0(const int aMode, const int aI, const int aJ, const int aK,
	const array_3D& aStat_policy, const array_3D& aStat_switch,
	array_4D& aVal_mat, array_4D_short& aSwitch_mat, array_4D& aPolicy_mat, const boostMap aEnum_map_for_steer_angle,
	array_4D_short& aVal_mat_out, array_4D_short& aPolicy_mat_out)
{
	aVal_mat[aMode][aI][aJ][aK] = 0.0;

	//for storage
	aVal_mat_out[aMode][aI][aJ][aK] = 0;

	if (aStat_switch[aMode][aJ][aK] == 1) {
		aSwitch_mat[aMode][aI][aJ][aK] = 2;
		aPolicy_mat[aMode][aI][aJ][aK] = -1000;
		//for storage
		aPolicy_mat_out[aMode][aI][aJ][aK] = -1000;
	}
	else {
		aSwitch_mat[aMode][aI][aJ][aK] = 3;
		aPolicy_mat[aMode][aI][aJ][aK] = aEnum_map_for_steer_angle.right.at(aStat_policy[aMode][aJ][aK]);
		//for storage
		aPolicy_mat_out[aMode][aI][aJ][aK] = (short)round(aStat_policy[aMode][aJ][aK] * 1e4);
	}
}

//This function implements the "Deadline-upgrade" policy when the value function computed is sufficiently close to 1.
//That is, we just need to assign the value function / optimal policy from the previous s-slice to the current s-slice.
void CSpecial_Treatment::Deadline_upgrade_policy(bool aIf_reached_top, const int aStart, const int aEnd,
	const int aI, const int aJ, const int aK, array_4D& aVal_mat, array_4D& aPolicy_mat, array_4D_short& aSwitch_mat,
	array_4D_short& aVal_mat_out, array_4D_short& aPolicy_mat_out)
{
	if (!aIf_reached_top){
		//If we haven't reached the max number (fNum_slices) of s-slices allowed in the RAM
		//Just assign everything from the previous s-slice to the current s-slice
		for (int Current_mode = aStart; Current_mode < aEnd; Current_mode++) {
			aVal_mat[Current_mode][aI][aJ][aK] = aVal_mat[Current_mode][aI - 1][aJ][aK];
			aSwitch_mat[Current_mode][aI][aJ][aK] = aSwitch_mat[Current_mode][aI - 1][aJ][aK];
			aPolicy_mat[Current_mode][aI][aJ][aK] = aPolicy_mat[Current_mode][aI - 1][aJ][aK];
			//for storage
			aVal_mat_out[Current_mode][aI][aJ][aK] = aVal_mat_out[Current_mode][aI - 1][aJ][aK];
			aPolicy_mat_out[Current_mode][aI][aJ][aK] = aPolicy_mat_out[Current_mode][aI - 1][aJ][aK];
		}
	}
	else {
		//If we have reached the max number (fNum_slices) of s-slices allowed in the RAM
		//Then assign everything from the (fNum_slices-2)-th s-slice to the top (fNum_slices-1)-th s-slice
		for (int Current_mode = aStart; Current_mode < aEnd; Current_mode++) {
			aVal_mat[Current_mode][fNum_slices -1][aJ][aK] = aVal_mat[Current_mode][fNum_slices - 2][aJ][aK];
			aSwitch_mat[Current_mode][fNum_slices - 1][aJ][aK] = aSwitch_mat[Current_mode][fNum_slices - 2][aJ][aK];
			aPolicy_mat[Current_mode][fNum_slices - 1][aJ][aK] = aPolicy_mat[Current_mode][fNum_slices - 2][aJ][aK];
			//for storage
			aVal_mat_out[Current_mode][fNum_slices - 1][aJ][aK] = aVal_mat_out[Current_mode][fNum_slices - 2][aJ][aK];
			aPolicy_mat_out[Current_mode][fNum_slices - 1][aJ][aK] = aPolicy_mat_out[Current_mode][fNum_slices - 2][aJ][aK];
		}
	}
}



//============================Methods for if we are forced to switch===============================
//This function checks if the value function for the case if we are forced to switch initially at slice i is close to 1. If so, assign "prob = 1" to all slices > i
void CSpecial_Treatment::Assign_if_Near_1_force_to_switch(const int aI, const int aActual_SliceNum,
	const int aJ, const int aK, const int aMode, const double aW_switch,
	array_3D& aForce_to_switch_Check_mat, array_4D& aWswitch_mat, array_4D_short& aWswitch_mat_out)
{
	if (abs(aW_switch - 1.0) < fTol_1) {
		//if it is the first time of reaching prob 1
		aForce_to_switch_Check_mat[aMode][aJ][aK] = aActual_SliceNum;
		aWswitch_mat[aMode][aI][aJ][aK] = 1.0;
		//for storage
		aWswitch_mat_out[aMode][aI][aJ][aK] = short(1.0*1e4);
	}
}

//This function assigns the risk-neutral policy to our risk-aware policy if the value function for the case if we are forced to switch initially is sufficiently close to 0
void CSpecial_Treatment::Assign_if_Near_0_force_to_switch(const int aI, const int aJ, const int aK, const int aMode, const double aW_switch, array_4D& aWswitch_mat,
	array_4D_short& aWswitch_mat_out)
{
	if (abs(aW_switch - 0.0) < fTol_zero) {
		aWswitch_mat[aMode][aI][aJ][aK] = 0.0;
		// for storage
		aWswitch_mat_out[aMode][aI][aJ][aK] = 0;
	}
}
