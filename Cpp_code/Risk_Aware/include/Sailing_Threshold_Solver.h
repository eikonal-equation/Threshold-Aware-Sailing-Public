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
  * File: Sailing_Threshold_Solver.h
  *
  * Author: MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu
  *
  * Description: This file contains the declarations of the class "Sailing_Threshold"
  * and its member functions that load speed profile and the corresponding risk-neutral
  * policy; implement the binary search; carry out the initializations and the main solver
  * that constructed based on a first-order semi-Lagrangian discretization with an "s-marching"
  * algorithm.
  *
  *
  * Details of all of these functions are found in Sailing_Threshold_Solver.cpp.
  *
  *============================================================================*/


#pragma once
#ifndef SAILING_THRESHOLD_SOLVER_H
#define SAILING_THRESHOLD_SOLVER_H
//----------------------Project specific header files---------------------------
#include "CGrid_2D.h"
#include "SimpleFunctions.h"
#include "CSpecial_Treatment.h"
#include "CMarch_Switch.h"

//-------------------------Libraires-------------------------------------------
#include<algorithm>
#include<numeric>
#include<tuple>
#include<boost/bimap.hpp>

//---------------------------Definitions---------------------------------------
enum class speed_method { INTERPOLATE, CEIL, SETCONTROL };

//---------------------------Main Class----------------------------------------
class Sailing_Threshold : public CMarch_Switch
{
public:

	// constructor(s)
	Sailing_Threshold() = default;
	Sailing_Threshold(int aRefinement_factor, int aSlice_factor, double aBudget,
		double aDiffConst, double aCost_to_Switch, double aTarget_radius,
		double aDrift, int aNumCtrls, string aSpeed_File, speed_method aSpeed_choice,
		string aFilepath_RN_policy, string aFilepath_RN_switch, int aWind_speed_index, double aMax_speed, bool aReturn_Wswitch,
		string aInterp_Choice, int aInterp_Order_of_ENO, double aTol_for_prob_near_1, double aTol_for_prob_near_0, int aNum_GH_pts, bool aIf_use_GoldenSS, double aTol_GoldenSS);

	//This function implements the binary search algorithm
	int binary_search(const vector<double> aVector, const double aXloc);

	// Loads in speed profile data
	void load_speed_data(string aSpeed_File);

	// Loads in risk-neutral policy
	void Read_Risk_Neutral_Policy();

	// Initialization of matrices which store value function values, optimal steering angles, and switchgrids
	void InitializeMat();
	// Initialization of control and speed profile
	void Initialize_Control_and_Speed_List(boostMap& aControl_List, array_1D& aSpeed_List, CGrid_2D& aGrid);

	//This function updates the value function / optimal policy (steering direction + switching grid) after computations of each "s"-slice
	void UpdateSlices();

	// This is the main solver that uses both time marching and value iterations. I.e., given a
	// maximum budget s, this function will compute the value function and the optimal policy for each slice of threshold
	// up to the maximum budget.
	void MainSolver();

	//This is the inner solver pointer where it points to which mode(s) are needed to be solved for efficiency
	tuple<int, int, int, int> InnerSolver_pointer(const int aI, const int aJ, const int aK, const array_3D& aCheck_mat);

	// This is the inner solver embedded in the main solver which only deals with the (r, theta, ctrl)-loop
	void InnerSolver(const int aStart, const int aEnd, const int aI, const int aJ, const int aK,
		const boostMap& control_list, const array_1D& speed_list, CGrid_2D& aGrid, CENO& aENOgrid,
		CSpecial_Treatment& aSpecial_treatment, const string aInterp_Choice, const bool aReturn_Wswitch,
		array<double, 2>& aW_update, array<double, 2>& aCtrl_update, array<double, 2>& aW_switch);

	// This is the Golden Section Search algorithm embedded in the Marching step (if decided to use)
	void Golden_section_search_in_Marching(const int aStart, const int aEnd, const int aIndx_current, const int aJ, const int aK,
		const boostMap& control_list, const array_1D& speed_list, CGrid_2D& aGrid, CENO& aENOgrid,
		const string aInterp_Choice, array<double, 2>& aW_update, array<double, 2>& aCtrl_update);


	//============================Methods for if we are forced to switch================================
	void InnerSolver_force_to_switch(const int aStart, const int aEnd, const int aI, const int aJ, const int aK, const array_3D& aCheck_mat,
		CSpecial_Treatment aSpecialTreatmentClass, CGrid_2D& aGrid, CENO& aENOgrid);

	/** Grid parameters */
	double getThetaGridSize() const; // number of theta points
	double getThetaMax() const; // largest point at theta-axis
	double getRadiusGridSize() const; // // number of r points on spatial grid
	double getRadiusMax() const; // largest point at r-axis
	double getGFactor() const; // for grid refinement studies
	double getTargetSize() const; // radius of target region

	// Writing Domain Parameters to file
	void writeDomainToFile(const std::string aFilename) const;

	~Sailing_Threshold() {};

//member variables of the class
private:
	int fNumCtrls; //discretization of steering control variable
	string fInterp_Choice; //choice of Interpolation method
	vector<vector<double>> fSpeed_Data; // speed profile data
	vector<double> fAngles; //list of base control values (steering angles)
	array_1D fSpeed_column; //list of base speed values
	string fSpeed_File; //filename of the speed profile
	speed_method fSpeed_Method; //method of interpreting/discretizing the speed profile
	int fWind_Speed_Index; //wind speed index in the .pol file

	double fMax_speed; //max speed of the sailboat
	array_3D fMyCheck; //3D array to check if "prob = 1"
	array_4D fMyValfn; //4D array for the value function
	array_4D fMyPolicy; //4D array for the optimal steering angles
	array_4D_short fMySwitch; //4D array for the switchgrid (data type: "short")

	array_4D fMyValfn_force_to_switch; //4D array for the value function if we are forced to switch initially
	array_3D fMyCheck_force_to_switch;  //3D array to check if "prob = 1" if we are forced to switch initially

	array_4D_short fMyPolicy_out; //4D array for the policy (output only; output data type: "short")
	array_4D_short fMyValfn_out; //4D array for the value function (output only; output data type: "short")
	array_4D_short fMyValfn_force_to_switch_out; //4D array for the value function if we are forced to switch intially (output only; output data type: "short")

	string fRN_policy_path; //filepath for the risk-neutral optimal steering angles
	string fRN_switch_path; //filepath for the risk-neutral switchgrid
	array_3D fRisk_Neutral_Policy; //3D array to store the risk-neutral optimal steering angles
	array_3D fRisk_Neutral_Switch; //3D array to store the risk-neutral switchgrid

	bool fReturn_Wswitch; //boolean variable to decide if we want to return the value function for the case if we are forced to switch initially
	int fInterp_Order; //interpolation order
	bool fIf_use_GoldenSS; //boolean variable to decied if we want to GSS or not
	double fDu; //Delta u (steering angle variable)
};


inline double Sailing_Threshold::getThetaGridSize() const {
	return fTheta;
}

inline double Sailing_Threshold::getThetaMax() const {
	return fThetamax;
}

inline double Sailing_Threshold::getRadiusGridSize() const {
	return 	fR;
}

inline double Sailing_Threshold::getRadiusMax() const {
	return fRmax;
}

inline double Sailing_Threshold::getGFactor() const {
	return fRefinement_factor;
}

inline double Sailing_Threshold::getTargetSize() const {
	return fR_target;
}

#endif // !SAILING_THRESHOLD_SOLVER_H
