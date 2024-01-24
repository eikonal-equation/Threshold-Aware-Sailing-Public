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
  * File: Sailing_Threshold_Solver.cpp
  *
  * Author: MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu
  *
  * Description: This file contains the actual implementations of the member functions
  * that load speed profile and the corresponding risk-neutral
  * policy; implement the binary search; carry out the initializations and the main solver
  * that constructed based on a first-order semi-Lagrangian discretization with an "s-marching"
  * algorithm.
  *
  *============================================================================*/

//----------------------Project specific header files---------------------------
#include "Sailing_Threshold_Solver.h"
#include "WriteToFile.h"

//------------------------------Libraries------------------------------------
#include<chrono>
#include<omp.h>
#include<numeric>

//Define the constructor
Sailing_Threshold::Sailing_Threshold(int aRefinement_factor, int aSlice_factor, double aBudget,
	double aDiffConst, double aCost_to_Switch, double aTarget_radius,
	double aDrift, int aNumCtrls, string aSpeed_File, speed_method aSpeed_choice,
	string aFilepath_RN_policy, string aFilepath_RN_switch, int aWind_speed_index, double aMax_speed, bool aReturn_Wswitch,
	string aInterp_Choice, int aInterp_Order_of_ENO, double aTol_for_prob_near_1, double aTol_for_prob_near_0, int aNum_GH_pts, bool aIf_use_GoldenSS, double aTol_GoldenSS)
	: CMarch_Switch(aRefinement_factor,
		aSlice_factor, aBudget, aDiffConst, aCost_to_Switch, aTarget_radius, aDrift, aInterp_Order_of_ENO,
		aTol_for_prob_near_1, aTol_for_prob_near_0, aNum_GH_pts, aIf_use_GoldenSS, aTol_GoldenSS)
{

	fSpeed_File = aSpeed_File;
	fSpeed_Method = aSpeed_choice; // 3 choices (based on what is given in .pol file): (1) SETCONTROL sets controls as angles ,
								  //                                                   (2) INTERPOLATE takes given control angles and interpolates between angles
								 //                                                    (3) CEIL takes given control angles and uses next highest speed

	if (fSpeed_Method == speed_method::SETCONTROL) {
		fNumCtrls = 24; //number of core controls (fixed for all of our numerical experiments according to the data file)
	}
	else {
		fNumCtrls = aNumCtrls; //number of controls
	}

	fWind_Speed_Index = aWind_speed_index;
	fMax_speed = aMax_speed;
	fReturn_Wswitch = aReturn_Wswitch;

	fMyCheck.resize(boost::extents[2][fR + 1][fTheta + 1]);
	fMyValfn.resize(boost::extents[2][fNum_slices][fR + 1][fTheta + 1]);
	fMyPolicy.resize(boost::extents[2][fNum_slices][fR + 1][fTheta + 1]);
	fMySwitch.resize(boost::extents[2][fNum_slices][fR + 1][fTheta + 1]);


	fMyValfn_out.resize(boost::extents[2][fNum_slices][fR + 1][fTheta + 1]);
	fMyPolicy_out.resize(boost::extents[2][fNum_slices][fR + 1][fTheta + 1]);


	if (fReturn_Wswitch) {
		fMyValfn_force_to_switch.resize(boost::extents[2][fNum_slices][fR + 1][fTheta + 1]);
		fMyCheck_force_to_switch.resize(boost::extents[2][fR + 1][fTheta + 1]);

		fMyValfn_force_to_switch_out.resize(boost::extents[2][fNum_slices][fR + 1][fTheta + 1]);
	}
	else {
		fMyValfn_force_to_switch.resize(boost::extents[1][1][1][1]);
		fMyCheck_force_to_switch.resize(boost::extents[1][1][1]);

		fMyValfn_force_to_switch_out.resize(boost::extents[1][1][1][1]);
	}


	fRN_policy_path = aFilepath_RN_policy;
	fRN_switch_path = aFilepath_RN_switch;

	fRisk_Neutral_Policy.resize(boost::extents[2][fR + 1][fTheta + 1]);
	fRisk_Neutral_Switch.resize(boost::extents[2][fR + 1][fTheta + 1]);

	fInterp_Choice = aInterp_Choice;
	fInterp_Order = aInterp_Order_of_ENO;
	fIf_use_GoldenSS = aIf_use_GoldenSS;
	fDu = PI / double(fNumCtrls - 1);
}


/*
					1d 2d   interation
row				x  r     j
col				y  theta k
sheet			s  s     i
*/
// Initialization of matrices which store value function values, optimal steering angles, and switchgrids
void Sailing_Threshold::InitializeMat() {
	for (int m = 0; m < 2; m++) {
		for (int i = 0; i < fNum_slices; i++) {
			for (int j = 0; j < fR + 1; j++) {
				for (int k = 0; k < fTheta + 1; k++) {
					if (k == fTheta) {//checking for being nonsense place (theta = 2pi)
						fMyValfn[m][i][j][k] = -1000;
						fMyPolicy[m][i][j][k] = -1000;
						fMySwitch[m][i][j][k] = -1000;

						fMyValfn_out[m][i][j][k] = -1000;
						fMyPolicy_out[m][i][j][k] = -1000;

						if (fReturn_Wswitch) {
							fMyValfn_force_to_switch[m][i][j][k] = -1000;

							fMyValfn_force_to_switch_out[m][i][j][k] = -1000;
						}

					}
					else {
						//if 0<=r<=e0
						if (fR0 + j * fDr <= fR_target)
						{  // at 0 or 1 in state space ...
							fMyValfn[m][i][j][k] = 1.0; // if we are already at 0 or 1 in the state space, probability of exiting = 1
							fMyPolicy[m][i][j][k] = 0; // if we are already at 0 or 1, no need to select a control
							fMySwitch[m][i][j][k] = 0;

							fMyValfn_out[m][i][j][k] = short(1.0*1e4);
							fMyPolicy_out[m][i][j][k] = 0;

							if (fReturn_Wswitch) {
								fMyValfn_force_to_switch[m][i][j][k] = 1.0;

								fMyValfn_force_to_switch_out[m][i][j][k] = 1;
							}

						}
						else
						{
							fMyValfn[m][i][j][k] = 0; // initializing remaining probability values = 0
							fMyPolicy[m][i][j][k] = 0; // initializing remaining control values = 0
							fMySwitch[m][i][j][k] = 0; // initializing remaining mode switch values = 0

							fMyValfn_out[m][i][j][k] = 0;
							fMyPolicy_out[m][i][j][k] = 0;

							if (fReturn_Wswitch) {
								fMyValfn_force_to_switch[m][i][j][k] = 0;
								fMyValfn_force_to_switch_out[m][i][j][k] = 0;
							}

						}
					}
				}
			}
		}
	}

	for (int m = 0; m < 2; m++) // mode
	{
		for (int j = 0; j < fR + 1; j++) // radius
		{
			for (int k = 0; k < fTheta + 1; k++) { //theta
				fMyCheck[m][j][k] = double(fM) + 1;
			}
		}
	}
}

//This function updates the value function / optimal policy (steering direction + switching grid) after computations of each "s"-slice
void Sailing_Threshold::UpdateSlices() {
	for (int m = 0; m < 2; m++) {
		for (int i = 0; i < fNum_slices - 1; i++) {
			for (int j = 0; j < fR + 1; j++) {
				for (int k = 0; k < fTheta + 1; k++) {

					fMyValfn[m][i][j][k] = fMyValfn[m][i + 1][j][k];
					fMyPolicy[m][i][j][k] = fMyPolicy[m][i + 1][j][k];
					fMySwitch[m][i][j][k] = fMySwitch[m][i + 1][j][k];

					fMyValfn_out[m][i][j][k] = fMyValfn_out[m][i + 1][j][k];
					fMyPolicy_out[m][i][j][k] = fMyPolicy_out[m][i + 1][j][k];

					if (fReturn_Wswitch) {
						fMyValfn_force_to_switch[m][i][j][k] = fMyValfn_force_to_switch[m][i + 1][j][k];

						fMyValfn_force_to_switch_out[m][i][j][k] = fMyValfn_force_to_switch_out[m][i + 1][j][k];
					}

				}
			}
		}
	}
}

//This function implements the binary search algorithm
int Sailing_Threshold::binary_search(const vector<double> aVector, const double aXloc)
{
	//cout << "binary search" << endl;
	int left = 0;
	int right = aVector.size() - 1;
	while (right > left + 1) {
		int mid = (right + left) / 2;
		if (aVector.at(mid) < aXloc) {
			left = mid;
		}
		else {
			right = mid;
		}
	}
	return right;
}

// Loads in speed profile data
void Sailing_Threshold::load_speed_data(string filename) {
	//read in the speed file
	ifstream infile(filename);
	string line;

	int i = 0;

	while (getline(infile, line)) // Read a line
	{
		if (i != 0) {
			fSpeed_Data.push_back(vector<double>()); //Speed data holds data from .pol file
			vector<double>& row = fSpeed_Data.back();
			istringstream iss(line);
			double value;
			int j = 0;
			while (iss >> value) {
				if (j != 0) {
					row.push_back(value);
				}
				if (j == 0) {
					fAngles.push_back(PI * value / 180); //angles holds list of angles given in .pol file
				}
				j++;
			}

		}
		i++;
	}
	infile.close();
}

//Defines list of controls to iterate over and their respective speeds
void Sailing_Threshold::Initialize_Control_and_Speed_List(boostMap& aControl_List, array_1D& aSpeed_List, CGrid_2D& aGrid)
{
	if (fSpeed_Method == speed_method::SETCONTROL)
	{

		aSpeed_List.resize(boost::extents[fAngles.size()]);
		fSpeed_column.resize(boost::extents[fAngles.size()]);
		for (int n = 0; n < fNumCtrls; n++)
		{
			//aControl_List[n] = fAngles.at(n);
			aControl_List.insert(mapEntry(n + 1, fAngles.at(n)));
			aSpeed_List[n] = fSpeed_Data[n][fWind_Speed_Index];
			fSpeed_column[n] = aSpeed_List[n];
		}
	}
	else if (fSpeed_Method == speed_method::INTERPOLATE) //Interpolate between angles in the speed plot
	{

		aSpeed_List.resize(boost::extents[fNumCtrls]);
		fSpeed_column.resize(boost::extents[fAngles.size()]);
		for (int n = 0; n < fAngles.size(); n++)
		{
			fSpeed_column[n] = fSpeed_Data[n][fWind_Speed_Index];
		}
		for (int n = 0; n < fNumCtrls; n++)
		{
			double ctrl = n * PI / (double(fNumCtrls) - 1.0);
			//aControl_List[n] = ctrl;
			aControl_List.insert(mapEntry(n + 1, ctrl));
			int index = binary_search(fAngles, ctrl);
			array<double, 2> speed_stencil = aGrid.Stencil_for_Linear_Interp(fSpeed_column, index);
			array<double, 2> angle_stencil;
			assert(index > 0);
			angle_stencil[0] = fAngles.at(index - 1);
			angle_stencil[1] = fAngles.at(index);
			aSpeed_List[n] = aGrid.Linear_Interp(speed_stencil, angle_stencil, ctrl);
		}
	}
	else if (fSpeed_Method == speed_method::CEIL) //take the speed corresponding to first angle greater than control
	{
		fSpeed_column.resize(boost::extents[fAngles.size()]);
		for (int n = 0; n < fAngles.size(); n++)
		{
			fSpeed_column[n] = fSpeed_Data[n][fWind_Speed_Index];
		}

		aSpeed_List.resize(boost::extents[fNumCtrls]);
		for (int n = 0; n < fNumCtrls; n++)
		{
			double ctrl = n * PI / (double(fNumCtrls) - 1.0);
			aControl_List.insert(mapEntry(n + 1, ctrl));
			int index = binary_search(fAngles, ctrl);
			aSpeed_List[n] = fSpeed_Data[index][fWind_Speed_Index];
		}
	}

	double max_speed = 0;
	for (int n = 0; n < aSpeed_List.size(); n++) {
		if (aSpeed_List[n] > max_speed) {
			max_speed = aSpeed_List[n];
		}
	}

	double max_speed_base = 0;
	for (int n = 0; n < fAngles.size(); n++) {
		if (fSpeed_column[n] > max_speed_base) {
			max_speed_base = fSpeed_column[n];
		}
	}

	for (int n = 0; n < fNumCtrls; n++) {
		aSpeed_List[n] = aSpeed_List[n] * fMax_speed / max_speed;
	}

	for (int n = 0; n < fAngles.size(); n++)
	{
		fSpeed_column[n] = fSpeed_column[n] * fMax_speed / max_speed_base;
	}

	// Writing control values and corresponding speeds to file
	io::writeMapToFile1D("ControlList.dat", aControl_List);
	io::writeToFile1D<double>("SpeedList.dat", aSpeed_List);
	io::writeVectorToFile<double>("Base_ControlList.dat", fAngles);
	io::writeToFile1D<double>("Base_SpeedList.dat", fSpeed_column);
}

// Loads in risk-neutral policy
void Sailing_Threshold::Read_Risk_Neutral_Policy()
{

	streampos size;

	char* memblock = NULL;

	streampos size2;
	char* memblock2 = NULL;


	for (int i = 0; i < 2; i++) {
		ifstream file(fRN_policy_path + to_string(i + 1) + ".dat", ios::in | ios::binary | ios::ate);
		size = file.tellg();
		memblock = new char[size];
		file.seekg(0, ios::beg);
		file.read(memblock, size);
		file.close();


		ifstream file2(fRN_switch_path + to_string(i + 1) + ".dat", ios::in | ios::binary | ios::ate);
		size2 = file2.tellg();
		memblock2 = new char[size2];
		file2.seekg(0, ios::beg);
		file2.read(memblock2, size2);
		file2.close();

		double* double_values = (double*)memblock;
		double* double_values2 = (double*)memblock2;

		//Should add assert here to check that we've read in the correct number of values
		//assert(double_values.size() == (fR+1)*fTheta);

		for (int j = 0; j < fR + 1; j++)
		{
			for (int k = 0; k < fTheta + 1; k++) {
				if (k == fTheta) {//checking for being nonsense place (theta = 2pi)
					fRisk_Neutral_Policy[i][j][k] = -1000;
					fRisk_Neutral_Switch[i][j][k] = -1000;
				}
				else {
					fRisk_Neutral_Policy[i][j][k] = double_values[fTheta * j + k];
					fRisk_Neutral_Switch[i][j][k] = double_values2[fTheta * j + k];
				}
			}
		}


	}
	delete[] memblock;
	delete[] memblock2;
}

// This is the Golden Section Search algorithm embedded in the Marching step (if decided to use)
inline void Sailing_Threshold::Golden_section_search_in_Marching(const int aStart, const int aEnd, const int aIndx_current, const int aJ, const int aK,
	const boostMap& control_list, const array_1D& speed_list, CGrid_2D& aGrid, CENO& aENOgrid,
	const string aInterp_Choice, array<double, 2>& aW_update, array<double, 2>& aCtrl_update)
{
	for (int Current_mode = aStart; Current_mode < aEnd; Current_mode++) {
		//Finding interval on which to do golden section search
		array<double, 2> lower_upper_ctrl_vals{ -1,-1 };
		array<double, 2> lower_upper_speed_vals{ -1,-1 };

		//Only interpolate using the core data points from the speed profile
		if (aCtrl_update[Current_mode] == control_list.left.at(1)) {
			for (int mm = 0; mm < 2; mm++) {
				lower_upper_ctrl_vals[mm] = fAngles[mm];
				lower_upper_speed_vals[mm] = fSpeed_column[mm];

			}
		}
		else if (aCtrl_update[Current_mode] == control_list.left.at(fNumCtrls)) {
			for (int mm = 0; mm < 2; mm++) {
				lower_upper_ctrl_vals[mm] = fAngles[fAngles.size() - 2 + mm];
				lower_upper_speed_vals[mm] = fSpeed_column[fAngles.size() - 2 + mm];

			}
		}
		else {
			int myIndex = binary_search(fAngles, aCtrl_update[Current_mode]);
			for (int mm = 0; mm < 2; mm++) {
				lower_upper_ctrl_vals[mm] = fAngles[myIndex - 1 + mm];
				lower_upper_speed_vals[mm] = fSpeed_column[myIndex - 1 + mm];

			}
		}


		//Actual implementation of the GSS algorithm
		if (aInterp_Choice == "ENO") {
			std::tie(aW_update[Current_mode], aCtrl_update[Current_mode]) =
				GoldenSection(lower_upper_ctrl_vals, lower_upper_speed_vals, aIndx_current, aJ, aK, Current_mode + 1, aENOgrid, aGrid, fInterp_Order, fMyValfn);
		}
		else if (aInterp_Choice == "Bilinear") {
			std::tie(aW_update[Current_mode], aCtrl_update[Current_mode]) =
				GoldenSection(lower_upper_ctrl_vals, lower_upper_speed_vals, aIndx_current, aJ, aK, Current_mode + 1, aGrid, fMyValfn);
		}
	}
}

//This is the inner solver pointer where it points to which mode(s) are needed to be solved for efficiency
inline tuple<int, int, int, int> Sailing_Threshold::InnerSolver_pointer(const int aI, const int aJ, const int aK, const array_3D& aCheck_mat)
{
	int aStart = -10;
	int aEnd = -10;
	int aStart_ddlupg = -10;
	int aEnd_ddlupg = -10;

	if ((aI < aCheck_mat[0][aJ][aK]) && (aI < aCheck_mat[1][aJ][aK]))
	{// neither reaches prob 1, neither upgrades the deadline
		aStart = 0;
		aEnd = 2;
		aStart_ddlupg = 0;
		aEnd_ddlupg = 0;
	}
	//if only mode 1 value function @ (j,k) has not yet reached 1 // mode 2 deadline-upgrade
	else if (aI < aCheck_mat[0][aJ][aK])
	{
		aStart = 0;
		aEnd = 1;
		aStart_ddlupg = 1;
		aEnd_ddlupg = 2;

	}
	//if only mode 2 value function @ (j,k) has not yet reached 1 // mode 1 deadline-upgrade
	else if (aI < aCheck_mat[1][aJ][aK])
	{
		aStart = 1;
		aEnd = 2;
		aStart_ddlupg = 0;
		aEnd_ddlupg = 1;
	}
	else { //both have reached prob 1, both need to upgrade the deadline
		aStart = 0;
		aEnd = 0;
		aStart_ddlupg = 0;
		aEnd_ddlupg = 2;
	}
	return make_tuple(aStart, aEnd, aStart_ddlupg, aEnd_ddlupg);
}

// Overloaded with ENO class
inline void Sailing_Threshold::InnerSolver(const int aStart, const int aEnd, const int aI, const int aJ, const int aK,
	const boostMap& control_list, const array_1D& speed_list, CGrid_2D& aGrid, CENO& aENOgrid,
	CSpecial_Treatment& aSpecial_treatment, const string aInterp_Choice, const bool aReturn_Wswitch,
	array<double, 2>& aW_update, array<double, 2>& aCtrl_update, array<double, 2>& aW_switch) {

	int aIndx_current = -10;
	array<bool, 2> Switch_flag_2modes{ false,false };

	if (aI * fDs > fC) {
		aIndx_current = fNum_slices - 1;
	}
	else {
		aIndx_current = aI;
	}

	for (int cntrli = 0; cntrli < fNumCtrls; cntrli++)
	{
		double  Current_Ctrl = control_list.left.at(cntrli + 1);
		double Current_speed = speed_list[cntrli];
		double aW_possible = -1;

		for (int Current_mode = aStart; Current_mode < aEnd; Current_mode++) {
			if (aInterp_Choice == "ENO") {
				aW_possible = Marching_Step(aIndx_current, aJ, aK, Current_mode + 1, Current_Ctrl, Current_speed, aENOgrid, fInterp_Order, fMyValfn);
			}
			else if (aInterp_Choice == "Bilinear") {
				aW_possible = Marching_Step(aIndx_current, aJ, aK, Current_mode + 1, Current_Ctrl, Current_speed, aGrid, fMyValfn);
			}

			if (aW_possible > aW_update[Current_mode])
			{
				aW_update[Current_mode] = min(aW_possible, 1.0);
				aCtrl_update[Current_mode] = Current_Ctrl;
			}
		}
	}
	//Call the Golden section search fucntion if we decide to use it
	if (fIf_use_GoldenSS) {

		Golden_section_search_in_Marching(aStart, aEnd, aIndx_current, aJ, aK, control_list, speed_list, aGrid, aENOgrid, aInterp_Choice, aW_update, aCtrl_update);
	}


	//if there is enough budget to switch
	if (aI * fDs > fC)
	{
		for (int Current_mode = aStart; Current_mode < aEnd; Current_mode++) {

			if (aInterp_Choice == "ENO") {
				aW_switch[Current_mode] = Compute_Switch_Value(aI, aJ, aK, aENOgrid, fInterp_Order, Current_mode + 1, fMyValfn);
			}
			else if (aInterp_Choice == "Bilinear") {
				aW_switch[Current_mode] = Compute_Switch_Value(aI, aJ, aK, aGrid, Current_mode + 1, fMyValfn);
			}

			Switch_Policy_Update(aIndx_current, aJ, aK,
				aW_switch[Current_mode], aW_update[Current_mode], aCtrl_update[Current_mode], Switch_flag_2modes[Current_mode],
				Current_mode + 1, fMySwitch, fMyPolicy, fMyPolicy_out);
		}
	}
	else {
		for (int Current_mode = aStart; Current_mode < aEnd; Current_mode++) {
			Switch_Policy_Update(aIndx_current, aJ, aK,
				aW_switch[Current_mode], aW_update[Current_mode], aCtrl_update[Current_mode], Switch_flag_2modes[Current_mode],
				Current_mode + 1, fMySwitch, fMyPolicy, fMyPolicy_out);
		}
	}

	for (int Current_mode = aStart; Current_mode < aEnd; Current_mode++) {
		//Write the maximized probability value into the 4D value function array
		fMyValfn[Current_mode][aIndx_current][aJ][aK] = (aW_update[Current_mode] > fTol_zero) ? aW_update[Current_mode] : 0.0;

		//Convert it to short type for output arrays
		fMyValfn_out[Current_mode][aIndx_current][aJ][aK] = (short)floor(fMyValfn[Current_mode][aIndx_current][aJ][aK] * 1e4);

		//Special treatments for probabilities near 1 (prepared for "Deadline-upgrade")
		aSpecial_treatment.Assign_if_Near_1(aIndx_current, aI, aJ, aK, Current_mode,
			aW_update[Current_mode], aCtrl_update[Current_mode], Switch_flag_2modes[Current_mode], fMyCheck, fMyValfn, fMyPolicy, fMySwitch, fMyValfn_out, fMyPolicy_out);

		//Special treatments for probabilities near 0 (need to read in Risk-Neutral policy)
		aSpecial_treatment.Assign_if_Near_0(aIndx_current, aJ, aK, Current_mode,
			aW_update[Current_mode], fMyValfn, fMyPolicy, fMySwitch, fRisk_Neutral_Policy, fRisk_Neutral_Switch, control_list, fMyValfn_out, fMyPolicy_out);

		//Record W_switch if we need it
		if (aReturn_Wswitch)
		{
			fMyValfn_force_to_switch[Current_mode][aIndx_current][aJ][aK] = (aW_switch[Current_mode] > fTol_zero) ? aW_switch[Current_mode] : 0.0;

			//Convert it to the "short" type for output arrays
			fMyValfn_force_to_switch_out[Current_mode][aIndx_current][aJ][aK] = (short)floor(fMyValfn_force_to_switch[Current_mode][aIndx_current][aJ][aK] * 1e4);

			aSpecial_treatment.Assign_if_Near_1_force_to_switch(aIndx_current, aI, aJ, aK, Current_mode,
				aW_switch[Current_mode], fMyCheck_force_to_switch, fMyValfn_force_to_switch, fMyValfn_force_to_switch_out);

			aSpecial_treatment.Assign_if_Near_0_force_to_switch(aIndx_current, aJ, aK, Current_mode, aW_switch[Current_mode],
				fMyValfn_force_to_switch, fMyValfn_force_to_switch_out);
		}
	}
}


//============================Methods for if we are forced to switch================================
void Sailing_Threshold::InnerSolver_force_to_switch(const int aStart, const int aEnd, const int aI, const int aJ, const int aK, const array_3D& aCheck_mat,
	CSpecial_Treatment aSpecialTreatmentClass, CGrid_2D& aGrid, CENO& aENOgrid)
{

	for (int Current_mode = aStart; Current_mode < aEnd; Current_mode++) {
		if ((aI > aCheck_mat[0][aJ][aK]) || (fMySwitch[Current_mode][fNum_slices - 1][aJ][aK] == 1))
		{
			fMyValfn_force_to_switch[Current_mode][fNum_slices - 1][aJ][aK] = 1.0;

			//Convert it to short type for output arrays
			fMyValfn_force_to_switch_out[Current_mode][fNum_slices - 1][aJ][aK] = (short)(1.0 * 1e4);
		}
		else {

			double aW_force_to_switch = Compute_Switch_Value(aI, aJ, aK, aENOgrid, fInterp_Order, Current_mode + 1, fMyValfn);

			fMyValfn_force_to_switch[Current_mode][fNum_slices - 1][aJ][aK] = (aW_force_to_switch > fTol_zero) ? aW_force_to_switch : 0.0;

			fMyValfn_force_to_switch_out[Current_mode][fNum_slices - 1][aJ][aK] = (short)round(fMyValfn_force_to_switch[Current_mode][fNum_slices - 1][aJ][aK] * 1e4);

			aSpecialTreatmentClass.Assign_if_Near_1_force_to_switch(fNum_slices - 1, aI, aJ, aK, Current_mode,
				aW_force_to_switch, fMyCheck_force_to_switch, fMyValfn_force_to_switch, fMyValfn_force_to_switch_out);

			aSpecialTreatmentClass.Assign_if_Near_0_force_to_switch(fNum_slices - 1, aJ, aK, Current_mode,
				aW_force_to_switch, fMyValfn_force_to_switch, fMyValfn_force_to_switch_out);
		}
	}
}


void Sailing_Threshold::MainSolver() {

	// Read in speed profile data
	load_speed_data(fSpeed_File);

	// Create grid object
	CGrid_2D myGrid(fDr, fR0, fR, fDtheta, fTheta0, fTheta, fDs);
	CENO myEnoGrid(fDr, fR0, fR, fDtheta, fTheta0, fTheta, fDs);
	CSpecial_Treatment mySpecial_Treatment(fNum_slices, fR, fTheta, fTol_1, fTol_zero);
	boostMap control_list; //list of possible controls
	array_1D speed_list(boost::extents[fNumCtrls]); //speed_list[i] is speed corresponding to control_list[i]

	// Initialize 1d arrays that store corresponding speed for each control
	Initialize_Control_and_Speed_List(control_list, speed_list, myGrid);

	// Initialization of matrices which store value function values, optimal policies, and switchgrids
	InitializeMat();

	// Read in Risk Neutral policy
	Read_Risk_Neutral_Policy();

	// Create date files and read the initial ones into them
	io::writeToFile4D_init_out<short>("test_ValueFunc.dat", fMyValfn_out, 0);
	io::writeToFile4D_init_out<short>("test_Policies.dat", fMyPolicy_out, 0);
	io::writeToFile4D_init_out<short>("test_Switchgrid.dat", fMySwitch, 0);
	if (fReturn_Wswitch) {
		io::writeToFile4D_init_out<short>("test_ForceSwitch.dat", fMyValfn_force_to_switch_out, 0);
	}

	int target_indx = int(fR_target / fDr);

	//budget loop
	for (int i = 1; i < fM + 1; i++)
	{
		auto start = std::chrono::steady_clock::now();
		bool ReachedTopSlice = false;
		if (i >= fNum_slices) {
			UpdateSlices();
			ReachedTopSlice = true;
		}


		#pragma omp parallel for schedule(static,1)
		for (int j = 0; j < fR + 1; j++)
		{//radius loop
			double Current_R = fR0 + j * fDr;
			bool Impossible_to_reach = false;
			//if we are not in target
			if (Current_R > fR_target)
			{
				// if the distance divided by the maximum speed (default to 2 above) is larger than the current budget,
				// it is impossible to reach the target within the current budget
				if ((Current_R - fR_target) / fMax_speed > i * fDs) {
					Impossible_to_reach = true;
				}
				//theta loop
				for (int k = 0; k < fTheta; k++)
				{
					if (Impossible_to_reach) {
						for (int mode = 0; mode < 2; mode++) {
							//if impossible to reach, the probability of success in either mode is 0
							if (ReachedTopSlice) {
								mySpecial_Treatment.Assign_if_Exactly_0(mode, fNum_slices - 1, j, k, fRisk_Neutral_Policy, fRisk_Neutral_Switch,
									fMyValfn, fMySwitch, fMyPolicy, control_list, fMyValfn_out, fMyPolicy_out);
							}
							else {
								mySpecial_Treatment.Assign_if_Exactly_0(mode, i, j, k, fRisk_Neutral_Policy, fRisk_Neutral_Switch,
									fMyValfn, fMySwitch, fMyPolicy, control_list, fMyValfn_out, fMyPolicy_out);
							}
						}
					}
					else {
						//Initilization for the inner solver (localize it be make parallelization easier)
						array<double, 2> W_update{ -1,-1 };
						array<double, 2> Ctrl_update{ -1,-1 };
						array<double, 2> W_switch{ -1,-1 };

						int Start_here = -10;
						int End_here = -10;
						int Start_here_ddlupg = -10;
						int End_here_ddlupg = -10;

						//InnerSolver_pointer(ReachedTopSlice, Start_here, End_here, Start_here_ddlupg, End_here_ddlupg, i, j, k, fMyCheck);
						std::tie(Start_here, End_here, Start_here_ddlupg, End_here_ddlupg) = InnerSolver_pointer(i, j, k, fMyCheck);
						InnerSolver(Start_here, End_here, i, j, k, control_list, speed_list, myGrid, myEnoGrid, mySpecial_Treatment,
							fInterp_Choice, fReturn_Wswitch, W_update, Ctrl_update, W_switch);

						mySpecial_Treatment.Deadline_upgrade_policy(ReachedTopSlice, Start_here_ddlupg, End_here_ddlupg, i, j, k,
							fMyValfn, fMyPolicy, fMySwitch, fMyValfn_out, fMyPolicy_out);

						if (fReturn_Wswitch)
						{
							InnerSolver_force_to_switch(Start_here_ddlupg, End_here_ddlupg, i, j, k, fMyCheck_force_to_switch, mySpecial_Treatment, myGrid, myEnoGrid);
						}
					}
				}
			}
		}

		if (i % fRefinement_factor == 0)
		{
			//Reshaffle slices only if the current budget is beyond the switching cost
			if (i >= fNum_slices) {
				//Dump the newest slice to the disk
				io::AppendToFile4D_out<short>("test_ValueFunc.dat", fMyValfn_out, fNum_slices - 1);
				io::AppendToFile4D_out<short>("test_Policies.dat", fMyPolicy_out, fNum_slices - 1);
				io::AppendToFile4D_out<short>("test_Switchgrid.dat", fMySwitch, fNum_slices - 1);
				if (fReturn_Wswitch) {
					io::AppendToFile4D_out<short>("test_ForceSwitch.dat", fMyValfn_force_to_switch_out, fNum_slices - 1);
				}
			}
			else {
				io::AppendToFile4D_out<short>("test_ValueFunc.dat", fMyValfn_out, i);
				io::AppendToFile4D_out<short>("test_Policies.dat", fMyPolicy_out, i);
				io::AppendToFile4D_out<short>("test_Switchgrid.dat", fMySwitch, i);
				if (fReturn_Wswitch) {
					io::AppendToFile4D_out<short>("test_ForceSwitch.dat", fMyValfn_force_to_switch_out, i);
				}
			}
		}

		cout << "The current slice is i = " << i << endl;
		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

	}

}


void Sailing_Threshold::writeDomainToFile(const string aFilename) const {

	const double Theta_dim = getThetaGridSize(); // index of last theta point
	const double theta_max = getThetaMax(); // largest point at theta-axis
	const double fR = getRadiusGridSize(); // index of last r point
	const double r_max = getRadiusMax(); // largest point at r-axis
	const double g_factor = getGFactor(); // for grid refinement studies
	const double target_radius = getTargetSize(); // radius of target

	array_2D params(boost::extents[1][18]); // Creating boostarray "row vector"

	params[0][0] = Theta_dim;
	cout << "fTheta = " << params[0][0] << endl;

	params[0][1] = theta_max;
	cout << "theta_max = " << params[0][1] << endl;

	params[0][2] = fR;
	cout << "fR = " << params[0][2] << endl;

	params[0][3] = r_max;
	cout << "r_max = " << params[0][3] << endl;

	params[0][4] = g_factor;
	cout << "g_factor = " << params[0][4] << endl;

	params[0][5] = target_radius;
	cout << "Target Radius = " << params[0][5] << endl;

	params[0][6] = fM;
	cout << "fM = " << params[0][6] << endl;

	params[0][7] = fBudget;
	cout << "Budget = " << params[0][7] << endl;

	params[0][8] = fR0;
	cout << "r_min = " << params[0][8] << endl;

	params[0][9] = fTheta0;
	cout << "theta_min = " << params[0][9] << endl;

	params[0][10] = fDr;
	cout << "Delta r = " << params[0][10] << endl;

	params[0][11] = fDtheta;
	cout << "Delta theta = " << params[0][11] << endl;

	params[0][12] = fDs;
	cout << "Delta budget = " << params[0][12] << endl;

	params[0][13] = fSigma;
	cout << "Wind diffusion = " << params[0][13] << endl;

	params[0][14] = fTau;
	cout << "Tau = " << params[0][14] << endl;

	params[0][15] = fC;
	cout << "Cost to Switch = " << params[0][15] << endl;

	params[0][16] = fDrift;
	cout << "Wind Drift = " << params[0][16] << endl;

	params[0][17] = fNumCtrls;
	cout << "Number of Controls = " << params[0][17] << endl;

	// Writing vector of parameters to file
	io::writeToFile2D<double>(aFilename + "DomainParameters", params);
}
