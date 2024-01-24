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
 * File: main.cpp
 *
 * Author: MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu
 *
 * Description: This file initializes all the global variables and executes the
 * the corresponding example from the command line.
 *
 *============================================================================*/

//-----------------------Project specific header files---------------------------
#include "Sailing_Threshold_Solver.h"

//----------------------------Libraries------------------------------------------
#include <chrono>
using namespace std;

int main()
{
	//=============================================Initilization=====================================================
	const int gRefinement_factor = 2; //Global refinement factor with respect to the spatial grid
	const int gslice = 5; //Additional global refinement factor with respect to the threshold dimension
	const double gBudget = 55; //Global max budget we are integrating up to
	const double gDiffConst = 0.05; //Global Diffusion constant ("sigma")
	const double gDrift = 0.05;  //Global drift constant ("a")
	const double gCost_to_switch = 2.0; //Global switching cost
	const double gR_target = 0.1; //Global target radius
	const int gNumCtrls = 49; //Global number of core (base) controls
	const string gInterpChoice = "ENO"; //Two choices: (1) "Bilinear" and (2) "ENO"
	const int gInterp_Order_ENO = 3;   //Two choices:  (1) 2 (ENO bi-quadratic) and (2) 3 (ENO bi-cubic)
	const int gNum_GH_pts = 3; //Global number of points of the Gauss Hermite Quadrature
	const double gTol_prob_near_1 = 1e-13; //Global tolerance for prob near 1
	const double gTol_prob_near_0 = 1e-13; //Global tolerance for prob near 0
	string gSpeed_File = "Sunodyssey40.pol"; //Global speed profile data file (Fixed for all our numerical experiments)
	const int gWind_speed_index = 3; //Global wind speed index (Fixed for all our numerical experiments)
	const double gMax_speed = 0.05; //Global max speed of the sailboat (Fixed for all our numerical experiments)
	string gRN_Policy_path = "../Risk_Neutral/output/PathTracer_policy_"; //filepath for the corresponding risk-neutral optimal steering angles
	string gRN_Switch_path = "../Risk_Neutral/output/PathTracer_switch_"; //filepath for the corresponding risk-neutral switchgrid
	speed_method gSpeed_choice = speed_method::INTERPOLATE; //Global method of interpreting/discretizing the speed profile. Three choices: (1) SETCONTROL, (2) INTERPOLATE, and (3) CEIL
	bool gReturn_Wswitch = false; //Decide whether to return the Value function if we force the boat to switch "initially"
	bool gIf_use_GoldenSS = true; //Decide whether to use Golded section search (GSS) or not
	double gTol_GoldenSS = 1e-3; //Tolerance for GSS

	//=============================================Start of Our PDE Solver=====================================================
	auto start = chrono::high_resolution_clock::now();
	cout << "Running  Risk-Aware Sailing Solver for our 2D example" << endl;

	//allocation heap/stack
	//Build the main class of our PDE solver
	Sailing_Threshold Example2D(gRefinement_factor, gslice, gBudget, gDiffConst, gCost_to_switch, gR_target, gDrift, gNumCtrls,
		gSpeed_File, gSpeed_choice, gRN_Policy_path, gRN_Switch_path, gWind_speed_index, gMax_speed, gReturn_Wswitch,
		gInterpChoice, gInterp_Order_ENO, gTol_prob_near_1, gTol_prob_near_0, gNum_GH_pts, gIf_use_GoldenSS, gTol_GoldenSS);

	//Calling the main PDE solver
	Example2D.MainSolver();

	// Writing Grid Parameters to file
	std::string file_label = std::to_string(gRefinement_factor);
	Example2D.writeDomainToFile("PathTracer_");

	auto end = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Total runtime: " << elapsed_seconds.count() << "s\n";
	cout << "Successfully completed!" << endl;
	//=============================================END of Our PDE Solver=====================================================
	return 0;
}
