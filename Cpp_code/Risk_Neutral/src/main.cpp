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
#include "CGrid_2D.h"
#include "Sailing_Solver.h"

//----------------------------Libraries------------------------------------------
#include <chrono>

using namespace std;

int main()
{
	const int gFactor = 2; //Global refinement factor with respect to the spatial grid
	const double  gDiffConst = 0.05;  //Global Diffusion constant ("sigma")
	const double  gDrift = 0.05;  //Global drift constant ("a")
	const double  gCost_to_switch = 2.0;  //Global switching cost
	const double  gR_target = 0.1; //Global target radius
	const int gNumCtrls = 49; //Global number of core (base) controls
	const double  gTol = 1e-7; //Global tolerance to check the convergence of Value Iterations
	speed_method gSpeed_choice = speed_method::INTERPOLATE; //Global method of interpreting/discretizing the speed profile. Three choices: (1) SETCONTROL, (2) INTERPOLATE, and (3) CEIL
	string gSpeed_File = "Sunodyssey40.pol"; //Global speed profile data file (Fixed for all our numerical experiments)


	//=============================================Start of Our PDE Solver=====================================================
	auto start = chrono::high_resolution_clock::now();
	cout << "Running Risk-Neutral Sailing Solver for our 2D example" << endl;
	//allocation heap/stack
	Sailing_Solver Example2D(gFactor, gDiffConst, gR_target, gDrift, gNumCtrls, gTol, gSpeed_File, gSpeed_choice, gCost_to_switch);
	//Calling the main PDE solver
	Example2D.solve_value_sweep(true);

	// Writing Grid Parameters to file
	Example2D.writeDomainToFile("PathTracer_");

	//CGrid_1D mygrid(dx, x0, N);
	auto end = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Total runtime: " << elapsed_seconds.count() << "s\n";
	cout << "Successfully completed!" << endl;
	//=============================================END of Our PDE Solver=====================================================
	return 0;
}
