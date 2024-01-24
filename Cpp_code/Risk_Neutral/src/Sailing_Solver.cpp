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
 * File: Sailing_Solver.cpp
 *
 * Author: MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu
 *
 * Description: This file contains the declarations of functions that compute the
 * drift / diffusive function of the dynamics;
 * and carry out (Row-wise) Value Iterations in the main solver based on a
 * first-order semi-Lagrangian discretization.
 *
 * Details of all of these functions are found in Sailing_Solver.cpp.
 *
 *============================================================================*/

 //----------------------Project specific header files---------------------------
#include "Sailing_Solver.h"
#include "WriteToFile.h"

//------------------------------Libraries------------------------------------
#include<chrono>
#include<omp.h>
#include<numeric>


//----------------------------Definitions----------------------------------------
using namespace std;

//Define the constructor
Sailing_Solver::Sailing_Solver(int aFactor, double aDiffConst, double aTarget_radius, double aDrift,
	int aNumCtrls, double aTol_converge, string aSpeed_File, speed_method aSpeed_choice, double aCost_to_switch) {
	fPi = 3.14159265358979323846 ; // constant for pi
	fMulti_factor = aFactor; // An refinment factor with respect to the spatial grid (Base case: 100*2)
	fR0 = 0; // starting point at r-axis
	fTheta0 = 0; // starting point at theta-axis
	fRmax = 2.0; // largest point at r-axis
	fThetamax = 2 * fPi; // largest point at theta-axis
	fR = 200 * fMulti_factor; // index of the last r point
	fTheta = 200 * fMulti_factor; // index of the last theta point
	fDr = (fRmax - fR0)/(double)fR; // Delta r
	fDtheta = (fThetamax - fTheta0)/(double)fTheta; // Delta theta
	fSigma = aDiffConst;  // Diffusion constant
	fR_target = aTarget_radius; //radius of target set
	fDrift = aDrift; //drift ("a")

	fC = aCost_to_switch; //Cost to switch
	fTol = aTol_converge; //Convergence error tolerance
	fMax_iters = 500; //Max number of iterations
	fLambda = 0.0; // discounting factor (always set to 0 for all of our numerical experiments)
	fSpeed_File = aSpeed_File;  //filename of the speed profile
	fMax_speed = 0.05;  //max speed of the sailboat
	fSpeed_Method = aSpeed_choice; //3 choices: (1) SETCONTROL: sets controls as angles given in .pol file,
								   //           (2) INTERPOLATE: takes given control angles and interpolates between angles given in .pol file to find speed
								   //           (3) CEIL: takes given control angles and uses next highest speed
	if (fSpeed_Method == speed_method::SETCONTROL) {
		fNumCtrls = 24; //number of core controls (fixed for all of our numerical experiments according to the data file)
	}
	else {
		fNumCtrls = aNumCtrls; //number of controls
	}
	fWind_Speed_Index = 3; //wind speed index in the .pol file
 	fMyValfn.resize(boost::extents[2][fR + 1][fTheta + 1]);
	fMyPolicy.resize(boost::extents[2][fR + 1][fTheta + 1]);
	fMySwitch.resize(boost::extents[2][fR + 1][fTheta + 1]);


	//Roots of 3rd Hermite Polynomial
	fHermite_root = { -1.0 * sqrt(3.0 / 2.0), 0.0, sqrt(3.0 / 2.0)};
	//Gauss-Hermite Quadrature weights (absorbed 1/sqrt(pi) into it)
	fHermite_weight = {1.0/6.0, 2.0/3.0, 1.0/6.0};
}

// Initialization of matrices which store value function values, optimal policies, and switchgrids
void Sailing_Solver::InitializeMat() {

	int  target_idx1 = int(floor((fR_target - fR0) / fDr));

	for (int q = 0; q < 2; q++) {
		for (int ir = 0; ir < fR + 1; ir++) {
			for (int itheta = 0; itheta < fTheta + 1; itheta++) {
				if (itheta == fTheta) {//checking for being nonsense place (theta = 2pi)
					fMyValfn[q][ir][itheta] = -1000;
					fMyPolicy[q][ir][itheta] = -1000;
					fMySwitch[q][ir][itheta] = -1000;
				}
				else
				{
					if (ir <= target_idx1)
					{  // at 0 or 1 in state space ...
						fMyValfn[q][ir][itheta] = 0;
						fMyPolicy[q][ir][itheta] = 0;
						fMySwitch[q][ir][itheta] = 0;
					}
					else
					{
						fMyValfn[q][ir][itheta] = 1e6;
						fMyPolicy[q][ir][itheta] = 0;
						fMySwitch[q][ir][itheta] = 0;
					}
				}
			}
		}
	}
}


//This function implements the binary search algorithm
int Sailing_Solver::binary_search(const vector<double > aVector, const double  aXloc)
{
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
void Sailing_Solver::load_speed_data(string filename) {
	ifstream infile(filename);
	string line;

	int i = 0;
	while (getline(infile, line)) // Read a line
	{
		if (i != 0) {
			fSpeed_Data.push_back(vector<double >()); //Speed data holds data from .pol file
			vector<double >& row = fSpeed_Data.back();
			istringstream iss(line);
			double  value;
			int j = 0;
			while (iss >> value) {
				if (j != 0) {
					row.push_back(value);
				}
				if (j == 0) {
					fAngles.push_back(fPi * value / 180); //angles holds list of angles given in .pol file
				}
				j++;
			}

		}
		i++;

	}

	infile.close();
}


// These functions compute the time derivatives in the differential equation, missing a factor of the boat speed
// r - dist, θ - angle, q - tack, u - steering angle
inline double  Sailing_Solver::dr(const double  aTheta, const double  aQ, const double  aU)
{
	return -1.0 * cos(aTheta - pow(-1.0, aQ) * aU);
}

inline double  Sailing_Solver::dtheta(const double  aR, const double  aTheta, const double  aQ, const double  aU)
{
	return (1.0 / aR) * sin(aTheta - pow(-1.0, aQ) * aU);
}

//Determine the size of the time step
inline double  Sailing_Solver::dt_disc(const double  aSr, const double  aStheta, const bool aRowwise)
{
	if (aRowwise == true) {
		return 1.5 * fDr / (abs(aSr) + 1e-8);
	}
	else {
		return min(1.5 * fDr / (abs(aSr) + 1e-8), 1.5 * fDtheta / (abs(aStheta) + 1e-8));
	}

}



//Defines list of controls to iterate over and their respective speeds
void Sailing_Solver::Initialize_Control_and_Speed_List(array_1D& aControl_List, array_1D& aSpeed_List, CGrid_2D& aGrid)
{
	if (fSpeed_Method == speed_method::SETCONTROL)
	{
		aControl_List.resize(boost::extents[fAngles.size()]);
		aSpeed_List.resize(boost::extents[fAngles.size()]);
		for (int n = 0; n < aControl_List.size(); n++)
		{
			aControl_List[n] = fAngles.at(n);
			aSpeed_List[n] = fSpeed_Data[n][fWind_Speed_Index];
		}
	}
	else if (fSpeed_Method == speed_method::INTERPOLATE) //Interpolate between angles in the speed plot
	{
		aControl_List.resize(boost::extents[fNumCtrls]);
		aSpeed_List.resize(boost::extents[fNumCtrls]);
		array_1D speed_column(boost::extents[fAngles.size()]);
		for (int n = 0; n < fAngles.size(); n++)
		{
			speed_column[n] = fSpeed_Data[n][fWind_Speed_Index];
		}
		for (int n = 0; n < fNumCtrls; n++)
		{
			double  ctrl = n * fPi / (fNumCtrls - 1);
			aControl_List[n] = ctrl;
			int index = binary_search(fAngles, ctrl);

			array<double, 2> speed_stencil = aGrid.Stencil_for_Linear_Interp(speed_column, index);
			array<double, 2> angle_stencil;
			assert(index > 0);
			angle_stencil[0] = fAngles.at(index - 1);
			angle_stencil[1] = fAngles.at(index);
			aSpeed_List[n] = aGrid.Linear_Interp(speed_stencil, angle_stencil, ctrl);
		}
	}
	else if (fSpeed_Method == speed_method::CEIL) //take the speed corresponding to first angle greater than control
	{
		aControl_List.resize(boost::extents[fNumCtrls]);
		aSpeed_List.resize(boost::extents[fNumCtrls]);
		for (int n = 0; n < fNumCtrls; n++)
		{
			double  ctrl = n * fPi / (fNumCtrls - 1);
			aControl_List[n] = ctrl;
			int index = binary_search(fAngles, ctrl);
			aSpeed_List[n] = fSpeed_Data[index][fWind_Speed_Index];
		}
	}

	double  max_speed = 0;
	for (int n = 0; n < aSpeed_List.size(); n++) {
		if (aSpeed_List[n] > max_speed) {
			max_speed = aSpeed_List[n];
		}
	}

	for (int n = 0; n < fNumCtrls; n++) {
		aSpeed_List[n] = aSpeed_List[n] * fMax_speed / max_speed; // CHANGED HERE
	}

	// Writing control values and corresponding speeds to file
	io::writeToFile1D<double >("ControlList.dat", aControl_List);
	io::writeToFile1D<double >("SpeedList.dat", aSpeed_List);

}

// Main function that computes the expected value if we switched mode using Gauss Hermite Quadrature
double  Sailing_Solver::gauss_hermite_switch(int aIr, int aItheta, int aQ, CGrid_2D& aGrid) {
	int  qprime = 3 - aQ;
	//current theta location
	double  theta = fTheta0 + aItheta * fDtheta;

	double  value_estimate = 0.0;
	for (int m = 0; m < 3; m++) {
		double  theta_prime = sqrt(2.0 * fC) * fSigma * fHermite_root[m] + theta + fDrift * fC;
		theta_prime = mod_2pi(theta_prime);
		int  theta_new_ind = find_index(theta_prime, fDtheta, fTheta0, fThetamax);
		array<double, 2> v3_stencil = aGrid.Stencil_for_Linear_Interp_Theta(fMyValfn, qprime, aIr, theta_new_ind);
		double  switch_value_possible = aGrid.Linear_Interp_Theta(v3_stencil, theta_new_ind, theta_prime);
		value_estimate += fHermite_weight[m] * switch_value_possible;
	}
	return value_estimate;
}

//This is the main solver that uses (Row-wise) Value Iterations to compute the value function and the corresponding optimal policy
void Sailing_Solver::solve_value_sweep(const bool aRowwise){

	bool rowwise = aRowwise; //Decide whether to use "Row-wise" algorithm or not

	// Read in speed profile data
	load_speed_data(fSpeed_File);

	// Create grid object
	CGrid_2D myGrid(fDr, fR0, fR, fDtheta, fTheta0, fTheta);

	array_1D control_list; //list of possible controls
	array_1D speed_list; //speed_list[i] is speed corresponding to control_list[i]

	// Initialize 1d arrays that store corresponding speed for each control
	Initialize_Control_and_Speed_List(control_list, speed_list, myGrid);

	// Initialization of matrices which store value function values, optimal policies, and switchgrids
	InitializeMat();

	//At less than this index, we have hit the target
	int  target_idx = int(floor(fR_target/fDr));

	//Value iterations!
	double max_diff = INFINITY;
	int it = 0;
	array<double , 2> step_lengths {double (1.5 * fDr), double (1.5 * fDtheta)};
	while (!(max_diff < fTol)){
		max_diff = 0.0;
		for(int ir = target_idx+1; ir < fR + 1; ir++){
			for(int q = 1; q < 3; q++){
				for(int itheta = 0; itheta < fTheta; itheta++){
					double  r = fR0 + ir * fDr;
					double  theta = fTheta0 + itheta * fDtheta;
					double  stay_value_possible = INFINITY;
					double  stay_value_min = INFINITY;
					double  bestu = 0;
					for(int i = 0; i < control_list.size(); i++){
						//compute the time derivatives
						double  speed_r = speed_list[i] * dr(theta, q, control_list[i]);
						double  speed_theta = speed_list[i] * dtheta(r, theta, q, control_list[i]) + fDrift;
						double  dt = dt_disc(speed_r, speed_theta, rowwise); //decide the time step

						//foot of characteristics
						double  newr = min(r + speed_r * dt,fRmax);
						double  newtheta = theta + speed_theta * dt;
						//foot of sample points (+- with equal probability representing a weak approximation of the 1D Brownian motion)
						double  newtheta_plus = mod_2pi(newtheta + sqrt(dt) * fSigma);
						double  newtheta_minus = mod_2pi(newtheta - sqrt(dt) * fSigma);

						//compute the right indices in the theta-direction
						int newtheta_plus_ind = find_index(newtheta_plus, fDtheta, fTheta0, fThetamax);
						int newtheta_minus_ind = find_index(newtheta_minus, fDtheta, fTheta0, fThetamax);
						//compute the right index in the r-direction
						int newr_ind = find_index(newr, fDr, fR0, fRmax);

						//Construct stencils and approximate the value function by bi-linear interpolation
						array<double , 4> v0_stencil = myGrid.Stencil_for_Bilinear_Interp(fMyValfn, q, newr_ind, newtheta_plus_ind);
						array<double , 4> v1_stencil = myGrid.Stencil_for_Bilinear_Interp(fMyValfn, q, newr_ind, newtheta_minus_ind);
						stay_value_possible = dt + exp(-1.0 *fLambda * dt) * 0.5 * (
							 myGrid.Bilinear_Interp(v0_stencil, newr_ind, newtheta_plus_ind, newr, newtheta_plus)
							 + myGrid.Bilinear_Interp(v1_stencil, newr_ind, newtheta_minus_ind, newr, newtheta_minus)
							 );

						 if(stay_value_possible < stay_value_min){
							 stay_value_min = stay_value_possible; //update the value function if smaller
							 bestu = control_list[i];
						 }
					}
					double  V_delta = gauss_hermite_switch(ir, itheta, q, myGrid) + fC; //compute the value function if immediately switching tack

 				 	fMyPolicy[q-1][ir][itheta] = bestu; //record the "best" steering angle

 				 	if(V_delta < stay_value_min){
						//if the value function for switching is smaller, we switch tacks
					 	max_diff = max(max_diff, double(abs(fMyValfn[q-1][ir][itheta] - V_delta)));
					 	fMySwitch[q-1][ir][itheta] = 1; // "1" means switching mode
					 	fMyValfn[q-1][ir][itheta] = V_delta;
						fMyPolicy[q - 1][ir][itheta] = -1000; //policy is -1000 where optimal to switch
					 }
 				  else{
						//Otherwise, we stay on the current tack
					 	if(!isinf(stay_value_min)){
							max_diff = max(max_diff, double(abs(fMyValfn[q-1][ir][itheta] - stay_value_min)));
						}
					 	fMySwitch[q-1][ir][itheta] = 0;  // "0" means staying on the current tack
					 	fMyValfn[q-1][ir][itheta] = stay_value_min;
						fMyPolicy[q - 1][ir][itheta] = bestu;
				 	}
				}
			}
		}
		cout << "Max_diff: " << max_diff << endl;
		if(fmod(it, 10) == 0){
			cout << "Iteration " << it << ", Maximum Difference " << max_diff << endl;
		}
		if(it == fMax_iters){
			cout << "Stopped after " << fMax_iters << " iterations";
			return;
		}
		it = it + 1;
	}
	cout << "Finished successfully in " << it << " iterations, with final max_diff " << max_diff << endl;

	writeMatstoFile();
	return;
}



void Sailing_Solver::writeMatstoFile()
{
    // Writing optimal policy matrix for Mode 1 to File
	array_2D myPolicy_1(boost::extents[fR + 1][fTheta]);
	for(int j = 0; j < fR + 1; j++){
		for(int k = 0; k < fTheta; k++){
			myPolicy_1[j][k] = fMyPolicy[0][j][k];
		}
	}
	string filename3 = "PathTracer_policy_1.dat";
	io::writeToFile2D<double >(filename3, myPolicy_1);

    // Writing optimal policy matrix for Mode 2 to File
	array_2D myPolicy_2(boost::extents[fR + 1][fTheta]);
	for(int j = 0; j < fR + 1; j++){
		for(int k = 0; k < fTheta; k++){
			myPolicy_2[j][k] = fMyPolicy[1][j][k];
		}
	}
	string filename4 = "PathTracer_policy_2.dat";
	io::writeToFile2D<double >(filename4, myPolicy_2);

    // Writing switchgrid matrix for Mode 1 to File
	array_2D mySwitch_1(boost::extents[fR + 1][fTheta]);
	for(int j = 0; j < fR + 1; j++){
		for(int k = 0; k < fTheta; k++){
			mySwitch_1[j][k] = fMySwitch[0][j][k];
		}
	}
	string filename5 = "PathTracer_switch_1.dat";
	io::writeToFile2D<double >(filename5, mySwitch_1);

    // Writing switchgrid matrix for Mode 2 to File
	 array_2D mySwitch_2(boost::extents[fR + 1][fTheta]);
	for(int j = 0; j < fR + 1; j++){
		for(int k = 0; k < fTheta; k++){
			mySwitch_2[j][k] = fMySwitch[1][j][k];
		}
	}
	string filename6 = "PathTracer_switch_2.dat";
	io::writeToFile2D<double >(filename6, mySwitch_2);

	// Writing value function matrix for Grid Refinement Mode 1 to File
	// string filename7 = "Control_Refinement/" + to_string(fNumCtrls) + "Controls_valuefn.dat";
	// io::writeToFile3D<double >(filename7, fMyValfn);


	io::writeToFile3D<double >("PathTracer_ValueFunc.dat", fMyValfn);
 	io::writeToFile3D<double >("PathTracer_Policies.dat", fMyPolicy);
 	io::writeToFile3D<double >("PathTracer_Switchgrid.dat", fMySwitch);

}


void Sailing_Solver::writeDomainToFile(const string aFilename) const {

		const double  fTheta = getThetaGridSize(); // index of last theta point
		const double  theta_max = getThetaMax(); // largest point at theta-axis
		const double  fR = getRadiusGridSize(); // index of last r point
		const double  r_max = getRadiusMax(); // largest point at r-axis
		const double  g_factor = getGFactor(); // for grid refinement studies
		const double  target_radius = getTargetSize(); // radius of target

		array_2D params(boost::extents[1][14]); // Creating boostarray "row vector"

  	params[0][0] = fTheta;
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

		params[0][6] = fR0;
		cout << "r_min = " << params[0][6] << endl;

		params[0][7] = fTheta0;
		cout << "theta_min = " << params[0][7] << endl;

		params[0][8] = fDr;
		cout << "Delta r = " << params[0][8] << endl;

		params[0][9] = fDtheta;
		cout << "Delta theta = " << params[0][9] << endl;

		params[0][10] = fSigma;
		cout << "Wind diffusion = " << params[0][10] << endl;

		params[0][11] = fC;
		cout << "Cost to Switch = " << params[0][11] << endl;

		params[0][12] = fDrift;
		cout << "Wind Drift = " << params[0][12] << endl;

		params[0][13] = fNumCtrls;
		cout << "Number of Controls = " << params[0][13] << endl;

		// Writing vector of parameters to file
   		io::writeToFile2D<double >(aFilename + "DomainParameters",  params);
 }
