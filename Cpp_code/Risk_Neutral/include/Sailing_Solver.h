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
  * File: Sailing_Solver.h
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


#pragma once
#ifndef Sailing_Solver_H
#define Sailing_Solver_H

//----------------------Project specific header files---------------------------
#include "CGrid_2D.h"

//-------------------------Libraires-------------------------------------------
#include<algorithm>
#include<numeric>
#include<tuple>

//---------------------------Definitions---------------------------------------
typedef boost::multi_array<double , 1> array_1D;
typedef boost::multi_array<double , 2> array_2D;
typedef boost::multi_array<double , 3> array_3D;
typedef boost::multi_array<double , 4> array_4D;
enum class speed_method { INTERPOLATE, CEIL, SETCONTROL };

//---------------------------Main Class----------------------------------------
class Sailing_Solver
{
public:
	// constructor(s)
	Sailing_Solver() = default;
	Sailing_Solver(int aFactor, double aDiffConst, double aTarget_radius, double aDrift, int aNumCtrls, double aTol_converge, string aSpeed_File, speed_method aSpeed_choice, double aCost_to_switch);

	// These functions compute the time derivatives in the differential equation, missing a factor of the boat speed
	// r - dist, theta - angle, q - tack, u - steering angle
	double  dr(const double  aTheta, const double  aQ, const double  aU);
	double  dtheta(const double  aR, const double  aTheta, const double  aQ, const double  aU);
	//Determine the size of the time step
	double  dt_disc(const double  aSr, const double  aStheta, const bool aRowwise);

	// Finds the first index s.t. xloc <= aArray
	int find_index(const double  aXloc, const double  aDeltaX, const double  aX0, const double  aXmax);

	// Loads in speed profile data
	void load_speed_data(string aSpeed_File);

	// Performs operations mode 2pi
	double  mod_2pi(const double  aValue);

	// Initialization of matrices which store value function values, optimal policies, and switchgrids
	void InitializeMat();

	// Write Matrices to file
	void writeMatstoFile();

	// Initialization of control and speed profile
	void Initialize_Control_and_Speed_List(array_1D& aControl_List, array_1D& aSpeed_List, CGrid_2D& aGrid);

	//This function implements the binary search algorithm
	int binary_search(const vector<double > aVector, const double  aXloc);

	// Main function that computes the expected value if we switched mode using Gauss Hermite Quadrature
	double  gauss_hermite_switch(int aIr, int aItheta, int aQ, CGrid_2D& aGrid);

	//This is the main solver that uses (Row-wise) Value Iterations to compute the value function and the corresponding optimal policy
	void solve_value_sweep(const bool aRowwise);


	/** Grid parameters */
	double  getThetaGridSize() const; // number of theta points
	double  getThetaMax() const; // largest point at theta-axis
	double  getRadiusGridSize() const; // // number of r points on spatial grid
	double  getRadiusMax() const; // largest point at r-axis
	double  getGFactor() const; // for grid refinement studies
	double  getTargetSize() const; // radius of target region

	// Writing Domain Parameters to file
	void writeDomainToFile(const std::string aFilename) const;

//member variables of the class
private:
	int fMulti_factor; // An input factor that would be multiplied by 100 used as number of points on each side of spatial grid. (Base case: 100*2)
	double  fPi; // constant for pi
	int fR; // number of points on spatial grid = fR + 1
	int fTheta; // number of theta points = fTheta + 1
	double  fR0; // starting point at r-axis
	double  fTheta0; // starting point at theta-axis
	double  fRmax; // largest point at r-axis
	double  fThetamax; // largest point at theta-axis
	double  fDr; // Delta r
	double  fDtheta; // Delta theta
	double  fSigma; // Diffusion constant
	double  fC; // cost to switch
	double  fTol; // Value Iterations (VI) convergence tolerance
	double  fR_target; // radius of target
	double  fDrift; // drift ("a")
	int fNumCtrls; // discretization of steering control variable

	vector<vector<double >> fSpeed_Data; // speed profile data
	vector<double > fAngles;  //list of base control values (steering angles)
	string fSpeed_File;  //filename of the speed profile
	double fMax_speed; //max speed of the sailboat
	int fMax_iters; // Maximum iteration numbers for VI
	double  fLambda; // discounting factor (always set to 0 for all of our numerical experiments)
	array_3D fMyValfn; // 3D array for the value function
	array_3D fMyPolicy; // 3D array for the optimal steering angles
	array_3D fMySwitch; // 3D array for the switchgrid

	speed_method fSpeed_Method; //method of interpreting/discretizing the speed profile
	int fWind_Speed_Index; //wind speed index in the .pol file

	array<double, 3> fHermite_root; //roots of 3rd-order Hermite polynomial
	array<double, 3> fHermite_weight;  //weights of 3rd-order Gauss-Hermite Quadrature
};


inline double  Sailing_Solver::getThetaGridSize() const{
	return fTheta;
}

inline double  Sailing_Solver::getThetaMax() const{
	return fThetamax;
}

inline double  Sailing_Solver::getRadiusGridSize() const{
	return 	fR;
}

inline double  Sailing_Solver::getRadiusMax() const{
	return fRmax;
}

inline double  Sailing_Solver::getGFactor() const{
	return fMulti_factor;
}

inline double  Sailing_Solver::getTargetSize() const{
	return fR_target;
}

inline int Sailing_Solver::find_index(const double  aXloc, const double  aDeltaX, const double  aX0, const double  aXmax)
{

	//assert(aXloc <= aXmax);
	//assert(aXloc >= aX0);

	int right = int(ceil((aXloc - aX0)/aDeltaX));
	return right;
}

inline double  Sailing_Solver::mod_2pi(const double  aValue){
	return aValue - 2.0*fPi*floor(aValue / (2.0*fPi));
}

#endif // !SEMICAUSAL_THRESHOLD_SOLVER_H
