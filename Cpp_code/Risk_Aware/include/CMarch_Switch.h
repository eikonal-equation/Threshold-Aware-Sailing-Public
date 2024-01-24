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
  * File: CMarch_Switch.h
  *
  * Author: MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu
  *
  * Description: This file contains the declarations of the class "CMarch_Switch"
  * and its member functions that compute the drift functions of the system dynamics;
  * "s-marching" step for a single gridpoint (i,j,k), a control value, and a mode;
  * best probability if switching; and Golden section search for the best steering angle.
  *
  *
  * Details of all of these functions are found in CMarch_Switch.cpp.
  *
  *============================================================================*/

#pragma once
#ifndef CMARCH_SWITCH_H
#define CMARCH_SWITCH_H

//----------------------Project specific header files---------------------------
#include "CGrid_2D.h"
#include "CENO.h"
#include "CSpecial_Treatment.h"
#include "SimpleFunctions.h"

class CMarch_Switch { // Class declaration begins here
public:

	// Constructors
	CMarch_Switch() = default;
	CMarch_Switch(int aFactor, int aSlice_factor, double aBudget, double aDiffConst, double aCost_to_Switch, double aTarget_radius,
		double aDrift, int aInterp_Order_of_ENO, double aTol_for_prob_near_1, double aTol_for_prob_near_0, int aNum_GH_pts, bool aIf_use_GoldenSS, double aTol_GoldenSS);

	// Drift function of dr =  -s(u)*cos(theta-(-1)^q*u)dt
	double r_drift(const int aMode, const double aCtrl, const double aTheta, const double aSpeed);

	// Drift function of dtheta = ((s(u)/r)*sin(theta - (-1)^q*u) +a)dt
	double theta_drift(const int aMode, const double aCtrl, const double aR, const double aTheta, const double aSpeed);

	// aI is the s slice
	// aJ is the radius index
	// aK is the theta index
	// Main function that implements the "s-marching" step for a single gridpoint (i,j,k), a control value, and a mode
	// with bi-linear interpolation
	double Marching_Step(const int aI, const int aJ, const int aK, const int aMode, const double aCtrl, const double aSpeed,
		CGrid_2D& aGrid, array_4D& aMyValfn);

	//overloded with ENO
	double Marching_Step(const int aI, const int aJ, const int aK, const int aMode, const double aCtrl, const double aSpeed,
		CENO& aENOGrid, const int aInterp_Order, array_4D& aMyValfn);

	// Main function that implements the Golden section search (GSS) on a given interval of the control value [u_lower, u_upper]
	// with bi-linear interpolation
	tuple<double, double> GoldenSection(const array<double, 2> aLower_upper_ctrl_vals, const array<double, 2> aLower_upper_speed_vals,
		const int aI, const int aJ, const int aK, const int aMode, CGrid_2D& aGrid, array_4D& aMyValfn);

	//overloded with ENO
	tuple<double, double> GoldenSection(const array<double,2> aLower_upper_ctrl_vals, const array<double,2> aLower_upper_speed_vals,
		const int aI, const int aJ, const int aK, const int aMode, CENO& aENOGrid, CGrid_2D& aGrid, const int aInterp_Order, array_4D& aMyValfn);

	// Main function that computes the expected value if we switched mode using Gauss Hermite Quadrature
	double Compute_Switch_Value(const int aI, const int aJ, const int aK, CGrid_2D& aGrid, const int aMode, array_4D& aMyValfn);

	//overloaded with ENO
	double Compute_Switch_Value(const int aI, const int aJ, const int aK, CENO& aENOgrid, const int aInterpOrder, const int aMode, array_4D& aMyValfn);

	// Main function that updates the optimal policy at a single gridpoint. I.e., whether to switch or not; if not, what is the best steering angle
	void Switch_Policy_Update(const int aI, const int aJ, const int aK,
		const double aW_switch, double& aW_best_steering, const double aBest_steer_angle, bool& aSwitch_flag,
		const int aMode, array_4D_short& aSwitch_mat, array_4D& aPolicy_mat, array_4D_short& aMyPolicy_out);

	~CMarch_Switch() {};

//member variables of the class
protected:
	int fRefinement_factor; // An input factor that would be multiplied by 100 used as number of points on each side of spatial grid. (Base case: 100*2.)
	double fBudget; // The maximum cost we are integrating up to.
	int fR; // number of points on spatial grid = fR + 1
	int fTheta; // number of theta points = fTheta + 1
	int fSlice_factor; // For grid refinement studies with respect to threshold dimension
	int fM; // number of slices
	double fR0; // starting point at r-axis
	double fTheta0; // starting point at theta-axis
	double fRmax; // largest point at r-axis
	double fThetamax; // largest point at theta-axis
	double fDr; // Delta r
	double fDtheta; // Delta theta
	double fDs; // Delta s (discretization of threshold)

	double fSigma; // Diffusion constant
	double fTau; // Time step for marching (tau = ds)
	double fRoot_Tau; //Square root of tau
	double fC; // cost to switch
	double fR_target; //radius of target
	double fDrift; //drift ("a")


	int fNum_slices; //Number of slices we keep in the RAM
	int fInterp_Order; //The interpolation order of ENO
	double fTol_1; // how close is close to 1?
	double fTol_zero;// how close is close to 0?

	int fNum_GH_pts; //The number of points of the Gauss Hermite Quadrature
	bool fIf_use_GoldenSS; //boolean variable that decides if we use Golden section search (GSS) or not
	double fTol_GoldenSS; //Tolerance of GSS

	array_1D fHermite_root; //roots of nth-order Hermite polynomial
	array_1D fHermite_weight; //weights of nth-order Gauss-Hermite Quadrature
};

// drift function of dr =  -s(u)*cos(theta-(-1)^qu)dt
inline double CMarch_Switch::r_drift(const int aMode, const double aCtrl, const double aTheta, const double aSpeed)
{
	return -1.0 * aSpeed * cos(aTheta - pow(-1, double(aMode)) * aCtrl);
}

// drift function of dtheta = ((s(u)/r)*sin(theta - (-1)^qu) +a)dt
inline double CMarch_Switch::theta_drift(const int aMode, const double aCtrl, const double aR, const double aTheta, const double aSpeed)
{
	return (aSpeed / aR) * sin(aTheta - pow(-1, double(aMode)) * aCtrl) + fDrift;
}

#endif // !CMARCH_SWITCH_H

