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
  * File: CMarch_Switch.cpp
  *
  * Author: MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu
  *
  * Description: This file contains the actual implementations of the member functions
  * that compute the drift functions of the system dynamics;
  * "s-marching" step for a single gridpoint (i,j,k), a control value, and a mode;
  * best probability if switching; and Golden section search for the best steering angle.
  *
  *
  *============================================================================*/

//----------------------Project specific header files---------------------------
#include "CMarch_Switch.h"

//Class constructor
CMarch_Switch::CMarch_Switch(int aRefinement_factor, int aSlice_factor, double aBudget, double aDiffConst, double aCost_to_Switch, double aTarget_radius,
	double aDrift, int aInterp_Order_of_ENO, double aTol_for_prob_near_1, double aTol_for_prob_near_0, int aNum_GH_pts, bool aIf_use_GoldenSS, double aTol_GoldenSS) {

	fRefinement_factor = aRefinement_factor; // An input factor that would be multiplied by 100 used as number of points on each side of spatial grid.  (base case: 201)
	fBudget = aBudget; // The max cost we are integrating up to.
	fR = 200 * fRefinement_factor; // index of the last r point
	fTheta = 200 * fRefinement_factor; // index of the last theta point
	fSlice_factor = aSlice_factor; // For grid refinement studies with respect to threshold dimension
	fM = int(aSlice_factor * fBudget * fRefinement_factor); // index of the last slice
	fR0 = 0; // starting point at r-axis
	fTheta0 = 0; // starting point at theta-axis
	fRmax = 2.0; // largest point at r-axis; default Rmax = 2
	fThetamax = 2 * PI; // largest point at theta-axis
	fDr = (fRmax - fR0) / (double)fR; // Delta r
	fDtheta = (fThetamax - fTheta0) / (double)fTheta; // Delta theta
	fDs = fBudget / fM; // Delta s (discretization of cost)
	fSigma = aDiffConst;  // Diffusion constant
	fTau = fDs; // time step for marching (tau = ds)
	fRoot_Tau = sqrt(fTau);
	fR_target = aTarget_radius; //radius of target set
	fDrift = aDrift; //drift ("a")
	fC = aCost_to_Switch; // cost to switch
	fNum_slices = int(fC / fDs) + 2; //Number of slices we keep in the RAM

	fNum_GH_pts = aNum_GH_pts; //The number of points of the Gauss Hermite Quadrature

	fHermite_root.resize(boost::extents[fNum_GH_pts]);  //roots of nth-order Hermite polynomial
	fHermite_weight.resize(boost::extents[fNum_GH_pts]); //weights of nth-order Gauss-Hermite Quadrature

	if (fNum_GH_pts == 3) { //3rd-order Hermite polynomial
		fHermite_root[0] = -1.0 * sqrt(3.0 / 2.0);
		fHermite_root[1] = 0.0;
		fHermite_root[2] = sqrt(3.0 / 2.0);
		//Gauss-Hermite Quadrature weights (absorbed 1/sqrt(pi) into it)
		fHermite_weight[0] = 1.0 / 6.0;
		fHermite_weight[1] = 2.0 / 3.0;
		fHermite_weight[2] = 1.0 / 6.0;
	}
	else if (fNum_GH_pts == 5) { //5th-order Hermite polynomial
		fHermite_root[0] = -1.0 * sqrt(5.0 / 2.0 + sqrt(5.0 / 2.0));
		fHermite_root[1] = -1.0 * sqrt(5.0 / 2.0 - sqrt(5.0 / 2.0));
		fHermite_root[2] = 0.0;
		fHermite_root[3] = sqrt(5.0 / 2.0 - sqrt(5.0 / 2.0));
		fHermite_root[4] = sqrt(5.0 / 2.0 + sqrt(5.0 / 2.0));
		//Gauss-Hermite Quadrature weights (absorbed 1/sqrt(pi) into it)
		fHermite_weight[0] = 0.011257411327721;
		fHermite_weight[1] = 0.222075922005613;
		fHermite_weight[2] = 0.533333333333333;
		fHermite_weight[3] = 0.222075922005613;
		fHermite_weight[4] = 0.011257411327721;
	}
	else if (fNum_GH_pts == 7) { //7th-order Hermite polynomial
		fHermite_root[0] = -2.651961356835236;
		fHermite_root[1] = -1.673551628767473;
		fHermite_root[2] = -0.816287882858964;
		fHermite_root[3] = 0.0;
		fHermite_root[4] = 0.816287882858964;
		fHermite_root[5] = 1.673551628767473;
		fHermite_root[6] = 2.651961356835236;
		//Gauss-Hermite Quadrature weights (absorbed 1/sqrt(pi) into it)
		fHermite_weight[0] = 0.000548268855972;
		fHermite_weight[1] = 0.030757123967586;
		fHermite_weight[2] = 0.240123178605013;
		fHermite_weight[3] = 0.457142857142857;
		fHermite_weight[4] = 0.240123178605013;
		fHermite_weight[5] = 0.030757123967586;
		fHermite_weight[6] = 0.000548268855972;
	}

	fTol_1 = aTol_for_prob_near_1; // how close is close to 1?
	fTol_zero = aTol_for_prob_near_0; // how close is close to 0?
	fInterp_Order = aInterp_Order_of_ENO; //The interpolation order of ENO
	fIf_use_GoldenSS = aIf_use_GoldenSS; //boolean variable that decides if we use GSS or not
	fTol_GoldenSS = aTol_GoldenSS; //Tolerance of GSS
}



// This function implements the "marching step" of the simplified
// problem. For each (aRowIndex, aColIndex)-pair, given an interpolation
// shceme, we compute the feet of two sample points:
// y+ = fx0 + f_i(aCtrl,aMode)*fTau + fEps*sqrt(fTau);
// y- = fx0 + f_i(aCtrl,aMode)*fTau - fEps*sqrt(fTau);
// And compute the interpolated value function at those two points.
// We compute their average and update the value function/optimal
// policy only if it gives a higher value
double CMarch_Switch::Marching_Step(const int aI, const int aJ, const int aK, const int aMode, const double aCtrl, const double aSpeed,
	CGrid_2D& aGrid, array_4D& aMyValfn)
{
	Where_to_interp HERE = Where_to_interp::marching;
	double r_0 = fR0 + double(aJ) * fDr;
	double theta_0 = fTheta0 + double(aK) * fDtheta;

	// There is NO stochasticity with respect to r state space
	double rloc = r_0 + fTau * r_drift(aMode, aCtrl, theta_0, aSpeed);

	// There is stochasticity with respect to theta state space
	array<double, 2> theta_loc{ 0.0,0.0 }; //storing the two possible locations of theta_loc, index [0] is plus while [1] is minus
	for (int ii = 0; ii < 2; ii++) {
		theta_loc[ii] = mod_2pi(theta_0 + fTau * theta_drift(aMode, aCtrl, r_0, theta_0, aSpeed) + pow(-1, ii) * (fSigma * fRoot_Tau));
	}


	double w_possible_avg = 0.0;
	//if the new location takes us to the target
	if (rloc <= fR_target)
	{
		w_possible_avg = 1.0;
	}
	else if (rloc > fRmax) {
		w_possible_avg = 0.0;
	}
	else {
		int right_index_r = fR + 1;
		if (rloc == fRmax) {
			right_index_r = fR;
		}
		else {
			right_index_r = min(find_index(rloc, fDr, fR0, fRmax), fR);
		}
		array<int, 2> right_index_theta{ -100,-100 }; //store the right index of theta_loc, [0] is "+" while [1] is "-"
		for (int ii = 0; ii < 2; ii++) {
			right_index_theta[ii] = find_index(theta_loc[ii], fDtheta, fTheta0, fThetamax);
		}
		array<double, 4> v_stencil;
		double right_Rloc, right_ThetaLoc;
		// avg of value function vals computed with theta+,theta-
		for (int ii = 0; ii < 2; ii++) {
			std::tie(v_stencil, right_ThetaLoc, right_Rloc) = aGrid.Stencil_for_Bilinear_Interp(aMyValfn, aMode, aI, right_index_r, right_index_theta[ii], HERE);
			w_possible_avg += min(aGrid.Bilinear_Interp(v_stencil, right_ThetaLoc, right_Rloc, theta_loc[ii], rloc, fDtheta, fDr), 1.0) / 2.0;
		}
	}
	return w_possible_avg;
}

//overloaded with ENO
double CMarch_Switch::Marching_Step(const int aI, const int aJ, const int aK, const int aMode, const double aCtrl, const double aSpeed,
	CENO& aENOGrid, const int aInterp_Order, array_4D& aMyValfn)
{
	double w_possible_avg = 0.0;
	Where_to_interp HERE = Where_to_interp::marching;
	double r_0 = fR0 + double(aJ) * fDr;
	double theta_0 = fTheta0 + double(aK) * fDtheta;

	// There is NO stochasticity with respect to r state space
	double rloc = r_0 + fTau * r_drift(aMode, aCtrl, theta_0, aSpeed);

	// There is stochasticity with respect to theta state space
	array<double, 2> theta_loc{ 0.0,0.0 }; //storing the two possible locations of theta_loc, index [0] is plus while [1] is minus
	for (int ii = 0; ii < 2; ii++) {
		theta_loc[ii] = mod_2pi(theta_0 + fTau * theta_drift(aMode, aCtrl, r_0, theta_0, aSpeed) + pow(-1, ii) * (fSigma * fRoot_Tau));
	}

	//if the new location takes us to the target
	if (rloc <= fR_target)
	{
		w_possible_avg = 1.0;
	}
	else if (rloc > fRmax) {
		w_possible_avg = 0.0;
	}
	else {
		int right_index_r = min(find_index(rloc, fDr, fR0, fRmax), fR);

		array<int, 2> right_index_theta{ -100,-100 }; //store the right index of theta_loc, [0] is "+" while [1] is "-"

		for (int ii = 0; ii < 2; ii++) {
			right_index_theta[ii] = find_index(theta_loc[ii], fDtheta, fTheta0, fThetamax);
		}

		for (int ii = 0; ii < 2; ii++) {
			//Construct computing stencil for ENO
			array_2D ENO_stencil = aENOGrid.Matrix_for_ENO_Interp_RTheta(aMyValfn, aInterp_Order, aMode, aI, right_index_r, right_index_theta[ii]);

			//ENO interpolation
			if (aInterp_Order == 2) { //quadatric interpolation
				w_possible_avg += min(aENOGrid.ENO_quad_interp_RTheta(ENO_stencil, right_index_r, right_index_theta[ii], rloc, theta_loc[ii]), 1.0) / 2.0;
			}
			else if (aInterp_Order == 3) {//cubic interpolation
				w_possible_avg += min(aENOGrid.ENO_cubic_interp_RTheta(ENO_stencil, right_index_r, right_index_theta[ii], rloc, theta_loc[ii]), 1.0) / 2.0;
			}
		}
	}
	return w_possible_avg;
}


// Adapted from https ://github.com/bilginmu/LineSearchAlgorithms/blob/master/c%2B%2B/golden_search_method.cpp
//aCtrl_lower is the control before the optimal base control, aCtrl_upper is the control after
tuple<double, double> CMarch_Switch::GoldenSection(const array<double, 2> aLower_upper_ctrl_vals, const array<double, 2> aLower_upper_speed_vals,
	const int aI, const int aJ, const int aK, const int aMode, CGrid_2D& aGrid, array_4D& aMyValfn)
{
	//lower value at index 0 and upper value at index 1
	double error = aLower_upper_ctrl_vals[1] - aLower_upper_ctrl_vals[0];
	double aPrev = aLower_upper_ctrl_vals[0];
	double bPrev = aLower_upper_ctrl_vals[1];

	double aNew, bNew, val_at_aNew, val_at_bNew;
	while (error > fTol_GoldenSS) {
		aNew = aPrev + (bPrev - aPrev) * inv_GR_square; // New lower limit
		bNew = bPrev - (bPrev - aPrev) * inv_GR_square;

		//Need to interpolate to find the speed first
		double aSpeed_aNew = aGrid.Linear_Interp(aLower_upper_speed_vals, aLower_upper_ctrl_vals, aNew);
		val_at_bNew = Marching_Step(aI, aJ, aK, aMode, bNew, aSpeed_aNew, aGrid, aMyValfn);

		double aSpeed_bNew = aGrid.Linear_Interp(aLower_upper_speed_vals, aLower_upper_ctrl_vals, bNew);
		val_at_aNew = Marching_Step(aI, aJ, aK, aMode, aNew, aSpeed_bNew, aGrid, aMyValfn);

		// if val @ aNew > val @ bNew -> max exists in the interval [aPrev, bNew], so we update bPrev as bNew
		if (val_at_bNew < val_at_aNew) {
			bPrev = bNew;
		}
		else if (val_at_bNew >= val_at_aNew) {
			aPrev = aNew;
		}
		error = fabs(bPrev - aPrev);
	}

	//returning maximum of the two endpoints of the interval
	if (val_at_bNew < val_at_aNew) {
		return make_tuple(val_at_aNew, aNew);
	}
	else {
		return make_tuple(val_at_bNew, bNew);
	}
}

//Overloaded with ENO
tuple<double, double> CMarch_Switch::GoldenSection(const array<double, 2> aLower_upper_ctrl_vals, const array<double, 2> aLower_upper_speed_vals,
	const int aI, const int aJ, const int aK, const int aMode, CENO& aENOGrid, CGrid_2D& aGrid, const int aInterp_Order, array_4D& aMyValfn) {

	//lower value at index 0 and upper value at index 1
	double error = aLower_upper_ctrl_vals[1] - aLower_upper_ctrl_vals[0];
	double aPrev = aLower_upper_ctrl_vals[0];
	double bPrev = aLower_upper_ctrl_vals[1];

	double aNew, bNew, val_at_aNew, val_at_bNew;
	int ii = 1;
	while (error > fTol_GoldenSS) {
		aNew = aPrev + (bPrev - aPrev) * inv_GR_square; // New lower limit
		bNew = bPrev - (bPrev - aPrev) * inv_GR_square;

		//Need to interpolate to find the speed first
		double aSpeed_bNew = aGrid.Linear_Interp(aLower_upper_speed_vals, aLower_upper_ctrl_vals, bNew);

		val_at_bNew = Marching_Step(aI, aJ, aK, aMode, bNew, aSpeed_bNew, aENOGrid, aInterp_Order, aMyValfn);

		double aSpeed_aNew = aGrid.Linear_Interp(aLower_upper_speed_vals, aLower_upper_ctrl_vals, aNew);
		val_at_aNew = Marching_Step(aI, aJ, aK, aMode, aNew, aSpeed_aNew, aENOGrid, aInterp_Order, aMyValfn);

        // if val @ aNew > val @ bNew -> max exists in the interval [aPrev, bNew], so we update bPrev as bNew
		if (val_at_bNew < val_at_aNew) {
			bPrev = bNew;
		}
		else if (val_at_bNew >= val_at_aNew) {
			aPrev = aNew;
		}

		error = fabs(bPrev - aPrev);
	}
	//returning maximum of the two endpoints of the interval
	if (val_at_bNew < val_at_aNew) {
		return make_tuple(val_at_aNew, aNew);
	}
	else {
		return make_tuple(val_at_bNew, bNew);
	}
}


// Main function that computes the expected value if we switched mode using Gauss Hermite Quadrature
double CMarch_Switch::Compute_Switch_Value(const int aI, const int aJ, const int aK, CGrid_2D& aGrid, const int aMode, array_4D& aMyValfn) {

	Where_to_interp HERE = Where_to_interp::switching;
	double w_switch = 0.0;

	//current theta location
	double theta_current = fTheta0 + double(aK) * fDtheta;

	array_1D eta_list(boost::extents[fNum_GH_pts]);

	for (int m = 0; m < fNum_GH_pts; m++)
	{
		eta_list[m] = mod_2pi(theta_current + fDrift * fC + fSigma * sqrt(2 * fC) * fHermite_root[m]);
	}

	//find bottom s-slice after switch
	double s_new_loc = double(aI) * fDs - fC;
	double s_max_current = double(aI) * fDs;
	double s_min_current = s_max_current - double(fNum_slices - 1) * fDs;
	int s_top_ind = find_index(s_new_loc, fDs, s_min_current, s_max_current);

	//Computing value of switching
	for (int n = 0; n < fNum_GH_pts; n++)
	{
		int theta_new_ind = find_index(eta_list[n], fDtheta, fTheta0, fThetamax);

		array<double, 4> Switching_stencil{ 0.0,0.0,0.0,0.0 };

		double right_ThetaLoc, top_Sloc;

		std::tie(Switching_stencil, right_ThetaLoc, top_Sloc) = aGrid.Stencil_for_Bilinear_Interp(aMyValfn, 3 - aMode, s_top_ind, aJ, theta_new_ind, HERE);

		top_Sloc += s_min_current; //need to add s_min back since we only store limited number of slices

		double switch_value_possible = aGrid.Bilinear_Interp(Switching_stencil, right_ThetaLoc, top_Sloc, eta_list[n], s_new_loc, fDtheta, fDs);

		w_switch += fHermite_weight[n] * switch_value_possible;

	}
	return min(w_switch, 1.0);
}

//overloaded with ENO
double CMarch_Switch::Compute_Switch_Value(const int aI, const int aJ, const int aK, CENO& aENOgrid, const int aInterpOrder, const int aMode, array_4D& aMyValfn)
{
	//Initialization
	double w_switch = 0.0;
	array_1D eta_list(boost::extents[fNum_GH_pts]);

	//current theta location
	double theta_current = fTheta0 + double(aK) * fDtheta;

	for (int m = 0; m < fNum_GH_pts; m++)
	{
		eta_list[m] = mod_2pi(theta_current + fDrift * fC + fSigma * sqrt(2 * fC) * fHermite_root[m]);
	}

	//find bottom s-slice after switch
	double s_new_loc = double(aI) * fDs - fC;
	int s_new_slice_num = 1; //Based on our memory saving structure, we will always traverse to the 2nd s-slice (whose index is "1")

	//Computing value of switching
	for (int n = 0; n < fNum_GH_pts; n++)
	{
		double switch_value_possible = -10;

		int theta_new_ind = find_index(eta_list[n], fDtheta, fTheta0, fThetamax);
		if (aInterpOrder == 2) {
			const int aSize = 4;
			array<double, aSize> aENO_switch_stencil = aENOgrid.Array_for_ENO_Interp_Theta<aSize>(aMyValfn, aInterpOrder, 3 - aMode, s_new_slice_num, theta_new_ind, aJ);
			switch_value_possible = aENOgrid.ENO_quad_interp_Thetadim(aENO_switch_stencil, theta_new_ind, eta_list[n], fDtheta, fTheta0, fTheta);
		}
		else if (aInterpOrder == 3) {
			const int aSize = 6;
			array<double, aSize> aENO_switch_stencil = aENOgrid.Array_for_ENO_Interp_Theta<aSize>(aMyValfn, aInterpOrder, 3 - aMode, s_new_slice_num, theta_new_ind, aJ);
			switch_value_possible = aENOgrid.ENO_cubic_interp_Thetadim(aENO_switch_stencil, theta_new_ind, eta_list[n], fDtheta, fTheta0, fTheta);
		}

		w_switch += fHermite_weight[n] * switch_value_possible;
	}
	return min(w_switch, 1.0);
}

// Main function that updates the optimal policy at a single gridpoint. I.e., whether to switch or not; if not, what is the best steering angle
void CMarch_Switch::Switch_Policy_Update(const int aI, const int aJ, const int aK,
	const double aW_switch, double& aW_best_steering, const double aBest_steer_angle, bool& aSwitch_flag, const int aMode,
	array_4D_short& aMySwitch, array_4D& aMyPolicy, array_4D_short& aMyPolicy_out)
{
	if ((aW_switch > aW_best_steering) || ((aW_switch <= aW_best_steering) && (aW_best_steering - aW_switch < 0.65 * fTol_1)))
	{
		aW_best_steering = aW_switch;
		aSwitch_flag = true;
		aMySwitch[aMode - 1][aI][aJ][aK] = 1;
	}
	else {
		aMySwitch[aMode - 1][aI][aJ][aK] = 0;
	}
	aMyPolicy[aMode - 1][aI][aJ][aK] = aBest_steer_angle;
	//for storage
	aMyPolicy_out[aMode - 1][aI][aJ][aK] = (short)round(aBest_steer_angle * 1e4);
}

