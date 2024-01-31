# Risk-aware stochastic control of a sailboat
This repository contains the source code used to generate all examples presented in "Risk-aware stochastic control of a sailboat" manuscript by MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu, and Alexander Vladimirsky.

# License #
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

## Abstract of our manuscript ##
Sailboat path-planning is a natural hybrid control problem (due to continuous steering and occasional "tack-switching" maneuvers), 
with the actual path-to-target greatly affected by stochastically evolving wind conditions. 
Previous studies have focused on finding risk-neutral policies that minimize the expected time of arrival. 
In contrast, we present a robust control approach, which maximizes the probability of arriving before a specified deadline/threshold. 
Our numerical method recovers the optimal risk-aware (and threshold-specific) policies for all initial sailboat positions and a broad
range of thresholds simultaneously. This is accomplished by solving two quasi-variational inequalities based on second-order
Hamilton-Jacobi-Bellman (HJB) PDEs with degenerate parabolicity. Monte-Carlo simulations show that risk-awareness
in sailing is particularly useful when a carefully calculated bet on the evolving wind direction might yield a reduction in the
number of tack-switches.

# Manuscript #
The "Risk-aware stochastic control of a sailboat" manuscript can be found [here](https://arxiv.org/abs/2309.13436).

# Contributions & Acknowledgements # 
  * The problem statement and the numerical scheme were developed by MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu, and Alexander Vladimirsky.
  * The implementation was carried out by MingYi Wang, Natasha Patnaik, Anne Somalwar, and Jingyi Wu.
  * The manuscript was written by MingYi Wang, Natasha Patnaik, Anne Somalwar, Jingyi Wu, and Alexander Vladimirsky.
  * The authors acknowledge Cole Miles, Roberto Ferretti, and Adriano Festa for inspiring our work.
# Instructions #
  
## Requirements: ## 
* The C++ code requires the users to install the "[Boost](https://www.boost.org/)" library (external). 
    * We install it directly under the directory `Cpp_code`. Please feel free to install it anywhere you want on your end but you would have to modify its path in the ***Makefile*** described below.

* The CDFs and the deterministic-optimal policy are generated with Matlab code.

## Running the C++ Code: ##
The following instructions explain how to run the Solver for our risk(threshold)-aware optimal value functions/policies using the ***Makefile***. 

To change compilers, edit `CC=` by putting your desired compiler after the `=` sign. The default compiler is set to `g++`. 

To update the path/version of the "Boost" library you've installed, edit `BOOSTPATH` by putting your own path/version after the `=` sign. The current path/version is `./boost_1_79_0`.

To compile the code, type `make main` in the terminal in this folder. 
* Please note that you need to enter the `Risk-Neutral` folder to run the PDE solver for the risk-neutral value function/policy, and the `Risk-Aware` folder to run the PDE solver for the risk-aware value function/policy, respectively.
* In either case, after compilation, type `./main` to run the respective PDE solver.
  
To delete the executable and the data files, type `make clean`.

## Running the Matlab Code: ##
To generate a CDF with risk(threshold)-aware optimal policies:
  * `Sailing_Risk_Aware_CDF.m`
      * Produces a CDF $y(s) = {\mathbb{P}}(T \le s)$ that measures the probability of keeping the accumulative time-to-target $T$ under any threshold value $s$ but maximized at a given initial threshold/budget $\hat{s}$ using threshold-aware policies. It will produce a plot of the CDF and a plot of a sample path at the end of execution. 
      * Note: It requires folders containing data matrices of the risk-neutral policy (optimal steering angle/switchgrid) and data matrices of the risk-aware policy (optimal steering angle/switchgrid) as inputs. 
          * Please make sure the folder names under the `output` folder in both cases are the **same** with the **same** wind drift and diffusion constants.
          * In the simplest case, you can just run the C++ code without creating another folder under `output` and just define the argument as `DataFile_dir = ''`.
     
To generate a CDF with the risk-neutral policy:
   * `Sailing_Risk_Neutral_CDF.m`
      * Produces a CDF $y(s) = {\mathbb{P}}(T \le s)$ that measures the probability of success where the cumulative time-to-target $T$ is within any positive threshold value $s$ using the risk-neutral policy. It will produce a plot of the CDF and a plot of a sample path at the end of execution. 
      * Note: It requires a folder (under `output`) containing data matrices of the risk-neutral policy (optimal steering angle/switchgrid).
      * In the simplest case, you can just run the C++ code without creating another folder under `output` and just define the argument as `DataFile_dir = ''`.

Demonstration:
  * `demo.m`
      * This script is solely for demonstration. The user need to compile and run the C++ code first before running this script. (Please start with a coarse grid.)
 
