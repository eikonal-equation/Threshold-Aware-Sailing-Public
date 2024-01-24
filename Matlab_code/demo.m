%% "Risk-aware stochastic control of a sailboat" Numerical Experiments Demo file
%
% This script is usued to demonstrate how to parse date from running our
% C++ code into Matlab and run our Matlab code to generate ECDFs
%
% The user need to run the C++ code themselves to obtain the data files
% before running this "Demo" script.
%
% Author: MingYi Wang, Cornell University
% Last modified: 01/2024
%
clear all;
close all;
%% Generating CDFs with the deterministic-optimal policy
%Specify the parameters
DataFile_dir = ''
sample_size = 1e3
time_step = 0.005
r0_indx = 381
t0_indx = 1
choice = 'nearest'
% Call our function to generate the CDF
Xcost_stationary = Sailing_Risk_Neutral_CDF(DataFile_dir,sample_size,time_step,choice,r0_indx,t0_indx);

%% Generating CDFs with threshold-aware policies
% Using the same initial tumor configuration, the same policy determination
% startegy, and the same sample size as above.
%
DataFile_dir = ''
% Set up the initial budget sbar
Initial_budget = 44
Max_cost_from_PDE_solver = 55
% Call our function to generate the CDF
Xcost_thres = Sailing_Risk_Aware_CDF(DataFile_dir,sample_size,time_step,choice,...
 r0_indx,t0_indx,Initial_budget,Max_cost_from_PDE_solver);