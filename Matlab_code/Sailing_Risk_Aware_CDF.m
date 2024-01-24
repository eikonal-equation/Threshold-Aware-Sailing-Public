function Xcost = Sailing_Risk_Aware_CDF(DataFile_dir,sample_size,time_step,choice,...
 r0_indx,t0_indx,Initial_budget,Max_cost_from_PDE_solver)
%This function computes the CDF y = Pr(T <= s) that measures the
%probability of success where the cumulative time-to-target T is 
%within any positive threshold value "s" using the Risk-Aware policy
%computed from Wang et al.
%https://arxiv.org/abs/2309.13436
%with a specification of a targeting initial budget "shat".
%
%Author: MingYi Wang, Cornell University
%Last Modified: 01/2024
%
%DataFile_dir (input): The directory under the folder "output" where the
%                       user stores their data files (risk-aware)
%sample_size (input): sample size of Monte Carlo simulations
%time_step (input): time step "dt" used in the Euler-Maruyama scheme
%choice (input): choice of the policy determination strategy, 
%                4 options:  (i) 'nearest', 
%                            (ii) 'conservative'; 
%                            (iii) 'aggressive'; 
%                            (iv) 'majority'.

%r0_indx (input): the index (i) of the sailboat's starting r value on the
%                 spatial grid
%t0_indx (input): the index (j) of the sailboat's starting theta value on the
%                 spatial grid
%                 E.g., (r_i, theta_j); in Matlab, i, j both start from 1.
%Initial_budget (input): The user chosen starting threhold value "shat"
%Max_cost_from_PDE_solver (input) : The maximum cost we integrating up to from the
%                                   PDE solver, i.e., "sbar" in the paper
%
%Xcost (output): Array of the accumulative (random) time-to-target of each sample
%
%% Read in parameter data
% Define file names
filename_parameters = ['../Risk_Aware/output/',DataFile_dir,'/PathTracer_DomainParameters'];
filename_controls = ['../Risk_Aware/output/',DataFile_dir,'/ControlList.dat'];
filename_base_controls = ['../Risk_Aware/output/',DataFile_dir,'/Base_ControlList.dat'];
filename_speeds = ['../Risk_Aware/output/',DataFile_dir,'/SpeedList.dat'];
filename_base_speeds = ['../Risk_Aware/output/',DataFile_dir,'/Base_SpeedList.dat'];
filename_policy = ['../Risk_Aware/output/',DataFile_dir,'/test_Policies.dat'];
filename_switch = ['../Risk_Aware/output/',DataFile_dir,'/test_Switchgrid.dat'];
filename_valuefn = ['../Risk_Aware/output/',DataFile_dir,'/test_ValueFunc.dat'];

filename_policy_stat = ['../Risk_Neutral/output/',DataFile_dir,'/PathTracer_Policies.dat'];
filename_switch_stat = ['../Risk_Neutral/output/',DataFile_dir,'/PathTracer_Switchgrid.dat'];

% Getting Parameters
precision = 'double';
uFile_params = fopen(filename_parameters);
u_params = fread(uFile_params, 18, precision);

u_params = u_params';

theta_num = u_params(1) + 1;

theta_max = u_params(2); % maximum theta in domain = 2pi
r_num = u_params(3) + 1;
r_max = u_params(4); % maximum radius in domain
g_factor = u_params(5);
target_radius = u_params(6); % radius of target region
num_slices = u_params(7);
Max_budget = u_params(8);
r_min = u_params(9);
theta_min = u_params(10);
dr = u_params(11);
dtheta= u_params(12);
ds = u_params(13);
sig = u_params(14);
tau = u_params(15);
switch_cost = u_params(16);
drift = u_params(17);
num_ctrl = u_params(18);
target_indx = floor(target_radius/dr)+1;

% Getting Controls and Corresponding Speeds
uFile_cntrls = fopen(filename_controls);
controls = fread(uFile_cntrls, num_ctrl, precision);
uFile_speeds = fopen(filename_speeds);
speeds = fread(uFile_speeds, num_ctrl, precision);
fclose(uFile_params);

base_num_ctrl = 24;
% Getting COre Controls and Corresponding Speeds from the data
uFile_cntrls = fopen(filename_base_controls);
base_controls = fread(uFile_cntrls, base_num_ctrl, precision);
uFile_speeds = fopen(filename_base_speeds);
base_speeds = fread(uFile_speeds, base_num_ctrl, precision);
fclose(uFile_params);

base_controls = round(base_controls*1e4)/1e4;
%% Creating memory maps
num_slices = num_slices / g_factor;
mPolicy_mat_thres = memmapfile(filename_policy,'Format',{'int16',[theta_num, r_num,2, num_slices+1],'uu'});
mSwitchgrid_thres = memmapfile(filename_switch,'Format',{'int16',[theta_num, r_num,2, num_slices+1],'ss'});
%% Read in Optimal Policy Data and Switchgrids
tic
precision = 'double';
% ------ Read in optimal policies, store in a matrix
vFile = fopen(filename_policy_stat);
Policy_mat_stat = fread(vFile, 2*theta_num*r_num, precision);
fclose(vFile);
% Reshaping
Policy_mat_stat= reshape(Policy_mat_stat,[theta_num, r_num, 2]);

for mode = 1:2
    temp = Policy_mat_stat(:,:,mode);
    temp(end, :) = temp(1,:);
    Policy_mat_stat(:,:,mode) = temp;
end

% ------ Read in switch grids, store in a matrix
sFile = fopen(filename_switch_stat);
Switchgrid_stat = fread(sFile, 2*theta_num*r_num, precision);
fclose(sFile);
% Reshaping
Switchgrid_stat = reshape(Switchgrid_stat,[theta_num, r_num, 2]);
for mode = 1:2
    temp = Switchgrid_stat(:,:,mode);
    temp(end,:) = temp(1,:);
    Switchgrid_stat(:,:,mode) = temp;

end
load_time = toc


%% Monte Carlo Simulation Initialization
dt = time_step;
root_dt = sqrt(dt);

num_steps_stopped = floor(switch_cost / dt);

fr = @(speed,u,q,theta) -speed.*cos(theta - (-1)^q.*u);
ftheta = @(speed,u,q,theta,r,a) speed./r.*sin(theta - (-1)^q.*u) + a;

% Uniform discretization of thresholds on [0,Initial_budget]
curr_ds = Max_cost_from_PDE_solver / num_slices;
budget_list = 0:curr_ds:Initial_budget;

mB = length(budget_list);
store_size  = 1000;

% Initializing vector to store random time-to-target values for each sample path
Xcost = zeros(1,sample_size);
R_set = cell(1,store_size);
Theta_set = cell(1,store_size);
Switch_set = cell(1,store_size);
Time_set = cell(1,store_size);
Budget_set = cell(1,store_size);
Speed_set = cell(1,store_size);
Alpha_set = cell(1,store_size);
Switch_number_set = zeros(1,sample_size);

tt = linspace(0,2*pi,theta_num);
rr = linspace(0,r_max,r_num);

%% Initial state
r0 =  rr(r0_indx);
t0 =  tt(t0_indx);
%% MC Main loop
tic
parfor ii = 1 : sample_size
    %ICs
    switch_counter = 0;
    rloc = r0;
    theta_loc = t0;
    tack = 1; %tack mode
    phi = theta_loc;
    path_r = [rloc];
    path_theta = [theta_loc];
    path_alpha = [];
    path_wind_angle = [phi];
    next_budget = [Initial_budget];
    cost_set = [0]; %accumulative cost
    speed_list = [];
    %potential storage array of indices if running out of budget
    after_indx =[];
    switch_indx = [];
    j = 1;
    switch_flag = false;

    while path_r(end) > target_radius
        idx_r = ceil( (path_r(j)-r_min)/dr ) + 1;
        idx_theta = ceil( (path_theta(j) - theta_min)/dtheta ) + 1;
        if idx_theta == 1
            idx_theta = idx_theta + 1;
        end

        %our definition of the relative angle
        path_alpha(j) = mod(path_wind_angle(j) - path_theta(j), 2*pi);

        if ((next_budget(j) <= 0))
            %decide if switching occurs
            switch_neighbors = [Switchgrid_stat(idx_theta,idx_r, tack),...
                Switchgrid_stat(idx_theta-1,idx_r-1, tack),...
                Switchgrid_stat(idx_theta-1,idx_r, tack),...
                Switchgrid_stat(idx_theta,idx_r-1, tack)];

            %find the nearest indices to the current state (risk-neutral)
            [nearest_r, nearest_theta] = ...
                nearest_neighbor_stat(idx_r,idx_theta,dr,dtheta,path_r(j),path_theta(j));

            if nearest_r == target_indx
                nearest_r = nearest_r +1;
            end
            %Read the switchgrid (risk-neutral)
            Switch_nearest = Switchgrid_stat(nearest_theta,nearest_r,tack);
            %Determine whether to switch or not (risk-neutral)
            switch_flag = switch_or_not_stat(switch_neighbors,Switch_nearest,choice);
            
            if switch_flag
                switch_counter = switch_counter + 1;
                tack = 3 -tack;
                switch_indx = [switch_indx,j];
                % Evolve the wind further for this duration
                [j,path_theta,path_r,...
                    path_wind_angle, path_alpha,cost_set,next_budget,speed_list] ...
                    = update_in_switch(j,path_theta,path_r,path_wind_angle,path_alpha,...
                    cost_set,next_budget,speed_list,dt,root_dt,...
                    num_steps_stopped, drift, sig, Initial_budget);
                switch_flag = false;
            else
                % find the optimal ctrl
                if nearest_r == target_indx
                    nearest_r = nearest_r + 1;
                end
                u = Policy_mat_stat(nearest_theta,nearest_r,tack);
                % Using speed data from polar plots
                jj = find(controls==u);
                speed = speeds(jj);
                speed_list(j) = speed;

                % Euler-Maruyama to evolve the system state
                j = j+1;
                dW = root_dt * normrnd(0,1); %Brownian increment
                path_r(j) = path_r(j-1) + fr(speed,u,tack,path_theta(j-1))*dt;
                path_theta(j) = mod(path_theta(j-1) ...
                    + ftheta(speed,u,tack,path_theta(j-1),path_r(j-1),drift)*dt ...
                    + (sig*dW), 2*pi);
                path_wind_angle(j) = mod(path_wind_angle(j-1) + (drift*dt) + (sig*dW),2*pi);
                path_alpha(j) = mod(path_wind_angle(j) - path_theta(j), 2*pi);

                if path_r(end) > r_max
                    cost_set(j) = 1e6;
                    break
                elseif path_r(end) <= target_radius
                    break
                else
                    cost_set(j) = cost_set(j-1) + dt;
                end

            end
        else %within the budget
            idx_s = find(next_budget(j) <= budget_list,1);

            switch_neighbors_top = double([mSwitchgrid_thres.Data.ss(idx_theta,idx_r,tack,idx_s),...
                mSwitchgrid_thres.Data.ss(idx_theta-1,idx_r-1,tack,idx_s),...
                mSwitchgrid_thres.Data.ss(idx_theta-1,idx_r,tack,idx_s),...
                mSwitchgrid_thres.Data.ss(idx_theta,idx_r-1,tack,idx_s)]);

            switch_neighbors_bot = double([mSwitchgrid_thres.Data.ss(idx_theta,idx_r,tack,idx_s-1),...
                mSwitchgrid_thres.Data.ss(idx_theta-1,idx_r-1,tack,idx_s-1),...
                mSwitchgrid_thres.Data.ss(idx_theta-1,idx_r,tack,idx_s-1),...
                mSwitchgrid_thres.Data.ss(idx_theta,idx_r-1,tack,idx_s-1)]);


            switch_neighbors_top(switch_neighbors_top == 3) = 0;
            switch_neighbors_top(switch_neighbors_top == 2) = 1;

            switch_neighbors_bot(switch_neighbors_bot == 3) = 0;
            switch_neighbors_bot(switch_neighbors_bot == 2) = 1;
            %construct a switchgrid data-cube surrounding the current state 
            switch_neighbors = [switch_neighbors_top, switch_neighbors_bot];
            %find the nearest indices to the current state (risk-aware)
            [nearest_r, nearest_theta, nearest_s] = ...
                nearest_neighbor_thres(idx_r,idx_theta,idx_s, dr, dtheta,curr_ds,path_r(j),path_theta(j),next_budget);

            if nearest_r == target_indx
                nearest_r = nearest_r+1;
            end

            if nearest_s == 1
                nearest_s = nearest_s+1;
            end
            %Read the switchgrid (risk-aware)
            Switch_nearest = double(mSwitchgrid_thres.Data.ss(nearest_theta,nearest_r,tack,nearest_s));

            if Switch_nearest == 3
                Switch_nearest = 0;
            elseif Switch_nearest == 2
                Switch_nearest = 1;
            end
             %Determine whether to switch or not (risk-aware)
            switch_flag = switch_or_not_thres(switch_neighbors,Switch_nearest,choice);

            if switch_flag
                switch_counter = switch_counter + 1;
                tack = 3 -tack;
                switch_indx = [switch_indx,j];
                % Evolve the wind further for this duration
                [j,path_theta,path_r,...
                    path_wind_angle, path_alpha,cost_set,next_budget,speed_list] ...
                    = update_in_switch(j,path_theta,path_r,path_wind_angle,path_alpha,...
                    cost_set,next_budget,speed_list,dt,root_dt,...
                    num_steps_stopped, drift, sig, Initial_budget);
                switch_flag = false;
            else
                % find the optimal ctrl
                if nearest_r == target_indx
                    nearest_r = nearest_r+1;
                end
                u = double(mPolicy_mat_thres.Data.uu(nearest_theta,nearest_r,tack,nearest_s))/1e4;
                u_right_index = find(u <= base_controls,1);
                u_arr = [base_controls(u_right_index -1), base_controls(u_right_index)];
                speed_arr = [base_speeds(u_right_index -1), base_speeds(u_right_index)];
                speed = LinearInterp_speed(speed_arr,u_arr,u);
                speed_list(j) = speed;

                % Euler-Maruyama to evolve the system state
                j = j+1;
                dW = root_dt * normrnd(0,1); %Brownian increment
                path_r(j) = path_r(j-1) + fr(speed,u,tack,path_theta(j-1))*dt;
                path_theta(j) = mod(path_theta(j-1) ...
                    + ftheta(speed,u,tack,path_theta(j-1),path_r(j-1),drift)*dt ...
                    + (sig*dW), 2*pi);
                path_wind_angle(j) = mod(path_wind_angle(j-1) + (drift*dt) + (sig*dW),2*pi);
                path_alpha(j) = mod(path_wind_angle(j) - path_theta(j), 2*pi);
                
                if path_r(end) > r_max
                    cost_set(j) = 1e6;
                    break
                elseif path_r(end) <= target_radius
                    break
                else
                    cost_set(j) = cost_set(j-1) + dt;
                end

            end

        end
        next_budget(j) = Initial_budget - cost_set(j);

    end
    Xcost(ii) = cost_set(end);
    Switch_number_set(ii) = switch_counter;
    if ii <= store_size
        R_set{ii} = path_r;
        Theta_set{ii} = path_theta;
        Switch_set{ii} = switch_indx;
        Time_set{ii} = cost_set;
        Budget_set{ii} = next_budget;
        Speed_set{ii} = speed_list;
        Wind_set{ii} = path_wind_angle;
        Alpha_set{ii} = path_alpha;
    end
    ii

end
Total_runtime = toc


%% ECDF plot
% memory map of the value function
v_1 = memmapfile(filename_valuefn,'Format',{'int16',[theta_num, r_num,2, num_slices+1],'vv'});
%ECDF
[f_1,x_1] = ecdf(Xcost);
s_indx = find(Xcost <= Initial_budget,1);
val_PDE = double(v_1.Data.vv(t0_indx,r0_indx,1,length(budget_list)))/1e4
val_CDF = f_1(s_indx)
abs_diff = abs(val_CDF - val_PDE)
sample_mean = mean(Xcost)
figure
plot(x_1,f_1,'linewidth',1.75);
hold on
yline(val_PDE,'r:','linewidth',1.75);
xline(Initial_budget,'k:','linewidth',1.75);
hold off
xlabel('Threshold (s)')
ylabel('Probability of success')
str1 = sprintf(' = (%.2f, %.2f, %g)', r0, t0, Initial_budget);
str2 = 'ECDF starting from $(\hat{r}, \hat{\theta}, \hat{s})$';
title([str2,str1],'Interpreter','latex');
legend('ECDF','Value function','location','southeast','fontsize',14);
axis tight
grid on

%% A sample path visualization
indxx = 1;
Sample_cost = Xcost(indxx);
Sample_r = R_set{indxx};
Sample_theta = Theta_set{indxx};
Sample_switch = Switch_set{indxx};
Sample_alpha = Alpha_set{indxx};
y0 = Sample_r(1).*cos(Sample_alpha(1));
Sample_x = Sample_r.*sin(Sample_alpha);
Sample_y = y0 - Sample_r.*cos(Sample_alpha);


figure
hold on
DrawTarget(target_radius,y0);
plot(Sample_x,Sample_y,'-','linewidth',1.1);
plot(Sample_x(Sample_switch),Sample_y(Sample_switch),'r.','markersize',15);
hold off
axis equal
grid on
title(sprintf('Total time to target: %.2f',Sample_cost));
end
%% subfunctions
function Indx = nearest_index(right_indx, xloc, dx)
if ((right_indx - 1)*dx - xloc) <= 0.5*dx
    Indx = right_indx;
else
    Indx = right_indx - 1;
end
end

function [nearest_r, nearest_theta] = nearest_neighbor_stat(idx_r,idx_theta,dr, dtheta,r,theta)
nearest_r = nearest_index(idx_r, r, dr);
nearest_theta = nearest_index(idx_theta, theta,dtheta);
end

function [nearest_r, nearest_theta, nearest_s] ...
    = nearest_neighbor_thres(idx_r,idx_theta,idx_s, dr, dtheta, ds,r,theta,s)
nearest_r = nearest_index(idx_r, r, dr);
nearest_theta = nearest_index(idx_theta, theta,dtheta);
nearest_s = nearest_index(idx_s, s, ds);
end

function switch_flag = switch_or_not_stat(neighbors,Switch_nearest, choice)
switch choice
    case 'conservative'
        if sum(neighbors) == 4
            switch_flag = true;
        else
            switch_flag = false;
        end
    case 'aggressive'
        if sum(neighbors) > 0
            switch_flag = true;
        else
            switch_flag = false;
        end
    case 'majority'
        if sum(neighbors) >= 2
            switch_flag = true;
        else
            switch_flag = false;
        end
    case 'nearest'
        if Switch_nearest == 1
            switch_flag = true;
        else
            switch_flag = false;
        end
end

end


function switch_flag = switch_or_not_thres(neighbors,Switch_nearest,choice)
switch choice
    case 'conservative'
        if sum(neighbors) == 8
            switch_flag = true;
        else
            switch_flag = false;
        end
    case 'aggressive'
        if sum(neighbors) > 0
            switch_flag = true;
        else
            switch_flag = false;
        end
    case 'majority'
        if sum(neighbors) >= 5
            switch_flag = true;
        else
            switch_flag = false;
        end
    case 'nearest'
        if Switch_nearest == 1
            switch_flag = true;
        else
            switch_flag = false;
        end
end

end

function [j,path_theta,path_r,...
    path_wind_angle, path_alpha,cost_set,next_budget,speed_list] ...
    = update_in_switch(j,path_theta,path_r,path_wind_angle,path_alpha,...
    cost_set,next_budget,speed_list,dt,root_dt,...
    num_steps_stopped, drift, sig, Initial_budget)

for idx = 1:num_steps_stopped
    j = j+1;
    dW = root_dt * normrnd(0,1); % Euler-Maruyama
    path_theta(j) = mod(path_theta(j-1) + (drift*dt) + (sig*dW),2*pi);
    path_wind_angle(j) = mod(path_wind_angle(j-1) + (drift*dt) + (sig*dW),2*pi);
    path_r(j) = path_r(j-1);
    path_alpha(j) = path_alpha(j-1);
    cost_set(j) = cost_set(j-1) + dt;
    speed_list(j) = 0;
    next_budget(j) = Initial_budget-cost_set(j);
end
end

function val=LinearInterp_speed(val_arr,angle_arr,u_loc)
dx = angle_arr(2) - angle_arr(1);
val= val_arr(1) + ((val_arr(2)-val_arr(1))/dx)*(u_loc - angle_arr(1));
end

function DrawTarget(target_radius,target_dist)
rectangle('Position',[-target_radius, target_dist-target_radius, 2*target_radius, 2*target_radius],...
    'Curvature',[1 1],'Facecolor',[1 0 0.563]);
hold on;
end
