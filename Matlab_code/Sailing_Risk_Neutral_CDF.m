function Xcost = Sailing_Risk_Neutral_CDF(DataFile_dir,sample_size,time_step,choice,r0_indx,t0_indx)
%This function computes the CDF y = Pr(T <= s) that measures the
%probability of success where the cumulative time-to-target T is 
%within any positive threshold value "s" using the Risk-Neutral policy
%computed from Miles and Vladimirsky
%https://ieeexplore.ieee.org/document/9655258
%
%Author: MingYi Wang, Cornell University
%Last Modified: 01/2024
%
%DataFile_dir (input): The directory under the folder "output" where the
%                       user stores their data files (risk-neutral)
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
%
%Xcost (output): Array of the accumulative (random) time-to-target of each sample
%
%% Read in parameter data
% Define file names
filename_parameters = ['../Risk_Neutral/output/',DataFile_dir,'/PathTracer_DomainParameters'];
filename_controls = ['../Risk_Neutral/output/',DataFile_dir,'/ControlList.dat'];
filename_speeds = ['../Risk_Neutral/output/',DataFile_dir,'/SpeedList.dat'];

filename_policy_stat = ['../Risk_Neutral/output/',DataFile_dir,'/PathTracer_Policies.dat'];
filename_switch_stat = ['../Risk_Neutral/output/',DataFile_dir,'/PathTracer_Switchgrid.dat'];
filename_value_stat = ['../Risk_Neutral/output/',DataFile_dir,'/PathTracer_ValueFunc.dat'];

% Load in Parameters
precision = 'double';
uFile_params = fopen(filename_parameters);
u_params = fread(uFile_params, 14, precision);

u_params = u_params';

theta_num = u_params(1)+1 % = 100 * gfactor = number of theta points
theta_max = u_params(2) % maximum theta in domain = 2pi
r_num = u_params(3)+1 % = 100* g_factor = number of radius points
r_max = u_params(4) % maximum radius in domain
% g_factor = u_params(5) 
target_radius = u_params(6) % radius of target region 
r_min = u_params(7)
theta_min = u_params(8)
dr = u_params(9)
dtheta = u_params(10)

sig = u_params(11)
drift = u_params(13)

switch_cost = u_params(12)
num_ctrl = u_params(14)
fclose(uFile_params);
target_indx = floor(target_radius/dr)+1;

% Getting Controls and Corresponding Speeds
uFile_cntrls = fopen(filename_controls);
controls = fread(uFile_cntrls, num_ctrl, precision);
uFile_speeds = fopen(filename_speeds);
speeds = fread(uFile_speeds, num_ctrl, precision);

%% Read in Optimal Policy Data and Switchgrids
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

% ------ Read in value function, store in a matrix
vFile = fopen(filename_value_stat);
Value_stat = fread(vFile, 2*theta_num*r_num, precision);
fclose(vFile);
% Reshaping
Value_stat = reshape(Value_stat,[theta_num, r_num, 2]);
for mode = 1:2
    temp = Value_stat(:,:,mode);
    temp(end,:) = temp(1,:);
    Value_stat(:,:,mode) = temp;
    
end

Policy_mat_stat(Policy_mat_stat == -1000) = NaN;
vv = Value_stat(:,:,1);
%% Monte Carlo Simulation Initialization
dt = time_step;
root_dt = sqrt(dt);

num_steps_stopped = floor(switch_cost / dt);

fr = @(speed,u,q,theta) -speed.*cos(theta - (-1)^q.*u);
ftheta = @(speed,u,q,theta,r,a) speed./r.*sin(theta - (-1)^q.*u) + a;

% Initializing vector to store random time-to-target values for each sample path
Xcost = zeros(1,sample_size);
R_set = cell(1,sample_size);
Theta_set = cell(1,sample_size);
Switch_set = cell(1,sample_size);
Wind_set = cell(1,sample_size);
Alpha_set = cell(1,sample_size);

tt = linspace(0,2*pi,theta_num);
rr = linspace(0,r_max,r_num);
%% Initial state
r0 = rr(r0_indx);
t0 = tt(t0_indx);
%% MC Main loop
tic
parfor ii = 1 : sample_size
    %ICs
    rloc = r0;
    theta_loc = t0;
    tack = 1; %tack mode
    phi = theta_loc;
    path_r = [rloc];
    path_theta = [theta_loc];
    path_alpha = [];
    path_wind_angle = [phi];

    cost_set = [0]; %accumulative cost
    switch_indx = [];

    j = 1;
    while path_r(j) > target_radius
        idx_r = ceil( (path_r(j)-r_min)/dr ) + 1;
        idx_theta = ceil( (path_theta(j) - theta_min)/dtheta ) + 1;
        
        if idx_theta == 1
            idx_theta = idx_theta + 1;
        end

        %our definition of the relative angle
        path_alpha(j) = mod(path_wind_angle(j) - path_theta(j), 2*pi);

        %decide if switching occurs
        switch_neighbors = [Switchgrid_stat(idx_theta,idx_r, tack),...
            Switchgrid_stat(idx_theta-1,idx_r-1, tack),...
            Switchgrid_stat(idx_theta-1,idx_r, tack),...
            Switchgrid_stat(idx_theta,idx_r-1, tack)];
        
        %find the nearest indices to the current state
        [nearest_r, nearest_theta] = ...
            nearest_neighbor_stat(idx_r,idx_theta,dr,dtheta,path_r(j),path_theta(j));
        
        if nearest_r == target_indx
            nearest_r = nearest_r +1;
        end
        %Read the switchgrid
        Switch_nearest = Switchgrid_stat(nearest_theta,nearest_r,tack);
        %Determine whether to switch or not
        switch_flag = switch_or_not_stat(switch_neighbors,Switch_nearest,choice);
        if switch_flag
            tack = 3 -tack;
            switch_indx = [switch_indx,j];
            % Evolve the wind further for this duration
            [j,path_theta,path_r,cost_set,path_alpha,path_wind_angle]= ...
                update_in_switch(j,path_theta,path_r,cost_set,dt, root_dt, num_steps_stopped,...
                drift, sig,path_alpha,path_wind_angle);
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
            
            %Euler-Maruyama to evolve the system state
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
            else
                cost_set(j) = cost_set(j-1) + dt;
            end
        end
        
    end
    Xcost(ii) = cost_set(end);
    if ii < 1001
        R_set{ii} = path_r;
        Theta_set{ii} = path_theta;
        Switch_set{ii} = switch_indx;
        Wind_set{ii} = path_wind_angle;
        Alpha_set{ii} = path_alpha;
    end
    ii
    
end
Total_runtime = toc
%% ECDF plot
[f_1,x_1] = ecdf(Xcost);
figure
hold on
plot(x_1,f_1,'linewidth',2);
xline(mean(Xcost),'b-','linewidth',1.5);
%compare with the value function
xline(vv(t0_indx,r0_indx),'r:','linewidth',1.5); 
hold off
lgd = legend('ECDF','Sample mean','Value function','location','southeast');
axis tight
lgd.FontSize = 14;
ICname = sprintf(' = (%.2f, %.2f)',r0,t0);
title(['ECDF starting from ($\hat{r},\hat{\theta}$)',ICname],'Interpreter','latex');

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


function [j,path_theta,path_r,cost_set,path_alpha,path_wind_angle]...
    = update_in_switch(j,path_theta,path_r,cost_set,dt, root_dt, num_steps_stopped, drift, sig,path_alpha,path_wind_angle)
% global dt root_dt num_steps_stopped drift sig
for idx = 1:num_steps_stopped
    j = j+1;
    dW = root_dt * normrnd(0,1); % Euler-Maruyama
    path_theta(j) = mod(path_theta(j-1) + (drift*dt) + (sig*dW),2*pi);
    path_wind_angle(j) = mod(path_wind_angle(j-1) + (drift*dt) + (sig*dW),2*pi);
    path_r(j) = path_r(j-1);
    path_alpha(j) = path_alpha(j-1);
    cost_set(j) = cost_set(j-1) + dt;
end
end

function DrawTarget(target_radius,target_dist)
rectangle('Position',[-target_radius, target_dist-target_radius, 2*target_radius, 2*target_radius],...
    'Curvature',[1 1],'Facecolor',[1 0 0.563]);
hold on;
end
