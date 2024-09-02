% Example 1

% This script employs TZPC algorithm to compute optimal control inputs
% for the Example 1 provided in the referenced paper. 

% Prerequisites: 
%
% - Ensure the MPT toolbox and YALMIP solvers are installed and accessible.
% - Add 'RequiredFiles_TZPC' folder to the MATLAB path.
% - Add 'Saved Workspace' folder to the MATLAB path.

% Date of creation: 2024-01-29
% Date of update:   2024-08-14


% Clear command window and variables
clc;
clear all;
close all;

addpath("Saved Workspace\")
addpath(genpath("..\RequiredFiles_TZPC\"))

% Set random seed for reproducibility
rand('seed',4500);

%% Load System Information 
load('sys.mat');
%% Load Process Noise Zonotope (W)
load('W_const.mat');
W  = zonotopempt([WCen*ones(nx,1),WGen]); % process noise zonotope
%% Generate Noise Matrix Zonotope (M_w)
load('M_w.mat');
% noise matrix zonotpe (M_w)
M_w = matZonotope(CW,GW);

%% Load Control and State Constraints
load('U_const');
load('X_const');
U = intervalmpt(zonotopempt([UCen,UGen])); % control input zonotope
X = intervalmpt(zonotopempt([XCen,XGen]));  % state constraint zonotope

%% Load Data
load('data.mat')
% Check rank condition assumption
Hankel = [X_0;U_0];
if rank(Hankel) == nx + nu 
    disp('The rank condition of the Hankel matrix is satisfied.')
else
    disp('The rank condition of the Hankel matrix is NOT satisfied.')
end

%% Compute System Matrices (M_\Sigma)
M_Sigma = matZonotope(X_1 + -1* M_w.center,M_w.generator) * pinv(Hankel);

%% Define Nominal System

load('intM_Sigma_intervalMatrix.mat');
% extract the infimum and supremum of the interval matrix
intM_Sigma_inf = intM_Sigma.Inf;
intM_Sigma_sup = intM_Sigma.Sup;

A_nominal = intM_Sigma_inf(:, 1 : nx);
B_nominal = intM_Sigma_inf(:, nx+1 : end);

% compute the controllability matrix
Co = ctrb(A_nominal,B_nominal);

% check and print the controllability condition
if length(A_nominal) == rank(Co)
    disp('The nominal system is controllable.')
else
    disp('The nominal system is NOT controllable.')
end

% define Combined Matrix M_bar
M_bar = [A_nominal B_nominal]; 

%% Compute the Closed-Loop System Matrix (A_K)
AK = A_nominal + B_nominal * sys.K;

%% Load the previously saved invariant set data
load('S_info.mat');
S = intervalmpt(S_inf,S_sup);
%% Load Alpha Value
load('alpha.mat');
% fprintf('alpha = %.3f\n',alpha)

%% Load Tightened Constraints
load('tight_Uconst.mat');
load('tight_Xconst.mat');
tight_U_cons = intervalmpt(tight_Ucons_inf,tight_Ucons_sup);
tight_X_cons = intervalmpt(tight_Xcons_inf,tight_Xcons_sup);

%% Solve the Optimization Problem
N = 7;                     % prediction horizon
Simulation_time = 15;      % simulation time
run_num = 5;               % number of runs to compute the average execution time

for kk = 1 : run_num
    x_t(:,1)   = sys.x0;       % initial state
    w_point    = randPoint(W); % initial process noise
    Cost = 0;               
    Constraints = [];

    for timesteps = 1 : Simulation_time
        % Define decision variables
        u = sdpvar(nu*ones(1,N),ones(1,N));
        x = sdpvar(nx*ones(1,N+1),ones(1,N+1));
        % Define constraints
        for q = 1 : N
            Constraints = [Constraints, x{:,q + 1} == A_nominal * x{:,q} + B_nominal * u{:,q},...
                u{:,q} <= tight_U_cons.sup,...
                u{:,q} >= tight_U_cons.inf,...
                 x{:,q} <= tight_X_cons.sup,...
                 x{:,q} >= tight_X_cons.inf,...
                x_t(:,timesteps) <= intervalmpt(S).sup + x{:,1},...
                x_t(:,timesteps) >= intervalmpt(S).inf + x{:,1}];
                % Accumulate cost
                Cost = Cost + norm(x{:,q},2) + 10^(-2) * abs(u{:,q});
        end

        % Final cost and constaint
        Cost = Cost + x{:,N+1}' * P * x{:,N+1};
        Constraints = [Constraints, x{:,N+1}' * P * x{:,N+1} <= alpha];
        options = sdpsettings('verbose', 0);
        % Solve the optimization problem
        tic;
        Problem   = optimize(Constraints,Cost,options)
        T(timesteps) = toc;

        % Extract the results
        Objective             = double(Cost);
        Status_ZPC{timesteps} = Problem.info;
        uPred(:,timesteps)    = double(u{1});
        xPred(:,timesteps)    = double(x{1});

        % Calculate the real control input and apply it to the system
        realu(:,timesteps)    = uPred(:,timesteps) + sys.K * ( x_t(:,timesteps) - xPred(:,timesteps));
        w_point = randPoint(W);
        x_t(:,timesteps+1)   = sys.A * x_t(:,timesteps) + sys.B * realu(:,timesteps) + w_point;
    end
    % Store the average execution time for this run
    TZPCtime(kk) = mean(T);
    disp(['Run ', num2str(kk), ': Average execution time = ', num2str(TZPCtime(kk)), ' seconds']);
end
% Calculate and display the overall average execution time
executiontime = mean(TZPCtime);
disp(['Overall average execution time: ', num2str(executiontime), ' seconds']);

%% Reachable sets 
R{1} = zonotopempt([x_t(:,1)]);
for i = 1 : Simulation_time
    cntr  = [x_t(:,i); uPred(:,i) + sys.K * (x_t(:,i)-xPred(:,i))];
    gener = zeros(3,1);
    R{i+1} = M_Sigma * zonotopempt([cntr,gener]) + W;
end

% Save the workspace for plotting purposes
clearvars -except Simulation_time x_t  R realu N executiontime xPred
save(['Tobeplotted-TZPC-N' num2str(N) '.mat']);
