% Example 2

% This script employs the TZPC algorithm to compute optimal control inputs
% for Example 2 provided in the referenced paper. 

% Prerequisites: 
%
% - Ensure the MPT and Yalmip toolbox are installed, along with the MOSEK solver.
% - Add 'RequiredFiles_TZPC' folder to the MATLAB path.
% - Add 'Saved Workspace' folder to the MATLAB path.

% Date of creation: 2024-01-29
% Date of update:   2024-08-17

% Set random seed for reproducibility
rand('seed',4500);

% Clear command window and variables
clc; clear variables;close all;

addpath("Saved Workspace\")
addpath(genpath("..\RequiredFiles_TZPC\"))

%% Number of Samples
initpoints = 20;                 % number of trajectories
steps = 5;                       % number of steps for each trajectory
totalsamples = initpoints*steps; % total number of samples

%% Load Weather Data
load('24Jan2010.mat');
weather = table2array(T_weath(:,2));

%% System Information 
load('building_sys');
sys.A  = Ad;
sys.B  = Bd(:,1);
sys.x0 = [17;15];          % initial state
[nx,nu] = size(sys.B);     % dimension of the system state and the control input

load('Localcontrol_info'); % load local control information (K,P)
sys.K = K;


%% Process Noise Zonotope (B*Tout)
WCen    = 0;                                           % center of process noise zonotope
WGen    = 2;                                           % generator of process noise zonotope
W       = Bd(:,2) * zonotopempt([WCen,WGen]);          % process noise zonotope which represents the scaled outside temperature

%% Generate Matrix Zonotope (M_w)

% Generators of noise matrix zonotope
for i = 1 : size(W.generators,2)
    vec = W.Z(:,i+1);
    GW{1+(i-1)*totalsamples} = [ vec, zeros(nx,totalsamples-1) ];
    GW{totalsamples+(i-1)*totalsamples} = [zeros(nx,totalsamples-1)  vec];
    for j = 2 : totalsamples-1
        GW{j+(i-1)*totalsamples} = [zeros(nx,j-1) vec zeros(nx,totalsamples-j)];
    end
end

% Center of noise matrix zonotope
CW = WCen * ones(nx,totalsamples);

% Noise matrix zonotpe (M_w)
M_w = matZonotope(CW,GW);

%% State and Input Setpoint

Tin_setpoint = 22;   % indoor air setpoint temperature

% Compute setpoint for wall temperature based on indoor air setpoint
Tw_star = sys.A(2,1) * Tin_setpoint / (1 - sys.A(2,2));

% Compute optimal control input based on setpoints
U_star = (Tin_setpoint - sys.A(1,1) * Tin_setpoint - sys.A(1,2) * Tw_star) / sys.B(1,:);

% Define state and control input setpoints 
rx_real = [Tin_setpoint;Tw_star];
ru_real = U_star;

%% State Constraints
XCen_org = [22; 20.25];                              % center of state constraint zonotope
XGen = diag([2 1.25]);                               % generator of state constraint zonotope
X_org = intervalmpt(zonotopempt([XCen_org,XGen]));   % state constraint zonotope

% State constraints for generating data (used in data generation)
Xcon_data = zonotopempt([[17 15]',diag([0.5 0.5])]);

% Setpoint normalization of the constraints 
XCen = XCen_org - rx_real;                        % adjust center based on the state setpoint
X = intervalmpt(zonotopempt([XCen,XGen]));        % normalized state constraint zonotope

%% Control Input Constraints
UCen_org = 35;                                     % center of control input zonotope
UGen = 7;                                          % generator of control input zonotope
U_org = intervalmpt(zonotopempt([UCen_org,UGen])); % control input zonotope

% Setpoint normalization of the constraints 
UCen = UCen_org - ru_real;                         % adjust center based on the control input setpoint 
U = intervalmpt(zonotopempt([UCen,UGen]));         % normalized control input constraint zonotope

%% Generate Data

x(:,1) = sys.x0;
index = 1;
for i = 1 : totalsamples
    u(i) = 0.01 * rand(1) + 25;
end

for j = 1 : nx : initpoints * nx
    x(j:j+nx-1,1) = randPoint(Xcon_data);
    for i = 1 : steps
        utraj(j,i) = u(index);
        x(j:j+nx-1,i+1) = sys.A * x(j:j+nx-1,i) + sys.B * u(index) + randPoint(W);
        index = index + 1;
    end
end

index_0 = 1;
index_1 = 1;

for j = 1 : nx : initpoints * nx
    for i = 2 : steps + 1
        x_meas_vec_1(:,index_1) = x(j:j+nx-1,i);
        index_1 = index_1 + 1;
    end
    for i = 1 : steps
        u_mean_vec_0(:,index_0) = utraj(j,i);
        x_meas_vec_0(:,index_0) = x(j:j+nx-1,i);
        index_0 = index_0 + 1;
    end
end

U_0 = u_mean_vec_0(:,1:totalsamples);
X_0 = x_meas_vec_0(:,1:totalsamples);
X_1 = x_meas_vec_1(:,1:totalsamples);

% Check rank condition assumption
Hankel = [X_0;U_0];

if rank(Hankel) == nx + nu 
    disp('The rank condition of the Hankel matrix is satisfied.');
else
    disp('The rank condition of the Hankel matrix is NOT satisfied.');
end

%% Compute System Matrices (M_\Sigma)
M_Sigma = matZonotope(X_1 + -1* M_w.center,M_w.generator) * pinv(Hankel);
% M_Sigma_cen = M_Sigma.center;
% M_Sigma_gen = M_Sigma.generator;

%% Compute the Lipschitz Constant Matrix and the zonotope Z_L

intM_Sigma = intervalMatrix(M_Sigma);

% Compute the Lipschitz constants using the Frobenius norm
ii = intM_Sigma.Sup;
Lip1 = norm(ii(1,:),'fro');
Lip2 = norm(ii(2,:),'fro');
L = [Lip1;Lip2];
% L = norm(intM_Sigma,'fro');
% L = L*eye(nx);

% Find delta
gamma = 1000000;
for i = 1 : totalsamples
    z1 = [X_0(:,i) ; U_0(:,i)];
    for j = 1 : totalsamples
        z2 = [X_0(:,j) ; U_0(:,j)];  
        if z1~=z2
            newgamma = norm(z1-z2);
            if newgamma < gamma
                gamma = newgamma;
            end
        end
    end
    g(:,i) = gamma;
end
% Compute the maximum gamma
gamma = max(g);

% Define zonotope Z_L
Z_L   = zonotopempt([zeros(nx,1), L * gamma / 2]);

%% Define Nominal System

A_nominal = M_Sigma.center(:, 1 : nx);
B_nominal = M_Sigma.center(:, nx+1 : end);

% Define combined matrix M_bar
M_bar = [A_nominal B_nominal]; 

% compute the controllability matrix
Co = ctrb(A_nominal,B_nominal);

% check and print the controllability condition
if length(A_nominal) == rank(Co)
    disp('The nominal system is controllable.')
else
    disp('The nominal system is NOT controllable.')
end

%% Compute Model Mismatch Zonotope (Z_M)

mismatch = X_1 - M_bar * [X_0; U_0];
min_m = min(mismatch, [], 2);
max_m  = max(mismatch, [], 2);

% Create the model mismatch zonotope Z_M
Z_M = zonotopempt(intervalmpt(min_m,max_m)) + -1 * W;

%% Compute the Closed-Loop System Matrix (A_K)
AK = A_nominal + B_nominal * sys.K;

%% Load the previously saved invariant set data
load('S_info.mat');
S = intervalmpt(S_inf,S_sup);

%% Load Alpha value
load('alpha.mat');

%% Load Tightened Constraints
load('tight_Uconst.mat');
load('tight_Xconst.mat');
tight_U_cons = intervalmpt(tight_Ucons_inf,tight_Ucons_sup);
tight_X_cons = intervalmpt(tight_Xcons_inf,tight_Xcons_sup);

%% Solve the optimization problem
N = 11;                     % prediction horizon
Simulation_time = 188;      % simulation time


x_t(:,1)   = [-2;-1];       % initial state
w_point    = randPoint(W);  % initial process noise
Cost = 0;               
Constraints = [];

for timesteps = 1 : Simulation_time
    fprintf('timestep: %d\n',timesteps)
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
    options = sdpsettings('verbose', 0,'solver', 'mosek');
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
    realu(:,timesteps)   = uPred(:,timesteps) + sys.K * ( x_t(:,timesteps) - xPred(:,timesteps)); 
    x_t(:,timesteps+1)   = sys.A * x_t(:,timesteps) + sys.B * realu(:,timesteps) + Bd(:,2) * weather(timesteps);

    if strfind(Status_ZPC{timesteps},'Successfully solved (MOSEK')==0
        disp('THE PROBLEM IS INFEASIBLE')
    end
end
% store the average execution time
executiontime = mean(T);
disp(['The execution time is: ', num2str(executiontime), ' seconds']);

%% Reachable Sets 
R{1} = zonotopempt([x_t(:,1)]);
for i = 1 : Simulation_time
    cntr  = [x_t(:,i); uPred(:,i) + sys.K * (x_t(:,i)-xPred(:,i))];
    gener = zeros(3,1);
    R{i+1} = M_Sigma * zonotopempt([cntr,gener]) + W;
    R_int{i+1} = intM_Sigma * zonotopempt([cntr,gener]) + W;

end

% Save the workspace for plotting purposes
clearvars -except Simulation_time x_t  R  R_int realu N executiontime
save(['188_Tobeplotted-TZPC-Building--N' num2str(N) '.mat']);

