% Example 2

% This script implements a simplified thermal model for a single zone 
% as described in the referenced paper. The script sets up both continuous-time 
% and discrete-time system matrices.

% Date of creation: 2024-01-29
% Date of update:   2024-08-15

clc;
close all;
clear variables;

%% System Parameters
% Model Parameters
n = 2; % Number of states (indoor temperature, wall temperature)
m = 2; % Number of inputs and disturbances (mean radiant tempearute of radiators, outside temperature)

% The original model is x(t+1) = A_d x(t) + B_d [u(t);w(t)]. 
% However, for simplicity in the paper B_d(2,:)w(t) is considered as w(t).

% Thermal Parameters

Ci     = 34600;     % [J/K]       % Thermal capacitance of indoor air
Cw     = 560000;    % [J/K]       % Thermal capacitance of walls
Riw    = 0.01;      % [K/W]       % Thermal resistance between indoor air and walls
Riwin  = 0.162;     % [K/W]       % Thermal resistance between indoor air and windows
Rwe    = 0.34;      % [K/W]       % Thermal resistance between walls and outside air
CCOtwo = 420;       % [ppm]       % CO2 concentration of the air inlet
gCOtwo = 0.0000103; % [kgCO2/s]   % CO2 emitted per occupant
C      = 125;       % [W]         % Heat emission per occupant
G      = 0.61;      % [-]         % Solar heat gain coefficient
Awin   = 4.22;      % [m^2]       % Total window area
Cpa    = 1006;      % [J/kgK]     % Air specific heat capacity 
Arad   = 3;         % [m^2]       % Surface area of radiator 
hrad   = 10;        % [W/m^2K]    % Heat transfer coefficient of the radiator 

%% Continuous-Time System Matrices
C_d = diag([Ci Cw]);

a11 = -1*(1/Riw + 1/Riwin) - Arad * hrad;
a12 = 1/Riw;
a21 = a12;
a22 = -1*(1/Riw + 1/Rwe);

A = C_d\[a11 a12;a21 a22];

Bu = C_d\[Arad * hrad;0];

Bw = C_d\[1/Riwin G*Awin C; 1/Rwe   0     0];

Bw = Bw(:,1); % only the first column is used because we do not consider the impact of solar radiation
B = [Bu,Bw];

%% Continuous-Time State-Space System
C = eye(n); D = zeros(n,m);
G = ss(A,B,C,D);

%% Discrete-Time System Matrices
delta_t = 240; % discretization time (4 min)
Ad  = A*delta_t + eye(n);
Bd  = B*delta_t;
% Create discrete-time state-space system (optional)
% G_discrete = ss(Ad, Bd, C, D, delta_t);

%% Display the Discrete-Time System Matrices Used in the Paper

fprintf(['Discrete-time system matrices used in the paper:','\nA = \n']);
disp(Ad);
fprintf(['B = \n']);
disp(Bd(:, 1));

