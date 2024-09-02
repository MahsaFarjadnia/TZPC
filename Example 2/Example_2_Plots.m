% Example 2 plots

% This script generates plots of the final results for Example 2 as presented in the referenced paper. 

% Prerequisites: 
%
% - Add the 'Saved Workspace' folder to the MATLAB path.
% - Run 'Example_2_TZPC.m' before executing this script.
% - Add 'RequiredFiles_TZPC' folder to the MATLAB path.

% Date of creation: 2024-01-29
% Date of update:   2024-08-17

clear all;
clc;
close all;

% Parameters
savefig = 0;
N = 11;
col = lines(2);
Fontsizeax = 22;

% Load necessary data
load(['188_Tobeplotted-TZPC-Building--N' num2str(N) '.mat']);
load('X_const');
load('U_const');
load('sys.mat');
load('ref');
load('24Jan2010.mat');

% Set weather data
weather = table2array(T_weath(:,2));

% Set constraints
X_const = intervalmpt(zonotopempt([XCen,XGen]));
X_const_inf = X_const.inf;
X_const_sup = X_const.sup;

U_const = intervalmpt(zonotopempt([UCen,UGen]));
U_const_inf = U_const.inf;
U_const_sup = U_const.sup;

%% Plot Indoor Temperature Trajectories 
tc = 60;
step_new = (0:tc:(12*60))/4;
fig1 = figure('Renderer', 'painters', 'Position', [10 50 700 550]);
plot(1: length(x_t),x_t(1,:)+rx_real(1,:),'*-','color', col(1,: ),'LineWidth',0.6)
hold on
yline(rx_real(1,:),'k','LineWidth',1)
yline(X_const_inf(1)+rx_real(1,:),'k--','LineWidth',1)
yline(X_const_sup(1)+rx_real(1,:),'k--','LineWidth',1)
ax = gca;
ax.XTick = step_new;
ax.XTickLabel = strcat(num2str((0:11)'),'h');
ax.XLim(2) = 11*60/4;
xtickangle(0)
xlabel('time [min]','Interpreter','latex')
ylabel('$T_\mathrm{in}$ $[^{\circ}C]$','Interpreter','latex')
legend1 = legend('TZPC: $x_1$','setpoint','constraints','Interpreter','latex')
set(legend1,'Position',[0.648634302048456 0.771709169935527 0.285890706380208 0.181890907287598]);
ylim([19.9 24.1])

tc = 60;
step_new2 = (4*60:tc:(9*60))/4;

ax2 = axes('Position',[0.5 0.22 0.35 0.3]);
hold on 
yyaxis left
plot(1: length(x_t),x_t(1,:)+rx_real(1,:),'*-','color', col(1,: ),'LineWidth',0.6)
hold on 
yline(22,'LineWidth',1)
ylabel('$T_\mathrm{in}$','Interpreter','latex')
yyaxis right
plot(1:length(weather),weather(1:end),'*-','LineWidth',0.6)
ylabel('$T_\mathrm{out}$','Interpreter','latex')

ylim([-2,0])
ax2.XTick = step_new2;
ax2.XTickLabel = strcat(num2str((4:9)'),'h');
ax2.XLim = [60*4,9*60]/4;
xtickangle(0)


ax.FontSize = Fontsizeax;
ax.XLabel.FontSize= 27;
ax.YLabel.FontSize= 27;
outerpos = ax.OuterPosition;
tih = ax.TightInset;
left = outerpos(1) + tih(1)+0.03;
bottom = outerpos(2) + tih(2);
ax_width = outerpos(3) - tih(1) - tih(3)-0.05;
ax_height = outerpos(4) - tih(2) - tih(4);
ax.Position = [left bottom ax_width ax_height];


%% Plot Control Input 
tc = 60;
step_new = (0:tc:(12*60))/4;
fig2 = figure('Renderer', 'painters', 'Position', [10 50 700 550]);
plot(1: length(realu),realu(1,:)+ru_real,'*-','color', col(1,: ),'LineWidth',0.6)
hold on
yline(ru_real,'k','LineWidth',1)
yline(U_const_inf(1)+ru_real,'k--','LineWidth',1)
yline(U_const_sup(1)+ru_real,'k--','LineWidth',1)
ax = gca;
ax.XTick = step_new;
ax.XTickLabel = strcat(num2str((0:11)'),'h');
ax.XLim(2) = 11*60/4;
xtickangle(0)
xlabel('time [min]','Interpreter','latex')
ylabel('$T_\mathrm{mr}$ $[^{\circ}C]$','Interpreter','latex')
legend1 = legend('TZPC: $u$','setpoint','constraints','Interpreter','latex');
% set(legend1,'Position',[0.700560389927455 0.803127273212779 0.243842439197359 0.158327272588556]);
set(legend1,'Position',[0.674680604480562 0.763466667984471 0.270840105329241 0.181890907287598]);
ylim([27.9 42.6])
ax.FontSize = Fontsizeax;
ax.XLabel.FontSize= 27;
ax.YLabel.FontSize= 27;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4)-0.01;
ax.Position = [left bottom ax_width ax_height];

% Save the figures
if savefig 
    saveas(fig1, ['Tin-Tout-TZPC-building-N' num2str(N) '.fig']);
    saveas(fig1, ['Tin-Tout-TZPC-building-N' num2str(N) '.png']);
    saveas(fig1, ['Tin-Tout-TZPC-building-N' num2str(N) '.eps'],'epsc');

    saveas(fig2, ['U-hourly-TZPC-building-N' num2str(N) '.fig']);
    saveas(fig2, ['U-hourly-TZPC-building-N' num2str(N) '.png']);
    saveas(fig2, ['U-hourly-TZPC-building-N' num2str(N) '.eps'],'epsc');
end
