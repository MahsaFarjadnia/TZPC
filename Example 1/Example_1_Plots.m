% Example 1 plots

% This script generates plots of the final results for Example 1 as presented in the referenced paper. 

% Prerequisites: 
%
% - Add the 'Saved Workspace' folder to the MATLAB path.
% - Run 'Example_1_TZPC.m' before executing this script.
% - Add 'RequiredFiles_TZPC' folder to the MATLAB path.

% Date of creation: 2024-01-29
% Date of update:   2024-08-14

clear all;
clc;
close all;

% Parameters
savefig = 0;
N = 7;
endtime = 10;
Fontsizeax = 22;
col = lines(2);

% Load necessary data
load(['Tobeplotted-ZPC-N' num2str(N) '.mat']);
load(['Tobeplotted-TZPC-N' num2str(N) '.mat']);
load('intM_Sigma_intervalMatrix');
load('W_const');
load('X_const');
load('U_const');
load('sys.mat');

% Set interval constraints
X_const = intervalmpt(zonotopempt([XCen,XGen]));
X_const_inf = X_const.inf;
X_const_sup = X_const.sup;

U_const = intervalmpt(zonotopempt([UCen,UGen]));
U_const_inf = U_const.inf;
U_const_sup = U_const.sup;

%% plot u 

fig1 = figure('Renderer', 'painters', 'Position', [10 50 700 550]);
plot(1: endtime ,U_zpc(1,1:endtime),'v-','color', col(2,: ),'LineWidth',1);
hold on
plot(1: endtime,realu(1,1:endtime),'*-','color', col(1,: ),'LineWidth',1);
yline(0,'k','LineWidth',1);
yline(U_const_inf(1),'k--','LineWidth',1);
yline(U_const_sup(1),'k--','LineWidth',1);
ax = gca;
xlabel('time','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
legend1 = legend('ZPC','TZPC','setpoint','constraints','Interpreter','latex');
set(legend1,'Position',[0.696102411179315 0.662581819476504 0.263000357491629 0.240399997595585]);
ax.FontSize = Fontsizeax;
ax.XLabel.FontSize= 27;
ax.YLabel.FontSize= 27;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%% plot x1
fig2 = figure('Renderer', 'painters', 'Position', [10 50 700 550]);
plot(1: endtime,x_t_zpc(1,1:endtime),'v-','color', col(2,: ),'LineWidth',1)
hold on
plot(1: endtime,x_t(1,1:endtime),'*-','color', col(1,: ),'LineWidth',1)
yline(0,'k','LineWidth',1)
yline(X_const_inf(1),'k--','LineWidth',1)
yline(X_const_sup(1),'k--','LineWidth',1)
ax = gca;
xlabel('time','Interpreter','latex')
ylabel('$x_1$','Interpreter','latex')
legend1 = legend('ZPC: $x_1$','TZPC: $x_1$','setpoint','constraints','Interpreter','latex','Location','northeast');
set(legend1,'Position',[0.667193806966146 0.207309066310074 0.288044288271949 0.250872751871744]);
ax.FontSize = Fontsizeax;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%% plot x2
fig3 = figure('Renderer', 'painters', 'Position', [10 50 700 550]);
plot(1: endtime,x_t_zpc(2,1:endtime),'v-','color', col(2,: ),'LineWidth',1);
hold on;
plot(1: endtime,x_t(2,1:endtime),'*-','color', col(1,: ),'LineWidth',1);
yline(0,'k','LineWidth',1);
yline(X_const_inf(2),'k--','LineWidth',1);
yline(X_const_sup(2),'k--','LineWidth',1);
ax = gca;
xlabel('time','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
legend1 = legend('ZPC: $x_2$','TZPC: $x_2$','setpoint','constraints','Interpreter','latex');
set(legend1,'Position',[0.680737697056361 0.189127229797892 0.271643255324591 0.250872751871745]);
ax.FontSize = Fontsizeax;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
ylim([-2.1 2.1]);
%}
%% plot reachable sets
fig4 = figure('Renderer', 'painters', 'Position', [10 50 700 550]);
clear pLeg;
pLeg(1) = plot(0,0,'ok','LineWidth',2,'MarkerFaceColor',[1 0.5 0]);
hold on;
pLeg(3) = plot(x_t_zpc(1,1:endtime),x_t_zpc(2,1:endtime),'v-','color', col(2,: ),'LineWidth',1);
pLeg(2) = xline(X_const_inf(1,:),'--k','LineWidth',1);
yline(X_const_inf(2,:),'--k','LineWidth',1);
xline(X_const_sup(1,:),'--k','LineWidth',1);
yline(X_const_sup(2,:),'--k','LineWidth',1);
pLeg(4) = plot(x_t(1,1:endtime),x_t(2,1:endtime),'*-','color', col(1,: ),'LineWidth',1);

for i = 1 : endtime
    pl1 = plot(R_ZPC{i},[1 2],'r','LineWidth',1);
    pl1.Color = col(2,: );
    hold on
    pl2 =plot(R{i},[1 2],'b','LineWidth',1);
    pl2.Color = col(1,: );
end


ax = gca;
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
legend1 = legend(pLeg,{'setpoint','constraints','ZPC','TZPC'},'Interpreter','latex', 'NumColumns', 2);
set(legend1, 'Position',[0.154285262559308 0.829345455747662 0.42857187906901 0.12338181697961],...
    'NumColumns',2);

ylim([-2.1 3]);
xlim([-7.8 0.7]);
annotation(fig4,'arrow',[0.822380952380952 0.883333333333333],...
    [0.403242424242424 0.464848484848485],'LineWidth',1);
annotation(fig4,'rectangle',...
    [0.888142857142857 0.470909090909091 0.0447142857142856 0.0630303030303029],...
    'LineWidth',1);

ax2 = axes('Position',[0.45 0.21 0.35 0.30],'Color','none');
plot(0,0,'ok','LineWidth',2,'MarkerFaceColor',[1 0.5 0]);
hold on;
plot(x_t_zpc(1,1:endtime),x_t_zpc(2,1:endtime),'v-','color', col(2,: ),'LineWidth',1);
plot(x_t(1,1:endtime),x_t(2,1:endtime),'*-','color', col(1,: ),'LineWidth',1);
for i = 1 : endtime
    pl1 = plot(R_ZPC{i},[1 2],'r','LineWidth',1);
    pl1.Color = col(2,: );
    hold on;
    pl2 = plot(R{i},[1 2],'b','LineWidth',1);
    pl2.Color = col(1,: );
end
ax2.XLim = [-0.3,0.06];
ax2.YLim = [-0.03,0.17];

box off;
grid on;
ax.FontSize = Fontsizeax;
ax.XLabel.FontSize= 27;
ax.YLabel.FontSize= 27;

outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%%
if savefig == 1
    saveas(fig1, ['U-TZPC-ZPC-N' num2str(N) '.fig']);
    saveas(fig2, ['X1-TZPC-ZPC-N' num2str(N) '.fig']);
    saveas(fig3, ['X2-TZPC-ZPC-N' num2str(N) '.fig']);
    saveas(fig4, ['R-TZPC-ZPC-N' num2str(N) '.fig']);
    
    saveas(fig1, ['U-TZPC-ZPC-N' num2str(N) '.png']);
    saveas(fig2, ['X1-TZPC-ZPC-N' num2str(N) '.png']);
    saveas(fig3, ['X2-TZPC-ZPC-N' num2str(N) '.png']);
    saveas(fig4, ['R-TZPC-ZPC-N' num2str(N) '.png']);
    
    saveas(fig1, ['U-TZPC-ZPC-N' num2str(N) '.eps'],'epsc');
    saveas(fig2, ['X1-TZPC-ZPC-N' num2str(N) '.eps'],'epsc');
    saveas(fig3, ['X2-TZPC-ZPC-N' num2str(N) '.eps'],'epsc');
    saveas(fig4, ['R-TZPC-ZPC-N' num2str(N) '.eps'],'epsc');
end
