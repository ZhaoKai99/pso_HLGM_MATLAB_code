%% Instructions
% High-order Logistic grey model: HLGM(1,1)
% This code is the main function of the pso-HLGM(1,1) model.
% Case 3: Predicted electricity consumption in China.

%% Start
clear; clc; close all; warning off;
global error_style pf ph p M
p = 13;  % Sample size
M = [1480.8	1654	1910.5	2203.3	2500.2	2865.7	3281.5	3495.7	3714.6	4192.3	4692.8	4959.1	5322.3	5523.3]'; % Original non-negative observation sequence
pf = 1; % Out of sample predictions
ph = p+1-pf; % Training sets
error_style='MAPE';
[mape,M0_hat,a,b,mape_all]=HLGM([-0.1520 	0.5096]); % case3

%% Drawing
plot(M,'b-*','LineWidth',1.2,'MarkerSize',6,'DisplayName','Raw data');hold on;
plot(M0_hat,'r-o','LineWidth',1.2,'MarkerSize',8,'DisplayName','Fitting value');
legend('Raw data','Fitting value')

