%% Instructions
% High-order Logistic grey model: HLGM(1,1)
% This code is the main function of the pso-HLGM(1,1) model.
% Case 8: Predicted of renewable energy consumption in Poland.

%% Start
clear; clc; close all; warning off;
global error_style pf ph p M
p = 9;  % Sample size
M = [1	1.4	1.8	2.4	3.4	3.3	4	4.7	4.7	4.9	4.4]'; % Original non-negative observation sequence
pf = 2; % Out of sample predictions
ph = p+2-pf; % Training sets
error_style='MAPE';
[mape,M0_hat,a,b,mape_all]=HLGM([-0.946689549	-0.045953304]); 

%% Drawing
plot(M,'b-*','LineWidth',1.2,'MarkerSize',6,'DisplayName','Raw data');hold on;
plot(M0_hat,'r-o','LineWidth',1.2,'MarkerSize',8,'DisplayName','Fitting value');
% legend('Raw data','Fitting value')

