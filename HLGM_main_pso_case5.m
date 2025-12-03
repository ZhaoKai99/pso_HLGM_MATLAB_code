%% Instructions
% High-order Logistic grey model: HLGM(1,1)
% This code is the main function of the pso-HLGM(1,1) model.
% Case 5: The proportion of renewable energy consumption in China is predicted.


%% Start
clear; clc; close all; warning off;
global error_style pf ph p M
p = 6;  % Sample size
M = [13.43	12.26	11.34	11.54	11.52	12.06	12.25	12.59	12.86	13.12]'; % Original non-negative observation sequence
pf = 4; % Out of sample predictions
ph = p+4-pf; % Training sets
error_style='MAPE';
[mape,M0_hat,a,b,mape_all]=HLGM([0.64540465718	-0.3697116732]); % case3 

%% Drawing
plot(M,'b-*','LineWidth',1.2,'MarkerSize',6,'DisplayName','Raw data');hold on;
plot(M0_hat,'r-o','LineWidth',1.2,'MarkerSize',8,'DisplayName','Fitting value');
legend('Raw data','Fitting value')

