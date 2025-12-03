%% Instructions
% High-order Logistic grey model: HLGM(1,1)
% This code is the main function of the pso-HLGM(1,1) model.
% Case 7: Prediction of installed hydropower capacity in China.


%% Start
clear; clc; close all; warning off;
global error_style pf ph p M
p = 13;  % Sample size
M = [7935	8301	8607	9494	10524	11739	13029	14823	17260	19629	21606	23298	24947	28044	30486	31954	33207	34359	35259	35804]'; % Original non-negative observation sequence
pf = 7; % Out of sample predictions
ph = p+7-pf; % Training sets
error_style='MAPE';
[mape,M0_hat,a,b,mape_all]=HLGM([-0.346098397	0.231463323]); 

%% Drawing
plot(M,'b-*','LineWidth',1.2,'MarkerSize',6,'DisplayName','Raw data');hold on;
plot(M0_hat,'r-o','LineWidth',1.2,'MarkerSize',8,'DisplayName','Fitting value');
% legend('Raw data','Fitting value')

