%% Instructions
% High-order Logistic grey model: HLGM(1,1)
% This code is the main function of the pso-HLGM(1,1) model.
% Case 6: The forecasting results in publication output analysis for USA.


%% Start
clear; clc; close all; warning off;
global error_style pf ph p M
p = 8;  % Sample size
M = [485892	505691	505726	515346	551036	577686	600611	622045	626583	621753]'; % Original non-negative observation sequence
pf = 2; % Out of sample predictions
ph = p+2-pf; % Training sets
error_style='MAPE';
[mape,M0_hat,a,b,mape_all]=HLGM([0.9802	-0.717493956]); 

%% Drawing
plot(M,'b-*','LineWidth',1.2,'MarkerSize',6,'DisplayName','Raw data');hold on;
plot(M0_hat,'r-o','LineWidth',1.2,'MarkerSize',8,'DisplayName','Fitting value');
% legend('Raw data','Fitting value')

