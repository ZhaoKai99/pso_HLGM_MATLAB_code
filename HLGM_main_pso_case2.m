%% Instructions
% High-order Logistic grey model: HLGM(1,1)
% This code is the main function of the pso-HLGM(1,1) model.
% This code is used for sensitivity analysis experiments on sample size.
% The sensitivity of the HLGM(1,1) model can be analyzed by changing the size of the input sample size p (belonging to 5, 6, ... ,10).
% Case 2: Australia  carbon  emission data prediction and HLGM  sensitivity  analysis.
% Reference: [1] Cai, K., Wu, L., 2024. Grey prediction of carbon emission and carbon peak in several developing countries. Engineering Applications of Artificial Intelligence  133, 108210.

%% Start
clear; clc; close all; warning off;
global error_style pf ph p M
p = 9;  % Sample size
M = [363.8; 370.2; 383.2; 383.0; 398.6; 405.8; 415.2; 405.6; 396.5; 400.7; 393.1]; % Original non-negative observation sequence
pf = 6; % Out of sample predictions
ph = 5; % Training sets
error_style='MAPE';

%% PSO algorithm
tic % Start timing
n = 30; % Number of particles
narvs =2; % Number of variables
c1 = 2;  % Learning factors c1
c2 = 2;  % Learning factors c2
w = 0.9;  % Inertial weight
K = 400;  % Number of iterations
x_lb = [-1 -1]; % The lower bound of x
x_ub = [1 1]; % The upper bound of x
range = [x_lb',x_ub'];
vmax = [0.1,0.1]; % Maximum velocity of particles

%% Initialize the position and velocity of particles
x = zeros(n,narvs);
for i = 1: narvs
    x(:,i) = x_lb(i) + (x_ub(i)-x_lb(i))*rand(n,1);    % Randomly initialize the particle's location within the defined domain
end
v = -vmax + 2*vmax .* rand(n,narvs);  % Randomly initialize the velocity of particles (here we set it to [- vmax, vmax])

%% Calculate fitness
fit = zeros(n,1);  % Initialize the fitness of all n particles to 0
for i = 1:n  % Cycle the entire particle swarm and calculate the fitness of each particle
    fit(i) = HLGM(x(i,:));   % Call HLGM function to calculate fitness
end
pbest = x;   % Initialize the best position found so far for these n particles (which is a vector of n * narvs)
ind = find(fit == min(fit), 1);  % Find the index of the particle with the smallest fitness
gbest = x(ind,:);  % Define the best position found by all particles so far (which is a 1 * narvs vector)

%% Mark the positions of these n particles on the graph for demonstration purposes
% h = scatter3(x(:,1),x(:,2),fit,'*r');

%% Iterate K times to update speed and position
fitnessbest = ones(K,1);  % Initialize the best fitness obtained from each iteration
for d = 1:K  % Start iteration, a total of K iterations
    disp(d)
    for i = 1:n   % Update the velocity and position of the i-th particle in sequence
        v(i,:) = w*v(i,:) + c1*rand(1)*(pbest(i,:) - x(i,:)) + c2*rand(1)*(gbest - x(i,:));  % Update the velocity of the i-th particle
        % If the particle's velocity exceeds the maximum velocity limit, adjust it accordingly
        for j = 1: narvs
            if v(i,j) < -vmax(j)
                v(i,j) = -vmax(j);
            elseif v(i,j) > vmax(j)
                v(i,j) = vmax(j);
            end
        end
        x(i,:) = x(i,:) + v(i,:); % Update the velocity of the i-th particle
        % If the position of the particle exceeds the domain, adjust it
        for j = 1: narvs
            if x(i,j) < x_lb(j)
                x(i,j) = x_lb(j);
            elseif x(i,j) > x_ub(j)
                x(i,j) = x_ub(j);
            end
        end
        fit(i) = HLGM(x(i,:));  % Recalculate the fitness of the i-th particle
        if fit(i) < HLGM(pbest(i,:))   % If the fitness of the i-th particle is less than the fitness corresponding to the best position found by this particle so far
            pbest(i,:) = x(i,:);   % Then update the best position found so far for the i-th particle
        end
        if  fit(i) < HLGM(gbest)  % If the fitness of the i-th particle is less than the fitness corresponding to the best position found by all particles so far
            gbest = pbest(i,:);   % Then update the best positions found by all particles so far
        end
    end
    fitnessbest(d) = HLGM(gbest);  % Update the best fitness obtained from the d-th iteration
    h.XData = x(:,1);
    h.YData = x(:,2);
    h.ZData = fit;
end
toc % End timing

%% output result
disp('The best position is: '); disp(gbest)
disp('The optimal value at this point is: '); disp(HLGM(gbest))
[mape,M0_hat,a,b,mape_all]=HLGM(gbest(1:2));
% [mape,M0_hat,a,b,mape_all]=HLGM([-0.319082784	-0.019734954]); % Sample size p=10
% [mape,M0_hat,a,b,mape_all]=HLGM([-0.208662999	-0.034588604]); % Sample size p=9
% [mape,M0_hat,a,b,mape_all]=HLGM([-0.236615319	-0.029291137]); % Sample size p=8
% [mape,M0_hat,a,b,mape_all]=HLGM([-0.055828854	0.019855285]); % Sample size p=7
% [mape,M0_hat,a,b,mape_all]=HLGM([-0.041079875	0.35342762]); % Sample size p=6
% [mape,M0_hat,a,b,mape_all]=HLGM([1.0000	0.075561859]); % Sample size p=5

%% Drawing
figure
subplot(2,1,1) % fig1
plot(fitnessbest,LineWidth=2)  % Draw a graph of the best fitness change for each iteration
xlabel('Iteration count');
legend('PSO algorithm')
subplot(2,1,2) % fig2
plot(M,'b-*','LineWidth',1.2,'MarkerSize',6,'DisplayName','Raw data');hold on;
plot(M0_hat,'r-o','LineWidth',1.2,'MarkerSize',8,'DisplayName','Fitting value');
legend('Raw data','Fitting value')

