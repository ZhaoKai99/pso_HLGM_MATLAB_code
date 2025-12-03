function [SMAPE,M0_hat,a,b,mape_all]=HLGM(r)
global pf ph p M
M0=M(1:ph,:); % Original non-negative observation sequence
[n,~]=size(M0);
for i=1:n+pf 
    for j=1:n+pf
        L0(i,j) = (1/(1+exp(-r(1)*j)))^( (i-j)*r(2)); % HLAGO,r(1)=xi,r(2)=rho
    end
end
L=tril(L0); % HLAGO  matrix
Mr=L(1:n,1:n)*M0; % HLAGO  sequence
Z1=(Mr(1:end-1)+Mr(2:end))/2; % Background value
Y=Mr(2:end,:)-Mr(1:end-1,:); 
Q=[-Z1(:,1),ones(n-1,1)];     % Q matrix
P=Y(:,1);                     % P matrix
Par=regress(P,Q);             % OLS algorithm for parameter calculation
a=Par(1);
b=Par(2);
Mr_hat(1)=M0(1);
for c=1:n+pf-1
    Mr_hat(c+1,1)=(M0(1)-b/a)*exp(-a*c)+b/a; % Predicted value of time response function 
end
M0_hat = inv(L)*Mr_hat;       % Prediction value restoration

%% Calculate MAPE error 
SMAPE=calculate_error(M(2:p,1),M0_hat(2:p,1));         % SMAPE
FMAPE=calculate_error(M(p+1:end,1),M0_hat(p+1:end,1)); % FMAPE
mape_all = [SMAPE,FMAPE];
end
