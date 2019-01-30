function [ X ] = my_expsmooth( X, dt, tau )
% MY_EXPSMOOTH Exponential smoothing of vector X. 

% T(ms)
T = 1000*dt;

for m = 1:size(X,2)  
    for n = 2:size(X,1)        
        X(n,m) = exp(-T/tau)*X(n-1,m) + (1-exp(-T/tau))*X(n,m);
    end    
end