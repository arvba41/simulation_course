function [t,y,yh] = rk45(fy,y0,tvec,h_init,reltol,abstol)
% Function to determine y, given ODE fy = f(t,y) at initital value y0 
% for a given time vector tvec = [ts tf] and step size h.
% ODE Solver: Dormand and prince 4(5) - 7 step method
% y = rk45(fy,y0,tvec,h_init,{tol})
% outputs
% t         --> time vector 
% y         --> y_n
% yh        --> \hat y_n
% inputs
% fy        --> y'(t) = f(t,y), inline function
% y0        --> y(0)
% [ts tf]   --> time vector [start_time, end_time]
% h         --> initital step size
% tol       --> relative tolerence (optional), default 1
% --------------------------------------------------------------------
% example: y'(t) = -10*y @ y(0) = 1, 
% %%% code
% y0 = 1; % initital value 
% h = 0.01; % time step
% t = 0:h:1; % time vector
% y = euler_f(@(y) -10*y, y0, t, h)

if nargin < 6
    abstol = 1e-6;
    if nargin < 5
        reltol = 1e-6;
    end
end

% ode45 Dormand and prince 4(5) method constants 
b   = [5179/57600,  0,  7571/16695, 393/640,    -92097/339200,  187/2100,   1/40]; 
bh  = [35/384,      0,  500/1113,   125/192,    -2187/6784,     11/84,      0];
c = [0; 1/5; 3/10; 4/5; 8/9; 1; 1];
a = [0,         0,          0,          0,          0,          0,      0;...
    1/5,        0,          0,          0,          0,          0,      0;...
    3/40,       9/40,       0,          0,          0,          0,      0;...  
    44/45,      -56/15,     32/9,       0,          0,          0,      0;...
    19372/6561, -25360/2187,64448/6561, -212/729,   0,          0,      0;...
    9017/3168,  -355/33,    46732/5247, 49/176,     -5103/18656,0,      0;...
    35/384,     0,          500/1113,   125/192,    -2187/6784, 11/84,  0]; 


y(:,1) = y0; % initital condition 

h = h_init; % step size
ts = tvec(1); 
tf = tvec(2);

% if length(t) <= 1
%     error('Error. \n set the step size lower than tf.');
% end

err = 0; % error init
ii = 2; % initital counter value
t(ii-1) = ts; % initital time step 
t(ii) = t(ii-1) + h; % time incriment
while t(ii) <= tf
    while err < reltol
        K1 = fy(t(ii - 1) + c(1)*h,y(:,ii-1));
        K2 = fy(t(ii - 1) + c(2)*h,y(:,ii-1) + K1*a(2,1));
        K3 = fy(t(ii - 1) + c(3)*h,y(:,ii-1) + K1*a(3,1) + K2*a(3,2));
        K4 = fy(t(ii - 1) + c(4)*h,y(:,ii-1) + K1*a(4,1) + K2*a(4,2) + K3*a(4,3));
        K5 = fy(t(ii - 1) + c(5)*h,y(:,ii-1) + K1*a(5,1) + K2*a(5,2) + K3*a(5,3) + K4*a(5,4));
        K6 = fy(t(ii - 1) + c(6)*h,y(:,ii-1) + K1*a(6,1) + K2*a(6,2) + K3*a(6,3) + K4*a(6,4) + K5*a(6,5));
        K7 = fy(t(ii - 1) + c(7)*h,y(:,ii-1) + K1*a(7,1) + K2*a(7,2) + K3*a(7,3) + K4*a(7,4) + K5*a(7,5) + K6*a(7,6));
        K = [K1, K2, K3, K4, K5, K6, K7]; % K-matrix
        y(:,ii) = y(:,ii-1) + h*sum(b.*K); % yn
        z(:,ii) = y(:,ii-1) + h*sum(bh.*K); % \hat yn

        err = abs(y(:,ii) - z(:,ii)) ; % error calculation
        
        if err < reltol
            h = (h/(2*err))^(1/5);
        end
%         if err > tol 
%             h = max(eps,h/2); % half the step-size;
%             if h == eps
%                 warning('minimum step size reached at t = %f',t(ii));
%             end
%         else
%             err_ckh = 0; % exit the inner loop
%         end
%         t(ii) = t(ii-1) + h; % time incriment
    end
    ii = ii + 1; % increment size of time vector 
    h = h_init; % reset the step size
    t(ii) = t(ii-1) + h; % time incriment
    err_ckh = 1; % reset the error check
end