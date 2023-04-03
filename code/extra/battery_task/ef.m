function y = ef(fy,y0,t,h)
% Function to determine y, given ODE fy = f(t,y) at initital value y0 
% for a given time vector t with step size h.
% Forward Euler method
% y = euler_f(fy,t,h,y0)
% 
% y     --> y_n
% y0    --> y(0)
% fyn   --> y'(t) = f(t,y), inline function
% t     --> time vector
% h     --> step size
% 
% --------------------------------------------------------------------
% example: y'(t) = -10*y - sin(t), @ y(0) = 1, 
% 
% y0 = 1; % initital value 
% h = 0.01; % time step
% t = 0:h:1; % time vector
% y = euler_f(@(t,y) -10*y - sin(t), y0, t, h)

y(:,1) = y0; % initital condition 
% Euler forward method
for ii=2:length(t)
    y(:,ii) = y(:,ii-1) + h.*fy(t(ii-1),y(:,ii-1));
end
