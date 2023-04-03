function y = euler_b(fy,y0,t,h)
% Function to determine y, given ODE fy = f(y,t) at initital value y0 
% for a given time vector t with step size h.
% Backwards Euler method
% y = euler_f(fy,t,h,y0)
% 
% y     --> y_n
% y0    --> y(0)
% fyn   --> y'(t) = f(y,t), inline function
% t     --> time vector
% h     --> step size
% 
% --------------------------------------------------------------------
% example: y'(t) = -10*y @ y(0) = 1, 
% 
% y0 = 1; % initital value 
% h = 0.01; % time step
% t = 0:h:1; % time vector
% y = euler_f(@(y) -10*y, y0, t, h)

y(1) = y0; % initital condition 
% Euler forward method
for ii=2:length(t)
    y(ii) = fsolve(@(x) y(ii-1) + h.*fy(x) - x,0); % Trust-Region-Dogleg
end
