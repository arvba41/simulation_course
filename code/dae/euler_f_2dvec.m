function y = euler_f_2dvec(fy,t,y0)
% Function to determine y, given ODE fy = f(y,t) at initital value y0 
% for a given time vector t with step size h.
% Forward Euler method
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

h = mean(diff(t));
y(:,:,1) = y0; % initital condition 
% Euler forward method
for ii=2:length(t)
    y(:,:,ii) = y(:,:,ii-1) + h.*fy(t(ii),y(:,:,ii-1));
end
