function f = bdf(u,t1)
% Function to determine f, where f = u' = f'(t) for a given time vector, t
% Backward finite difference method
% f = bdf(u,t)
u(1) = 0; % initital du/dt
for ii=2:length(t)
    f = (u(ii-1) - u(ii-1))/(t(ii) - t(ii-1));
end
f = f(find())