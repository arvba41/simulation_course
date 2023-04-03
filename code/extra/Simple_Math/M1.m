clear all; clc

plot_vector_field = 1;

fy = @(t,y) -sign(y)*sqrt(abs(y));

y0 = 10;
ts = 0;
tf = 10;

tic; [t_ode45,y_ode45] = ode45(fy,[ts tf],y0);

figure(1)
clf 

if plot_vector_field
    [T_mtx,Y_mtx] = meshgrid(-10:1e-1:10,-10:1e-1:10);
    FY_mtx = fy(T_mtx,Y_mtx);
    dT = ones(size(T_mtx));
    L = sqrt(1+FY_mtx.^2);
    quiver(T_mtx, Y_mtx, dT./L, FY_mtx./L, 0.5)
    axis tight
end

hold on
plot(t_ode45,y_ode45)

