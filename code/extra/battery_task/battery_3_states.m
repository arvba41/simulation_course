clear all
% clc

% two cases are investigated, 
%   one, potentiostatic simulation [vbatt is known]
%   two, galvanostatic simulation [ibatt is known]

load battparams_22Age_0mesh_BOL.mat

% variables
vocv = 3.5; % open circuit voltage
iamp = 2; % output current amplitude 
vamp = 0.05; % output voltage magnitude from ocv

% inputs 
    % case 1
    vbatt = @(t) vamp/2*square(2*pi*1/1800*t,50) + vocv; 
    % case 2
    ibatt = @(t) iamp/2*square(2*pi*1/1800*t,50) + iamp/2; 
%         ibatt = @(t) mmccurr(t,iamp);
        
% function
fy21 = @(t,y) [(vocv - y(1)*batt.R0 - y(2) - vbatt(t))/batt.L0;...
              y(1)/batt.Csei - y(2)/(batt.Rsei*batt.Csei) + y(3)/(batt.Rsei*batt.Csei);...
              y(2)/(batt.Rsei*batt.Cdl) - y(3)/(batt.Cdl)*(1/batt.Rsei + 1/batt.Rct)];
          
fy22 = @(t,y) [-y(1)/(batt.Rsei*batt.Csei) + y(2)/(batt.Rsei*batt.Csei) + ibatt(t)/batt.Csei;...
              y(1)/(batt.Rsei*batt.Cdl) - y(2)/(batt.Cdl)*(1/batt.Rsei + 1/batt.Rct)];

%% solutions
ts = 0; tf = 3600;
y021 = [0; 0; 0]; y022 = [0; 0]; % initital conditions

% options = odeset('Stats','on','MaxStep',1e-3);
% options = odeset('Stats','on');
% tic; [t21_ode23,y21_ode23] = ode23(@(t,y) fy21(t,y),[ts tf],y021); toc;% ode23 non-stiff solver
% tic; [t22_ode23,y22_ode23] = ode23(@(t,y) fy22(t,y),[ts tf],y022); toc;% ode23 non-stiff solver
tic; [t21_ode15s,y21_ode15s] = ode15s(@(t,y) fy21(t,y),[ts tf],y021); toc;% ode15s stiff solver
tic; [t22_ode15s,y22_ode15s] = ode15s(@(t,y) fy22(t,y),[ts tf],y022); toc;% ode15s non-stiff solver

%% plots

load SPICE_output_3_states_vbatt
t_SPICE = SPICEoutput3statesvbatt(:,1);
Ibatt_SPICE = SPICEoutput3statesvbatt(:,2);
Vsei_SPICE = SPICEoutput3statesvbatt(:,3);
Vdl_SPICE = SPICEoutput3statesvbatt(:,4);
Vbatt_SPICE = SPICEoutput3statesvbatt(:,5);


figure(3) % case 1 plots
tiledlayout(4,2)

nexttile([2 1]);
plot(t21_ode15s,y21_ode15s(:,1),t_SPICE,Ibatt_SPICE,'--')
grid on
ylabel('$i_{batt}(t)$'); xlabel('$t$')
legend('ode15s','SPICE')
ax1 = figtex(gca,1);

nexttile([2 1]);
plot(t21_ode15s,vbatt(t21_ode15s),t_SPICE,Vbatt_SPICE,'--')
grid on
ylabel('$v_{batt}(t)$'); xlabel('$t$')
ax2 = figtex(gca);

nexttile;
plot(t21_ode15s,y21_ode15s(:,2),t_SPICE,Vsei_SPICE,'--')
grid on
ylabel('$v_{sel}(t)$'); xlabel('$t$')
ax3 = figtex(gca);

nexttile([1 2]);
semilogy(t21_ode15s(2:end),diff(t21_ode15s),t_SPICE(2:end),diff(t_SPICE),'--')
grid on
ylabel('$h$'); xlabel('$t$')
ax4 = figtex(gca);

nexttile;
plot(t21_ode15s,y21_ode15s(:,3),t_SPICE,Vdl_SPICE,'--')
grid on
ylabel('$v_{dl}(t)$'); xlabel('$t$')
ax5 = figtex(gca);

linkaxes([ax1 ax2 ax3 ax4 ax5],'x');

figure(4) % case 1 plots
tiledlayout(4,2)

nexttile([2 1]);
plot(t22_ode15s,ibatt(t22_ode15s))
grid on
ylabel('$i_{batt}(t)$'); xlabel('$t$')
ax1 = figtex(gca);

nexttile([2 1]);
dibattdt = zeros(size(ibatt(t22_ode15s)));
dibattdt(2:end) = diff(ibatt(t22_ode15s));
x = vocv - batt.R0*ibatt(t22_ode15s) - y22_ode15s(:,1);
plot(t22_ode15s,x)
grid on
ylabel('$v_{batt}(t)$'); xlabel('$t$')
ax2 = figtex(gca);

nexttile;
plot(t22_ode15s,y22_ode15s(:,1))
grid on
ylabel('$v_{sel}(t)$'); xlabel('$t$')
ax3 = figtex(gca);

nexttile([1 2]);
semilogy(t22_ode15s(2:end),diff(t22_ode15s))
grid on
ylabel('$h$'); xlabel('$t$')
ax4 = figtex(gca);

nexttile;
plot(t22_ode15s,y22_ode15s(:,2))
grid on
ylabel('$v_{dl}(t)$'); xlabel('$t$')
ax3 = figtex(gca);

linkaxes([ax1 ax2 ax3 ax4 ax5],'x');

%% functions
function it = mmccurr(t,iamp)
    it = iamp*(square(2*pi*1/3600*t,50) + sin(2*pi*100*t));
end
