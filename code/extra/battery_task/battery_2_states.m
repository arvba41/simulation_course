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
              y(1)/batt.Csei - y(2)/(batt.Rsei*batt.Csei)];
          
fy22 = @(t,y) -y/(batt.Rsei*batt.Csei) + ibatt(t)/batt.Csei;

%% solutions
ts = 0; tf = 3600;
y021 = [0; 0]; y022 = 0; % initital conditions

% options = odeset('Stats','on','MaxStep',1e-3);
% options = odeset('Stats','on');
% tic; [t21_ode23,y21_ode23] = ode23(@(t,y) fy21(t,y),[ts tf],y021); toc;% ode23 non-stiff solver
% tic; [t22_ode23,y22_ode23] = ode23(@(t,y) fy22(t,y),[ts tf],y022); toc;% ode23 non-stiff solver

tic; [t21_ode15s,y21_ode15s] = ode15s(@(t,y) fy21(t,y),[ts tf],y021); toc;% ode15s stiff solver
tic; [t22_ode15s,y22_ode15s] = ode15s(@(t,y) fy22(t,y),[ts tf],y022); toc;% ode15s stiff solver

%% plots

load  SPICE_output_2_states_vbatt
t_SPICE = SPICEoutput2statesvbatt(:,1);
Ibatt_SPICE = SPICEoutput2statesvbatt(:,2);
Vsei_SPICE = SPICEoutput2statesvbatt(:,3);
Vbatt_SPICE = SPICEoutput2statesvbatt(:,4);
figure(1) % case 1 plots
tiledlayout(2,2)

nexttile;
plot(t21_ode15s,y21_ode15s(:,1),t_SPICE,Ibatt_SPICE,'--')
grid on
ylabel('$i_{batt}(t)$'); xlabel('$t$')
legend('ode15s','SPICE')
ax1 = figtex(gca,1);

nexttile;
plot(t21_ode15s,vbatt(t21_ode15s),t_SPICE,Vbatt_SPICE,'--')
grid on
ylabel('$v_{batt}(t)$'); xlabel('$t$')
ax2 = figtex(gca);

nexttile;
plot(t21_ode15s,y21_ode15s(:,2),t_SPICE,Vsei_SPICE,'--')
grid on
ylabel('$v_{sel}(t)$'); xlabel('$t$')
ax3 = figtex(gca);

nexttile;
semilogy(t21_ode15s(2:end),diff(t21_ode15s),t_SPICE(2:end),diff(t_SPICE),'--')
grid on
ylabel('$h$'); xlabel('$t$')
ax4 = figtex(gca);

linkaxes([ax1 ax2 ax3 ax4],'x');

figure(2) % case 1 plots
tiledlayout(2,2)

nexttile;
plot(t22_ode15s,ibatt(t22_ode15s))
grid on
ylabel('$i_{batt}(t)$'); xlabel('$t$')
ax1 = figtex(gca);

nexttile;
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

nexttile;
semilogy(t22_ode15s(2:end),diff(t22_ode15s))
grid on
ylabel('$h$'); xlabel('$t$')
ax4 = figtex(gca);

linkaxes([ax1 ax2 ax3 ax4],'x');

%% functions
function it = mmccurr(t,iamp)
    it = iamp*(square(2*pi*1/3600*t,50) + sin(2*pi*100*t));
end