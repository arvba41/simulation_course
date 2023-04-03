clc
clear all

% Motor specifications
P = 300e3; % motor power rating [W]
Vs = 440;% stator phase-phase voltage [V]
J = 1.09; % induction motor inertia [kgm^2]
np = 20; % number of pole pairs
Nr = 1000; % output rotor speed [rpm] (nominal)
cosphi = 0.9; % power factor

phi = acos(cosphi); % phase angle
Sbase = P/cosphi; % apparent power

Ibase_mag = Sbase/(Vs*sqrt(3)); % output current magnitude
Ibase_ang = -phi; % output current phase angle
Ibase = Ibase_mag*exp(1j*Ibase_ang); % complex current

f = np*Nr/60; % electrical frequency

Zbase = Vs/(Ibase*sqrt(3)); % base impedance

% T-model paramteres
Ls1 = 0.4*imag(Zbase)/(2*pi*f); % Stator inductance
% Lr1 = 0.1*imag(Zbase)/(2*pi*f); % Rotor inductance
Lm = 3*imag(Zbase)/(2*pi*f); % magnatizing inductance
Rs = 0.1*real(Zbase); % stator resistance
% Rr = 0.1*real(Zbase); % rotor resistance

% Rs = 1.5;
% Rr = 1.04;
% Lm = 0.29;
% Ls1 = 0.0079;
% Lr1 = 0.0059;


% Inverse-\Gamma model parameters
L_sigma = Ls1; % total leakage inductance
% L_M = Lm; % transformed magnetizing inductance
% R_R = Rr; % transformed rotor resistance
% b = Lm/Lr1; % transformation factor

% Torque rated (electrical)
Te_rated = P/(Nr*pi/30);


%% Motor parameters from the Book pp 175
% Po = 1.6e3; % output power [W]
% Vs = 220; % stator nominal voltage [V]
% Nr = 1430; % nominal rotor speed [rpm]
% f = 50; % stator frequency [Hz]
% cosphi = 0.76; % power factor 
% R_R = 4.1;
% Rs = 4.6;
% L_sigma = 33e-3; 
% L_M = 0.32;
% 
% wm = pi*Nr/30;
% wn = 2*pi*f;
% np = double(int8(wn/wm));