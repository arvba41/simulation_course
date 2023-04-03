function [is,Te,wout] = em_readable_params(y,params,K)

is = abs(y(1,:) + 1j*y(2,:));

% iabc = dq2abc(y(1:2,:),t,freq);

% calculating the torques
Te = 3*params.np/(2*K^2)*((params.Ld - params.Lq)*y(1,:).*y(2,:) + y(2,:)*params.psi);

% speed conversion
wout = y(3,:)*30/pi./params.np; 