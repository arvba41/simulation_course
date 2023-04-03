function [dq,alphabeta] = abc2dq(abc,t,f,K)
% ABC2DQ Summary of this function goes here
%   Detailed explanation goes here
wt = 2*pi*f*t;
for ii = 1:length(t)
    alphabeta(:,ii) = 2/3*K*[1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2]*abc(:,ii);
    dq(:,ii) = [cos(wt(ii)) sin(wt(ii)); -sin(wt(ii)) cos(wt(ii))]*alphabeta(:,ii);
end


