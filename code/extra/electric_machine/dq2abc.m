function [abc,alphabeta] = dq2abc(dq,t,f)
% ABC2DQ Summary of this function goes here
%   Detailed explanation goes here
wt = 2*pi*f*t;
for ii = 1:length(t)
    alphabeta(:,ii) = [cos(wt(ii)) -sin(wt(ii)); sin(wt(ii)) cos(wt(ii))]*dq(:,ii);
    abc(:,ii) =[1 0;-1/2 sqrt(3)/2;-1/2 -sqrt(3)/2]*alphabeta(:,ii);
end


