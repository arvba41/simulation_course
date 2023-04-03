function dn = ef_diff(fy,t,y,h)
% Difference operator using the Euler forward method
dn = zeros(size(y));
for ii=2:length(t)
    dn(ii) = abs((y(ii) - y(ii-1))/h - fy(t(ii-1),y(ii-1))); % EF difference opeartor
end