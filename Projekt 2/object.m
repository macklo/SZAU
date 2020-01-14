function [x1,x2,y] = object(um4,x1m1,x2m1)
%OBJECT Calculates objects output
%   Detailed explanation goes here

alfa1 = -1.473409;
alfa2 = 0.525788;
beta1 = 0.026085;
beta2 = 0.021057;
noise = 0.00;

g1 = (exp(6*um4) - 1)/(exp(6*um4) + 1);
x1 = -alfa1*x1m1 + x2m1 + beta1*g1;
x2 = -alfa2*x1m1 + beta2*g1;
y = -0.5*(1-exp(-2*x1)) - noise + 2*noise*rand;

end

