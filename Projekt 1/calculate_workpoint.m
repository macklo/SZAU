function [workpoint] = calculate_workpoint(h2)
%CALCULATE_WORKPOINT Summary of this function goes here
%   Detailed explanation goes here
A1    = 540;
C2    = 0.85;
alfa1 = 26;
alfa2 = 20;
FD = 30;
h1=(alfa2/alfa1)^2*h2;
workpoint = struct('x', [A1*h1, C2*h2^2], 'u', alfa1*sqrt(h1)-FD, 'y', h2);
end

