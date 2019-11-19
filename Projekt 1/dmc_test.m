clear;
close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

workpoint = struct('x', [1.150295897591316e+04, 1.101600093823581e+03], 'u', 90, 'y', 36);

tanks = LinearTankSystem(workpoint);
tanks.resetToWorkPoint(workpoint);

umin = 0;
umax = 1000;
dumax = 200;

D = 1500;
N = 200;
Nu = 150;
lambda = 0.5;
psii = 1;
sim_length = 4000;
% load("data/setPoints.mat")
setPoints = build_random_setpoints_array(workpoint, sim_length, 500, 10 , 100);

load('./data/s.mat', 's');
reg = DMC_Regulator(tanks, workpoint, s, D, N, Nu, lambda, psii, umin, umax, dumax);
% reg = Numeric_DMC_Regulator(tanks, workpoint, s, D, N, Nu, lambda, psii, umin, umax, dumax);

u = workpoint.u.*ones(tanks.nu, sim_length);
y = workpoint.y.*ones(tanks.ny, sim_length);

for k = 1:sim_length
	k
    output = tanks.getOutput();
    y(:, k) = output;
    control = reg.calculate(output, setPoints(:, k));
    u(:, k) = control';
    tanks.setControl(control);
    tanks.nextIteration();
end

figure;
	stairs(u, 'r');

figure;
	hold on;
	stairs(setPoints, 'b');
	stairs(y, 'r');


e = (y - setPoints)*(y - setPoints)';