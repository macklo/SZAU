clear;
% close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

workpoint = struct('x', [1.150295897591316e+04, 1.101600093823581e+03], 'u', 90, 'y', 36);
workpointD = struct('y', 30);

tanks = TankSystem(workpoint);
tanks.resetToWorkPoint(workpoint);

umin = 0;
umax = 180;
dumax = 1;

D = 1500;
N = 500;
Nu = 500;
lambda = 10;
psii = 1;
sim_length = 10000;
load("data/setPointsY.mat")
load("data/d.mat")
setPoints = setPoints(1:sim_length);
% setPointsY = build_random_setpoints_array(workpoint, sim_length, 1000, 2000, 6 , 66);
% d = build_random_setpoints_array(workpoint, sim_length, 2000, 2000, 0 , 60);

load('./data/s.mat', 's');
% load('./data/fuzzyS.mat', 'fuzzyS')
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
	tanks.setDisturbance(d(k));
    tanks.nextIteration();
end

figure;
	stairs(u, 'r');

figure;
	stairs(d, 'r');

figure;
	hold on;
	stairs(setPoints, 'b');
	stairs(y, 'r');


e = (y - setPoints)*(y - setPoints)';