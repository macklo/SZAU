clear;
% close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

workpoint = struct('x', [1.150295897591316e+04, 1.101600093823581e+03], 'u', 90, 'y', 36);

tanks = TankSystem(workpoint);
tanks.resetToWorkPoint(workpoint);

umin = 30;
umax = 150;
dumax = 1;

D = 1500;
N = 500;
Nu = 500;
lambda = 100;
psii = 1;
sim_length = 15500;
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

figure
	subplot(3, 1, 1)
		hold on;
		plot(setPoints, 'b');
		plot(y, 'r');
		ylabel("h_2[cm]")
		xlabel("t[s]")
		ylim([0 60])
	subplot(3, 1, 2)
		plot(u, 'r');
		ylabel("F_{1in}[cm^3/s]")
		xlabel("t[s]")
		ylim([30 150])
	subplot(3, 1, 3)
		plot(d, 'r');
		ylabel("F_{D}[cm^3/s]")
		xlabel("t[s]")
		ylim([15 45])

e = (y - setPoints)*(y - setPoints)' / sim_length;