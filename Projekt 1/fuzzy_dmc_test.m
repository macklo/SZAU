clear
% close all
clc

addpath('./abstraction')
addpath('./classes')
addpath('./..')

workpoint = struct('x', [1.150295897591316e+04, 1.101600093823581e+03], 'u', 90, 'y', 36);

tanks = TankSystem(workpoint);
tanks.resetToWorkPoint(workpoint);

umin = 30;
umax = 150;
dumax = 1;

D = 2200;
N = 500;
Nu = 500;
lambda = 50;
psii = 1;
sim_length = 15500;

numberOfModels = 5;
ymin = 16;
ymax = 66;
a = 3;

load("data/setPointsY.mat")
load("data/d.mat")
setPoints = setPoints(1:sim_length);
% setPoints = build_random_setpoints_array(workpoint, sim_length, 1000, 20 , 120);

% load('./data/s.mat', 's');
% reg = DMC_Regulator(tanks, workpoint, s, D, N, Nu, lambda, psii, umin, umax, dumax);

% [mf, linPoints] = createMembershipFunction(numberOfModels, ymin, ymax, a);
[mf, linPoints] = createMembershipFunctionFromCuts([25, 38, 45, 60], ymin, ymax, a);
fuzzyS = createFuzzyS(linPoints);
localRegs = cell(numberOfModels, 1);

for  i = 1:numberOfModels
	localRegs{i} = Local_DMC_Regulator(tanks, workpoint, {fuzzyS{i}}, D, N, Nu, lambda, psii);
end

fuzzyReg = Fuzzy_Regulator(mf, localRegs, umin, umax, dumax);

u = workpoint.u.*ones(tanks.nu, sim_length);
y = workpoint.y.*ones(tanks.ny, sim_length);
localControl = zeros(numberOfModels, sim_length);
weights = zeros(numberOfModels, sim_length);

for k = 1:sim_length
	k
    output = tanks.getOutput();
    y(:, k) = output;
    [control, localControl(:, k), weights(:, k)] = fuzzyReg.calculate(output, setPoints(:, k));
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

figure
	hold on
	for i = 1:numberOfModels
		stairs(localControl(i, :));
	end
	
figure
	hold on
	for i = 1:numberOfModels
		stairs(weights(i, :));
	end


e = (y - setPoints)*(y - setPoints)' / sim_length;