clear
% close all
clc

addpath('./abstraction')
addpath('./classes')
addpath('./..')

workpoint = struct('x', [1.150295897591316e+04, 1.101600093823581e+03], 'u', 90, 'y', 36);

tanks = TankSystem(workpoint);
tanks.resetToWorkPoint(workpoint);

umin = 0;
umax = 300;
dumax = 200;

D = 3700;
N = 400;
Nu = 150;
lambda = 100;
psii = 1;
sim_length = 10000;

numberOfModels = 5;
ymin = 0;
ymax = 120;
dy = (ymax - ymin)/numberOfModels;
linPoints = dy/2:dy:ymax-dy/2;

load("data/setPoints.mat")
setPoints = setPoints(1:sim_length);
% setPoints = build_random_setpoints_array(workpoint, sim_length, 1000, 20 , 120);

% load('./data/s.mat', 's');
% reg = DMC_Regulator(tanks, workpoint, s, D, N, Nu, lambda, psii, umin, umax, dumax);

fuzzyS = createFuzzyS(linPoints);
localRegs = cell(numberOfModels, 1);
mf = createMembershipFunction(numberOfModels);

for  i = 1:numberOfModels
	localRegs{i} = DMC_Regulator(tanks, workpoint, {fuzzyS{i}}, D, N, Nu, lambda, psii, umin, umax, dumax);
end

fuzzyReg = Fuzzy_Regulator(mf, localRegs);

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
    tanks.nextIteration();
end

figure;
	stairs(u, 'r');

figure;
	hold on;
	stairs(setPoints, 'b');
	stairs(y, 'r');

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


e = (y - setPoints)*(y - setPoints)';