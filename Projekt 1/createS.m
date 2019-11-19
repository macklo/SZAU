clear;
close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

workpoint = struct('x', [1.150295897591316e+04, 1.101600093823581e+03], 'u', 90, 'y', 36);
tanks = LinearTankSystem(workpoint);
tanks.resetToWorkPoint(workpoint);

sim_length = 1510;

s = cell(tanks.ny, tanks.nu);

start = 10;

for n = 1:tanks.nu
    u = workpoint.u.*ones(tanks.nu, sim_length);
    y = workpoint.y.*ones(tanks.ny, sim_length);
	
    u(n, start:end) = workpoint.u(n) + 1;
    for k = 1:sim_length
        y(:, k) = tanks.getOutput();
        tanks.setControl(u(:, k));
        tanks.nextIteration();
    end
    tanks.resetToWorkPoint(workpoint);
    for m = 1:tanks.ny
        s{m, n} = y(m, start+1:end) - y(m, start);
    end
end

save('./data/s.mat', 's');
figure;
    stairs(s{1, 1});