close all
clear
clc

addpath("./classes")
addpath("./abstraction")

workpoint = struct('x', [1.150295897591316e+04, 1.101600093823581e+03], 'u', 90, 'y', 36);

sim_length = 4000;
jumpK = 100;

legends = [];
y = [];
ystat = [];
ystatlin = []';

jumps = -1:0.2:1;
uJumps = workpoint.u + jumps * workpoint.u;

figure
	hold on
	grid on

tanks       = TankSystem(workpoint);
linearTanks = LinearTankSystem(workpoint);
	
for uJump = uJumps
	uJump
	
	tanks.resetToWorkPoint(workpoint);
	linearTanks.resetToWorkPoint(workpoint);
	
	u = workpoint.u.*ones(1, sim_length);
	u(1, jumpK:end) = uJump;
	
	y = workpoint.y.*ones(1, sim_length);
	ylin = workpoint.y.*ones(1, sim_length);
	
	for k = 1:1:sim_length
		y(k) = tanks.getOutput();
		ylin(k) = linearTanks.getOutput();
		
		tanks.setControl(u(k));
		linearTanks.setControl(u(k));
		
		tanks.nextIteration();
		linearTanks.nextIteration();
	end

	stairs(y)
	stairs(ylin, '--')
	ystat = [ystat y(end)];
	ystatlin = [ystatlin ylin(end)];
end

figure
	grid on
	hold on
	plot(uJumps, ystat, '-o');
	plot(uJumps, ystatlin, '-o');
	xlabel("F_1 [cm^3/s]");
	ylabel("h_2 [cm]");
	title("Charakterystyka statyczna");
