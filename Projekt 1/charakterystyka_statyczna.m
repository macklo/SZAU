% close all
clear
clc

addpath("./classes")
addpath("./abstraction")

% workpoint = struct('x', [1.150295897591316e+04, 1.101600093823581e+03], 'u', 90, 'y', 36);
workpoint = calculate_workpoint(36)

sim_length = 3000;
jumpK      = 100;

jumps  = -1:0.2:1;
uJumps = workpoint.u + jumps * workpoint.u;

y       = cell(size(jumps));
ylin    = cell(size(jumps));
legends = cell(size(jumps));

ystat    = zeros(size(jumps));
ystatlin = zeros(size(jumps));

tanks       = TankSystem(workpoint);
linearTanks = LinearTankSystem3(workpoint);
	
for i = 1:size(uJumps, 2)
	uJump = uJumps(i)
	legends{i} = "F_{1in} = " + num2str(uJump);
	
	tanks.resetToWorkPoint(workpoint);
	linearTanks.resetToWorkPoint(workpoint);
	
	u = workpoint.u.*ones(1, sim_length);
	u(1, jumpK:end) = uJump;
	
	y{i}    = workpoint.y.*ones(1, sim_length);
	ylin{i} = workpoint.y.*ones(1, sim_length);
	
	for k = 1:1:sim_length
		y{i}(k)    = tanks.getOutput();
		ylin{i}(k) = linearTanks.getOutput();
		
		tanks.setControl(u(k));
		linearTanks.setControl(u(k));
		
		tanks.nextIteration();
		linearTanks.nextIteration();
	end
	
	ystat(i)    = y{i}(end);
	ystatlin(i) = ylin{i}(end);
end


figure
	grid on
	hold on
	for i = 1:size(uJumps, 2)
		stairs(y{i})
	end
	
	set(gca, 'ColorOrderIndex', 1)
	
	for i = 1:size(uJumps, 2)
		stairs(ylin{i}, '--')
	end
	title ("Przebiegi wyjœcia modeli dla ró¿nych skoków sterowania")
	xlabel("t [s]")
	ylabel("h_2 [cm]")
% 	legend(legends, 'Location', 'EastOutside')

figure
	grid on
	hold on
	plot(uJumps, ystat);
	plot(uJumps, ystatlin);
	xlabel("F_{1in} [cm^3/s]");
	ylabel("h_2 [cm]");
	title("Charakterystyka statyczna");
	legend("Model nieliniowy", "Model liniowy")
