close all
clear
clc

addpath("./classes")
addpath("./abstraction")

workpoint = struct('x', [1.150295897591316e+04, 1.101600093823581e+03], 'u', 90, 'y', 36);

workpoint1 = struct('x', 1.0e+03 * [1.8405    0.0282], 'u', 18, 'y', 5.76);

sim_length = 4000;
jumpK      = 100;

jumps  = -1:0.2:1;
uJumps = 0:18:180;

y       = cell(size(jumps));
ylin    = cell(size(jumps));
legends = cell(size(jumps));

ystat    = zeros(size(jumps));
ystatlin = zeros(size(jumps));

tanks       = TankSystem(workpoint);
linearTanks = LinearTankSystem(workpoint1);
	
for i = 1:size(uJumps, 2)
	uJump = uJumps(i)
	legends{i} = "F_{1in} = " + num2str(uJump);
	
	tanks.resetToWorkPoint(workpoint);
	linearTanks.resetToWorkPoint(workpoint1);
	
	u = workpoint.u.*ones(1, sim_length);
	u(1, jumpK:end) = uJump;
	
	u1 = workpoint1.u.*ones(1, sim_length);
	u1(1, jumpK:end) = uJump;
	
	y{i}    = workpoint.y.*ones(1, sim_length);
	ylin{i} = workpoint.y.*ones(1, sim_length);
	
	for k = 1:1:sim_length
		y{i}(k)    = tanks.getOutput();
		ylin{i}(k) = linearTanks.getOutput();
		
		tanks.setControl(u(k));
		linearTanks.setControl(u1(k));
		
		tanks.nextIteration();
		linearTanks.nextIteration();
	end
	tanks.x(end, :)
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
	title ("Przebiegi wyj�cia modeli dla r�nych skok�w sterowania")
	xlabel("t [s]")
	ylabel("h_2 [cm]")
	legend(legends, 'Location', 'EastOutside')

figure
	grid on
	hold on
	plot(uJumps, ystat);
	plot(uJumps, ystatlin);
	xlabel("F_{1in} [cm^3/s]");
	ylabel("h_2 [cm]");
	title("Charakterystyka statyczna");
	legend("Model nieliniowy", "Model liniowy")
