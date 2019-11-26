close all
clear
clc

addpath("./classes")
addpath("./abstraction")

workpoint = calculate_workpoint(36);

% workpoint1 = struct('x', 1.0e+03 * [1.8405    0.0282], 'u', 90, 'y', 33.43);
workpoint1 = workpoint;

sim_length = 2000;
jumpK      = 100;
numberOfModels = 5;

jumps  = -0.5:0.125:0.5;
uJumps = workpoint.u + jumps * workpoint.u;

y       = cell(size(jumps));
ylin    = cell(size(jumps));
yfuz    = cell(size(jumps));
legends = cell(size(jumps));

ystat    = zeros(size(jumps));
ystatlin = zeros(size(jumps));

tanks       = TankSystem(workpoint);
fuzzyTanks  = FuzzyTankSystem(workpoint1, numberOfModels);
linTanks    = LinearTankSystem3(workpoint);

fuzzyTanksOutputs = cell(size(jumps));
fuzzyTanksWeights = cell(size(jumps));
	
for i = 1:size(uJumps, 2)
	uJump = uJumps(i)
	legends{i} = "F_{1in} = " + num2str(uJump);
	
	tanks.resetToWorkPoint(workpoint);
	linTanks.resetToWorkPoint(workpoint);
	fuzzyTanks.resetToWorkPoint(workpoint1);
	
	u = workpoint.u.*ones(1, sim_length);
	u(1, jumpK:end) = uJump;
	
	u1 = workpoint1.u.*ones(1, sim_length);
	u1(1, jumpK:end) = uJump;
	
	y{i}    = workpoint.y.*ones(1, sim_length);
	ylin{i} = workpoint.y.*ones(1, sim_length);
	yfuz{i} = workpoint.y.*ones(1, sim_length);
	
	for k = 1:1:sim_length
		y{i}(k)    = tanks.getOutput();
		ylin{i}(k) = linTanks.getOutput();
		yfuz{i}(k) = fuzzyTanks.getOutput();
		
		fuzzyTanksOutputs{i}(:, k) = fuzzyTanks.getModelOutputs();
		fuzzyTanksWeights{i}(:, k) = fuzzyTanks.getWeights();
		
		tanks.setControl(u(k));
		linTanks.setControl(u(k));
		fuzzyTanks.setControl(u(k));
		
		tanks.nextIteration();
		linTanks.nextIteration();
		fuzzyTanks.nextIteration();
	end
	tanks.x(end, :)
	ystat(i)    = y{i}(end);
	ystatlin(i) = yfuz{i}(end);
end


figure
	grid on
	hold on
	for i = 1:size(uJumps, 2)
		plot(y{i}, 'b')
	end
	
	set(gca, 'ColorOrderIndex', 1)
	
	for i = 1:size(uJumps, 2)
		plot(yfuz{i}, 'g')
	end
	
	for i = 1:size(uJumps, 2)
		plot(ylin{i}, 'm')
	end
	xlabel("t [s]")
	ylabel("h_2 [cm]")

figure(1)
saveas(gcf, "./fig/bmf" + num2str(numberOfModels) + ".emf")

figure(2)
saveas(gcf, "./fig/2_1_" + num2str(numberOfModels) + ".emf")
% figure
% 	grid on
% 	hold on
% 	plot(uJumps, ystat, 'o');
% 	plot(uJumps, ystatlin, 'o');
% 	xlabel("F_{1in} [cm^3/s]");
% 	ylabel("h_2 [cm]");
% 	title("Charakterystyka statyczna");
% 	legend("Model nieliniowy", "Model liniowy")

% figure
% 	grid on
% 	hold on
% 	for i = 1:numberOfModels
% 		plot(fuzzyTanksOutputs{1}(i, :), 'b');
% 	end
% 	
% figure
% 	grid on
% 	hold on
% 	for i = 1:numberOfModels
% 		stairs(fuzzyTanksWeights{1}(i, :), 'g');
% 	end