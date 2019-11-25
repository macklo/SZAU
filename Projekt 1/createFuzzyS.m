function s = createFuzzyS(linPoints)
	numberOfModels = length(linPoints);
	s = cell(numberOfModels, 1); 

	sim_length = 8000;
	start = 10;

	for i = 1:numberOfModels
		workpoint = calculate_workpoint(linPoints(i));
		tanks = LinearTankSystem(workpoint);
		tanks.resetToWorkPoint(workpoint);

		u = workpoint.u.*ones(tanks.nu, sim_length);
		y = workpoint.y.*ones(tanks.ny, sim_length);

		u(1, start:end) = workpoint.u(1) + 1;

		for k = 1:sim_length
			y(:, k) = tanks.getOutput();
			tanks.setControl(u(:, k));
			tanks.nextIteration();
		end
		tanks.resetToWorkPoint(workpoint);
		for m = 1:tanks.ny
			s{i} = y(m, start+1:end) - y(m, start);
		end
	end

	figure
		hold on
		for i = 1:numberOfModels
			stairs(s{i})
		end
end

