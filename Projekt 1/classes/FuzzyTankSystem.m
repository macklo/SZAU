classdef FuzzyTankSystem < AbstractObject
	properties
		A1    = 540;
		C2    = 0.85;
		alfa1 = 26;
		alfa2 = 20;
		
		F1  = 90;
		FD  = 30;
		tau = 100;
		h2  = 36;
		
		u0 = 90;
		y0 = 36
		
		u, y, x
		uk, yk, xk
		
		ymin = 2.25;
		ymax = 110.2;
		
		x0 = []
		t = 0
		linearModels
		numberOfModels
		linearWorkpoints
		c
	end
	
	methods
		function self = FuzzyTankSystem(workpoint, numberOfModels)
			ny = 1; nu = 1; nd = 0; Ts = 1;
			self@AbstractObject(ny, nu, nd, Ts);
			
			self.x0 = workpoint.x;
			self.numberOfModels = numberOfModels;
			self.linearWorkpoints = cell(3, 1);
			dy = (self.ymax - self.ymin) / numberOfModels;
			a = 0.5;
			c = self.ymin+dy:dy:self.ymax-dy;
			h2r0 = ones(1, numberOfModels);
			h2r0(1) = (c(1)+self.ymin)/2-1;
			h2r0(numberOfModels) = min((self.ymax+c(numberOfModels-1))/2+1, self.ymax);
			if numberOfModels > 2
				h2r0(2:numberOfModels-1) = (c(2:numberOfModels-1)+c(1:numberOfModels-2))./2;
			end
			
			self.linearModels = cell(numberOfModels, 1);
			for i = 1:numberOfModels
				self.linearWorkpoints{i} = self.calculateWorkpoint(h2r0(i));
				self.linearModels{i} = LinearTankSystem3(self.linearWorkpoints{i});
				self.linearModels{i}.calculateWorkpoint(self.linearWorkpoints{i}.u);
				self.linearModels{i}.resetToWorkPoint(workpoint);
			end
			self.c = c;
			
			figure
			hold on
			for r = 1:numberOfModels
				if r == 1
					plot(self.ymin:0.1:self.ymax, sigmf(self.ymin:0.1:self.ymax , [-a c(1)]))
				elseif r == numberOfModels
					plot(self.ymin:0.1:self.ymax,sigmf(self.ymin:0.1:self.ymax , [a c(numberOfModels-1)]))
				else
					plot(self.ymin:0.1:self.ymax,dsigmf(self.ymin:0.1:self.ymax, [a c(r-1) a c(r)]))
				end
			end
			plot(h2r0, ones(1,numberOfModels), 'ko')
			xlabel('h2')
			ylabel('przynale?no??')
			title('Funkcje przynale?no?ci regulator?w lokalnych w regulacji rozmytej wzgl?dem warto?ci wyj?cia')
			
		end
		
		function output = getOutput(self)
			output = self.yk;
		end
		
		function setControl(self, control)
			for r = 1:self.numberOfModels
				self.linearModels{r}.setControl(control);
			end
			self.uk = control;
		end
			
		function nextIteration(self)
			for r = 1:self.numberOfModels
				self.linearModels{r}.shiftArrays();
				self.linearModels{r}.simulate();
			end
			self.shiftArrays();
			self.simulate();
		end
		
		function shiftArrays(self)
			self.u = circshift(self.u, [0 1]);
			self.y = circshift(self.y, [0 1]);
			
			self.u(:, 1) = self.uk;
			self.y(:, 1) = self.yk;
		end
		
% 		function simulate(self)
% 			lin_x0 = self.yk;
% 			lin_c0 = self.uk;
% 			
% 			x = self.A * lin_x0 + self.B*lin_c0;
% 			self.yk = self.C * x;
% 		end
		
		function resetToWorkPoint(self, workPoint)
			self.uk = workPoint.u;
			self.yk = workPoint.y; 
			self.xk = workPoint.x;

			self.u = self.uk*ones(1, self.tau);
			self.y = self.yk*ones(1, self.tau);
			self.x = (self.xk'.*ones(2, self.tau))';
			for r = 1:self.numberOfModels
				self.linearModels{r}.resetToWorkPoint(workPoint);
			end
		end
		
		function u = getU(self)
			time = size(self.u, 2) * self.Ts;
			if (time > self.tau)
				u = self.u(time - self.Ts);
			else
				u = self.u0;
			end
		end
		
		function w = getWeights(self)
			a = 0.5;
			w = zeros(self.numberOfModels, 1);

			for r = 1:self.numberOfModels
				if r == 1
                    w(r) = sigmf(self.yk , [-a self.c(1)]);
                elseif r == self.numberOfModels
                    w(r) = sigmf(self.yk , [a self.c(self.numberOfModels-1)]);
                else
                    w(r) = dsigmf(self.yk, [a self.c(r-1) a self.c(r)]);
				end
			end
		end
		
		function ykr = getModelOutputs(self)
			ykr = zeros(self.numberOfModels, 1);
			self.yk;
			for r = 1:self.numberOfModels
				ykr(r) = self.linearModels{r}.getOutput();
			end
		end
		
		function [y, t] = simulate(self)
			a = 0.5;
			w = zeros(self.numberOfModels, 1);
			ykr = zeros(self.numberOfModels, 1);
			self.yk;
			for r = 1:self.numberOfModels
				ykr(r) = self.linearModels{r}.getOutput();
				if r == 1
                    w(r) = sigmf(self.yk , [-a self.c(1)]);
                elseif r == self.numberOfModels
                    w(r) = sigmf(self.yk , [a self.c(self.numberOfModels-1)]);
                else
                    w(r) = dsigmf(self.yk, [a self.c(r-1) a self.c(r)]);
				end
			end
% 			ykr
			self.yk = (w' * ykr) / sum(w);
		end
		
% 		function xret = getX(obj, x)
% 			xret = x;
% 			if x(1) < 0
% 				xret(1) = 0;
% 			end
% 			if x(2) < 0
% 				xret(2) = 0;
% 			end
% 		end
		
		function [t, x] = simulateODE(self, x0, u0, t0, tfinal)
			[t, x] = ode45(@(t, x) self.differential(x, u0), t0:self.Ts:tfinal, x0);
		end
		
		function dx = differential(self, x, u)
			dx1 = u + self.FD - self.alfa1 * sqrt(x(1)/self.A1);
			dx2 = self.alfa1 * sqrt(x(1)/self.A1) - self.alfa2 * (x(2)/self.C2)^(1/4);
			
			dx = [dx1, dx2];
		end
		
		function workpoint = calculateWorkpoint(self, h2)
			h1 = (self.alfa2 / self.alfa1)^2 * h2;
			workpoint = struct('x', [self.A1 * h1, self.C2 * h2^2], 'u', self.alfa1 * sqrt(h1) - self.FD, 'y', h2, 'h1', h1, 'h2', h2);
		end
	end
end

