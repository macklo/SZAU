classdef TankSystem < AbstractObject
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
		uk, yk, xk, dk
		
		x0 = []
		t = 0
	end
	
	methods
		function self = TankSystem(workpoint)
			ny = 1; nu = 1; nd = 0; Ts = 1;
			self@AbstractObject(ny, nu, nd, Ts);
			
			self.x0 = workpoint.x;
			self.dk = self.FD;
		end
		
		function output = getOutput(self)
			output = self.yk;
		end
		
		function setControl(self, control)
			self.uk = control;
		end
		
		function setDisturbance(self, control)
			self.dk = control;
		end
			
		function nextIteration(self)
			self.shiftArrays();
			self.simulate();
		end
		
		function shiftArrays(self)
			self.u = circshift(self.u, [0 1]);
			self.y = circshift(self.y, [0 1]);
			
			self.u(:, 1) = self.uk;
			self.y(:, 1) = self.yk;
		end
		
		
		function resetToWorkPoint(self, workPoint)
			self.uk = workPoint.u;
			self.yk = workPoint.y; 
			self.xk = workPoint.x;

			self.u = self.uk*ones(1, self.tau);
			self.y = self.yk*ones(1, self.tau);
			self.x = (self.xk'.*ones(2, self.tau))';
		end
		
		
		function [y, t] = simulate(self)
% 			self.u = [self.u, self.uk];
			if size(self.x, 1) == 0
				x = self.x0;
			else
				x = self.x(end, :);
				u = self.u(self.tau);
				
				k1 = self.Ts * self.differential(x, u);
				k2 = self.Ts * self.differential(x + k1/2, u);
				k3 = self.Ts * self.differential(x + k2/2, u);
				k4 = self.Ts * self.differential(x + k3, u);

				x = x + (k1 + 2*k2 + 2*k3 + k4)/6;
			end
			if x(1) < 0
				x(1) = 0;
			end
			if x(2) < 0
				x(2) = 0;
			end
			self.x = [self.x ; x];
			self.yk = sqrt(x(2)/self.C2);
		end
		
		
		function [t, x] = simulateODE(self, x0, u0, t0, tfinal)
			[t, x] = ode45(@(t, x) self.differential(x, u0), t0:self.Ts:tfinal, x0);
		end
		
		function dx = differential(self, x, u)
			dx1 = u + self.dk - self.alfa1 * sqrt(x(1)/self.A1);
			dx2 = self.alfa1 * sqrt(x(1)/self.A1) - self.alfa2 * (x(2)/self.C2)^(1/4);
			
			dx = [dx1, dx2];
		end
	end
end

