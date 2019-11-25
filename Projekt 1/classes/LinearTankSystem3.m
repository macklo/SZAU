classdef LinearTankSystem3 < AbstractObject
	properties
		A1    = 540;
		C2    = 0.85;
		alfa1 = 26;
		alfa2 = 20;
		
		F1  = 90;
		FD  = 30;
		FD0 = 30;
		tau = 100;
		h2  = 36;
		
		u0;
		y0;
		h10;
		h20;
		
		V1 = 0;
		V2 = 0;
		
		x0 = []
		t = 0
		
		u, y, x
		uk, yk, xk, h1k, h2k
		
		A = [];
		B = [];
		C = [];
		D = [];
	end
	
	methods
		function self = LinearTankSystem3(workpoint)
			ny = 1; nu = 1; nd = 0; Ts = 1;
			self@AbstractObject(ny, nu, nd, Ts);
			
			self.x0 = workpoint.x;
			self.u0 = workpoint.u;
			self.h10 = workpoint.h1;
			self.h20 = workpoint.h2;
			self.y0 = workpoint.y;
			self.V1 = self.x0(1);
			self.V2 = self.x0(2);
			
			self.A = [
				-self.alfa1/(2*self.A1*(self.V1/self.A1)^(1/2)), 0;
				self.alfa1/(2*self.A1*(self.V1/self.A1)^(1/2)), -self.alfa2/(4*self.C2*(self.V2/self.C2)^(3/4))
				];

			self.B = [
				1;
				0
				];

			self.C = [
				0, 1/(2*self.C2*(self.V2/self.C2)^(1/2))
				];

			self.D = 0;
		end
		
		function output = getOutput(self)
			output = self.yk;
		end
		
		function setControl(self, control)
			self.uk = control;
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
			workpoint = self.calculateWorkpoint(workPoint.u);
			workPoint = workpoint;
			self.uk = workPoint.u;
			self.yk = workPoint.y;
			self.xk = workPoint.x;
			self.h1k = workPoint.h1;
			self.h2k = workPoint.h2;
			

			self.u = self.uk*ones(1, self.tau);
			self.y = self.yk*ones(1, self.tau);
			self.x = (self.xk'.*ones(2, self.tau))';
		end
		
		function [y, t] = simulate(self)
			x = self.x(end, :);
			u = self.u(self.tau);

			k1 = self.Ts * self.differential(u);
% 			k2 = self.Ts * self.differential(x + k1/2, u);
% 			k3 = self.Ts * self.differential(x + k2/2, u);
% 			k4 = self.Ts * self.differential(x + k3, u);

			x = x + k1;%(k1 + 2*k2 + 2*k3 + k4)/6;
			
			h1 = self.h10 + (x(1) - self.x0(1))/self.A1;
			h2 = self.h20 + 1/2*(self.C2*self.x0(2))^-0.5 * (x(2) - self.x0(2));
			
			self.h1k = h1;
			self.h2k = h2;
			self.yk = h2;

			self.x = [self.x ; x];
		end
		
		
		function [t, x] = simulateODE(self, x0, u0, t0, tfinal)
			[t, x] = ode45(@(t, x) self.differential(x, u0), t0:self.Ts:tfinal, x0);
		end
		
		function dx = differential(self, u)
			dV1 = (u - self.u0) + (self.FD - self.FD0) - self.alfa1/2*self.h10^-0.5 * (self.h1k - self.h10);
			dV2 = self.alfa1/2*self.h10^-0.5 * (self.h1k - self.h10) - self.alfa2/2*self.h20^-0.5 * (self.h2k - self.h20);
			dx = [dV1 dV2];
		end
		
		function workpoint = calculateWorkpoint(self, u)
			h1 = (((u - self.u0) + (self.FD - self.FD0)) / ((self.alfa1/2*self.h10^-0.5))) + self.h10;
			h2 = (self.alfa1/2*self.h10^-0.5 * (h1 - self.h10) / (self.alfa2/2*self.h20^-0.5)) + self.h20;
			
			x1 = self.A1*(h1 - self.h10) + self.x0(1);
			x2 = ((h2 - self.h20) / (1/2*(self.C2*self.x0(2))^-0.5)) + self.x0(2);
			
			x = [x1 x2];
			workpoint = struct('x', x, 'u', u, 'y', h2, 'h1', h1, 'h2', h2);
		end
	end
end

