classdef LinearTankSystem < AbstractObject
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
		
		V1 = 0;
		V2 = 0;
		
		x0 = []
		t = 0
		
		u, y, x
		uk, yk, xk
		
		A = [];
		B = [];
		C = [];
		D = [];
	end
	
	methods
		function self = LinearTankSystem(workpoint)
			ny = 1; nu = 1; nd = 0; Ts = 1;
			self@AbstractObject(ny, nu, nd, Ts);
			
			self.x0 = workpoint.x;
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
			self.uk = control - self.u0;
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
			self.uk = workPoint.u - workPoint.u;
			self.yk = workPoint.y;
			self.xk = workPoint.x - workPoint.x;

			self.u = self.uk*ones(1, self.tau);
			self.y = self.yk*ones(1, self.tau);
			self.x = (self.xk'.*ones(2, self.tau))';
		end
		
		function u = getU(self)
			time = size(self.u, 2) * self.Ts;
			if (time > self.tau)
				u = self.u(time - self.Ts);
			else
				u = self.u0;
			end
		end
		
		function [y, t] = simulate(self)
% 			self.u = [self.u, self.uk];
			if size(self.x, 1) == 0
				x = self.x0 - self.x0;
			else
				x = self.x(end, :);
				u = self.u(self.tau);
				
				k1 = self.Ts * self.differential(x, u);
				k2 = self.Ts * self.differential(x + k1/2, u);
				k3 = self.Ts * self.differential(x + k2/2, u);
				k4 = self.Ts * self.differential(x + k3, u);

				x = x + (k1 + 2*k2 + 2*k3 + k4)/6;
				self.yk = x(2)* self.C(1, 2) + self.y0;
			end
% 			if x(1) < 0
% 				x(1) = 0;
% 			end
% 			if x(2) < 0
% 				x(2) = 0;
% 			end
			self.x = [self.x ; x];
		end
		
		function xret = getX(self, x)
			xret = x;
			if x(1) < 0
				xret(1) = 0;
			end
			if x(2) < 0
				xret(2) = 0;
			end
		end
		
		function [t, x] = simulateODE(self, x0, u0, t0, tfinal)
			[t, x] = ode45(@(t, x) self.differential(x, u0), t0:self.Ts:tfinal, x0);
		end
		
		function dx = differential(self, x, u)
			dx = self.A*x' + self.B*u;
			dx = dx';
		end
	end
end

