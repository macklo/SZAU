classdef LinearTankSystem2 < AbstractObject
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
		function self = LinearTankSystem2(workpoint)
			ny = 1; nu = 1; nd = 0; Ts = 1;
			self@AbstractObject(ny, nu, nd, Ts);
			
			self.x0 = workpoint.x;
			self.u0 = workpoint.u;
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

			self.u = self.uk*ones(1, self.tau);
			self.y = self.yk*ones(1, self.tau);
			self.x = (self.xk'.*ones(2, self.tau))';
		end
		
		function [y, t] = simulate(self)
			x = self.x(end, :);
			u = self.u(self.tau);

			k1 = self.Ts * self.differential(x, u);
% 			k2 = self.Ts * self.differential(x + k1/2, u);
% 			k3 = self.Ts * self.differential(x + k2/2, u);
% 			k4 = self.Ts * self.differential(x + k3, u);

			x = x + k1;%(k1 + 2*k2 + 2*k3 + k4)/6;
			self.yk = (1/(2*self.C2*(self.x0(2)/self.C2)^(1/2))) * (x(2) - self.x0(2)) + self.y0;

			self.x = [self.x ; x];
		end
		
		
		function [t, x] = simulateODE(self, x0, u0, t0, tfinal)
			[t, x] = ode45(@(t, x) self.differential(x, u0), t0:self.Ts:tfinal, x0);
		end
		
		function dx = differential(self, x, u)
			dx1 = (u - self.u0) + (self.FD - self.FD0) - self.alfa1/(2*self.A1*(self.x0(1)/self.A1)^(1/2)) * (x(1) - self.x0(1));
            dx2 = (self.alfa1/(2*self.A1*(self.x0(1)/self.A1)^(1/2))) * (x(1) - self.x0(1)) - (self.alfa2/(4*self.C2*(self.x0(2)/self.C2)^(3/4))) * (x(2) - self.x0(2));
%             h1(2,t) = (x1(2,t))/A1;
%             h2(2,t) = h20 + 1/2*(C2*V20)^-0.5 * (x2(2,t) - V20);
% 			dx = self.A*x' + self.B*u;
			dx = [dx1 dx2];
		end
		
		function workpoint = calculateWorkpoint(self, u)
% 			u
			x1 = (((u - self.u0) + (self.FD - self.FD0)) / (self.alfa1/(2*self.A1*(self.x0(1)/self.A1)^(1/2)))) + self.x0(1);
			x2 = ((self.alfa1/(2*self.A1*(self.x0(1)/self.A1)^(1/2))) * (x1 - self.x0(1))) / ((self.alfa2/(4*self.C2*(self.x0(2)/self.C2)^(3/4)))) + self.x0(2);
			y = (1/(2*self.C2*(self.x0(2)/self.C2)^(1/2))) * (x2 - self.x0(2)) + self.y0;
			x = [x1 x2];
			workpoint = struct('x', x, 'u', u, 'y', y);
		end
	end
end

