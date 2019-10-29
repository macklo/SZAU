classdef LinearTankSystem < handle
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
		
		Ts = 1;
		
		x0 = []
		t = 0
		x = [];
		u = [];
		
		A = [];
		B = [];
		C = [];
		D = [];
	end
	
	methods
		function obj = LinearTankSystem(x0)
			obj.x0 = x0;
			obj.V1 = x0(1);
			obj.V2 = x0(2);
			
			obj.A = [
				-obj.alfa1/(2*obj.A1*(obj.V1/obj.A1)^(1/2)), 0;
				obj.alfa1/(2*obj.A1*(obj.V1/obj.A1)^(1/2)), -obj.alfa2/(4*obj.C2*(obj.V2/obj.C2)^(3/4))
				];

			obj.B = [
				1;
				0
				];

			obj.C = [
				0, 1/(2*obj.C2*(obj.V2/obj.C2)^(1/2))
				];

			obj.D = 0;
		end
		
		function u = getU(obj)
			time = size(obj.u, 2) * obj.Ts;
			if (time > obj.tau)
				u = obj.u(time - obj.Ts);
			else
				u = obj.u0;
			end
		end
		
		function [y, t] = simulate(obj, u)
			obj.u = [obj.u, u];
			if size(obj.x, 1) == 0
				x = obj.x0 - obj.x0;
			else
				x = obj.x(end, :);
				u = obj.getU() - obj.u0;
				
				k1 = obj.Ts * obj.differential(x, u);
				k2 = obj.Ts * obj.differential(x + k1/2, u);
				k3 = obj.Ts * obj.differential(x + k2/2, u);
				k4 = obj.Ts * obj.differential(x + k3, u);

				x = x + (k1 + 2*k2 + 2*k3 + k4)/6;
			end
% 			if x(1) < 0
% 				x(1) = 0;
% 			end
% 			if x(2) < 0
% 				x(2) = 0;
% 			end
			obj.x = [obj.x ; x];
		end
		
		function xret = getX(obj, x)
			xret = x;
			if x(1) < 0
				xret(1) = 0;
			end
			if x(2) < 0
				xret(2) = 0;
			end
		end
		
		function [t, x] = simulateODE(obj, x0, u0, t0, tfinal)
			[t, x] = ode45(@(t, x) obj.differential(x, u0), t0:obj.Ts:tfinal, x0);
		end
		
		function dx = differential(obj, x, u)
			dx = obj.A*x' + obj.B*u;
			dx = dx';
		end
	end
end

