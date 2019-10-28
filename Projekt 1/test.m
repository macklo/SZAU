close all
clear
clc

x0 = [5.787739978769097e+03; 0.022775759372083];
% x0 = [1 1 1 1];
u0 = 90;

y0 = 21.31;

workpoint = struct('x0', x0, 'u0', u0, 'y0', y0, 't0', 0);

t0 = 0;
tfinal = 100;

r = Reactor();

figure
hold on
grid on
legends = [];
y = [];

for mult = 0 %-0.8:0.1:1
	u = u0 + mult * u0;
	legends = [legends u];
	[t, x] = r.simulateODE(x0, u, t0, tfinal);
	y = sqrt(x(:, 2)/r.C2);
	stairs(t, y)
	
	react = Reactor();
	
	for t = t0:react.Ts:tfinal
		react.simulate(u);
	end
	y = sqrt(react.x(:, 2)/react.C2);
	stairs(t0:react.Ts:tfinal, y)
end



