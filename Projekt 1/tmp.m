close all
clear
clc

x0 = [11502; 385.7950];
% x0 = [1 1 1 1];
u0 = 90;

y0 = 21.31;

workpoint = struct('x0', x0, 'u0', u0, 'y0', y0, 't0', 0);

t0 = 0;
tfinal = 10000;

r = Reactor();

figure
hold on
grid on
legends = [];
y = [];

for t = t0:r.Ts:tfinal
	r.simulate(90);
end

y = sqrt(r.x(:, 2)/r.C2);
stairs(y)

% 
% for mult = -0.5:0.1:1
% 	u = u0 + mult * u0;
% 	legends = [legends u];
% 	[t, x] = r.simulateODE(x0, u, t0, tfinal);
% 	y = sqrt(x(:, 2)/r.C2);
% 	stairs(t, y)
% end