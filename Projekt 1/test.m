close all
clear
clc

x0 = [1.150295897591316e+04, 1.101600093823581e+03];
u0 = 90;

y0 = 36;
% 
% workpoint = struct('x0', x0, 'u0', u0, 'y0', y0, 't0', 0);

t0 = 0;
tfinal = 4000;

legends = [];
y = [];
ystat = [];
ystatlin = []'

jumps = -1:0.2:1;
uJumps = u0 + jumps * u0;

figure
	hold on
	grid on

for mult = jumps
	mult
	u = u0 + mult * u0;
	
	tanks = TankSystem(x0);
	linearTanks = LinearTankSystem(x0);
	
	for t = t0:tanks.Ts:tfinal
		tanks.simulate(u);
		linearTanks.simulate(u);
	end
	y = sqrt(tanks.x(:, 2)/tanks.C2);
	ylin = linearTanks.x(:, 2) * linearTanks.C(1, 2) + y0;
	stairs(t0:tanks.Ts:tfinal, y)
	stairs(t0:tanks.Ts:tfinal, ylin, '--')
	ystat = [ystat y(end)];
	ystatlin = [ystatlin ylin(end)];
end

figure
	grid on
	hold on
	plot(uJumps, ystat, '-o');
	plot(uJumps, ystatlin, '-o');
	xlabel("F_1 [cm^3/s]");
	ylabel("h_2 [cm]");
	title("Charakterystyka statyczna");
