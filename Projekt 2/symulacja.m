clc
close all
clear

alfa1 = -1.473409;
alfa2 = 0.525788;
beta1 = 0.026085;
beta2 = 0.021057;

sim_length = 4000;
tau = 4;
jumpInterval = 100;
u = build_random_setpoints_array(struct("y", 0), sim_length, jumpInterval, jumpInterval, -1, 1);
y = zeros(1, sim_length);
x1 = zeros(1, sim_length);
x2 = zeros(1, sim_length);

for k = (tau+1):sim_length
	g1 = (exp(6*u(k-4)) - 1)/(exp(6*u(k-4)) + 1);
	x1(k) = -alfa1*x1(k-1) + x2(k-1) + beta1*g1;
	x2(k) = -alfa2*x1(k-1) + beta2*g1;
	y(k) = -0.5*(1-exp(-2*x1(k)));
end

figure
	subplot(2, 1, 1)
		stairs(1:sim_length, y)
		xlabel("k")
		ylabel("y")
		title("y")
	subplot(2, 1, 2)
		stairs(1:sim_length, u)
		xlabel("k")
		ylabel("u")
		title("u")