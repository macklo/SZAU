clc
close all
clear

alfa1 = -1.473409;
alfa2 = 0.525788;
beta1 = 0.026085;
beta2 = 0.021057;
noise = 0.01;
		
%% Dane ucz¹ce

tau = 4;
jumpInterval = 100;
u_rand = build_random_setpoints_array(struct("y", 0), 2000, jumpInterval, jumpInterval, -1, 1);
u = generate_u(0, jumpInterval, [0.2 0.4 -0.2 -1 1 0.8 0.7  -1 0]);
u = [u u_rand];
sim_length = size(u, 2);
stairs(u)
y = zeros(1, sim_length);
x1 = zeros(1, sim_length);
x2 = zeros(1, sim_length);

for k = (tau+1):sim_length
	g1 = (exp(6*u(k-4)) - 1)/(exp(6*u(k-4)) + 1);
	x1(k) = -alfa1*x1(k-1) + x2(k-1) + beta1*g1;
	x2(k) = -alfa2*x1(k-1) + beta2*g1;
	y(k) = -0.5*(1-exp(-2*x1(k))) - noise + 2*noise*rand;
end

figure
	subplot(2, 1, 1)
		stairs(1:sim_length, y)
		xlabel("k")
		ylabel("y")
		title("Dane ucz¹ce")
	subplot(2, 1, 2)
		stairs(1:sim_length, u)
		xlabel("k")
		ylabel("u")
% writematrix([u' y'], './dane.txt', 'Delimiter', 'space')
% save("data/dane_ucz.mat", "u", "y");

%% Dane weryfikuj¹ce
sim_length = 3000;
tau = 4;
jumpInterval = 100;
u_rand = build_random_setpoints_array(struct("y", 0), 2000, jumpInterval, jumpInterval, -1, 1);
u = generate_u(0, jumpInterval, [0.8 -0.2 -1 0.5 1 0.8 -1  0.6 0]);
u = [u u_rand];
y = zeros(1, sim_length);
x1 = zeros(1, sim_length);
x2 = zeros(1, sim_length);

for k = (tau+1):sim_length
	g1 = (exp(6*u(k-4)) - 1)/(exp(6*u(k-4)) + 1);
	x1(k) = -alfa1*x1(k-1) + x2(k-1) + beta1*g1;
	x2(k) = -alfa2*x1(k-1) + beta2*g1;
	y(k) = -0.5*(1-exp(-2*x1(k))) - noise + 2*noise*rand;
end

figure
	subplot(2, 1, 1)
		stairs(1:sim_length, y)
		xlabel("k")
		ylabel("y")
		title("Dane weryfikuj¹ce")
	subplot(2, 1, 2)
		stairs(1:sim_length, u)
		xlabel("k")
		ylabel("u")
% writematrix([u' y'], './dane_wer.txt', 'Delimiter', 'space')
% save("data/dane_wer.mat", "u", "y");