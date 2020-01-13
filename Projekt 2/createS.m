clear;
close all;
clc;

tau = 4;
load("data/lin.mat")

sim_length = 500;

s = zeros(1, sim_length);

start = 10;

u = 0*ones(1, sim_length);
y = 0*ones(1, sim_length);

u(start:end) = 1;
for k = tau+2:sim_length
	q = [u(k-tau) u(k-tau-1) y(k-1) y(k-2)];
	y(k) = q*w;
end
s = y(start+1:end) - y(start);

save('./data/s.mat', 's');
figure;
    stairs(s);