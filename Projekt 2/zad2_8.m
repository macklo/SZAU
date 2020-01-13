close all
clear
clc

tau = 4;

load("data\dane_ucz.mat")

Y = y(tau+2:end)';
M = [u(2:end-tau)' u(1:end-tau-1)' y(tau+1:end-1)' y(tau:end-2)'];
w = M\Y;
save("data/lin.mat", "w");

sim_length= size(u, 2);
y_m = zeros(size(u));
for k = tau+2:size(u, 2)
	q = [u(k-tau) u(k-tau-1) y_m(k-1) y_m(k-2)];
	y_m(k) = q*w;
end

figure
	subplot(2, 1, 1)
		hold on
		stairs(1:sim_length, y)
		stairs(1:sim_length, y_m)
		title("Przebiegi dla zbioru ucz¹cego")
		xlabel("k")
		ylabel("y")
		legend("Dane", "Model")
	subplot(2, 1, 2)
		hold on
		stairs(1:sim_length, u)
		xlabel("k")
		ylabel("u")
print(gcf,'./fig/2_8_przebieg_ucz', '-dmeta')

figure
	scatter(y, y_m, '.');
	title("Relacja dla zbioru ucz¹cego")
	xlabel("Dane")
	ylabel("Model")
print(gcf,'./fig/2_8_relacja_ucz', '-dmeta')


load("./data/dane_wer.mat")
sim_length= size(u, 2);
y_m = zeros(size(u));
for k = tau+2:size(u, 2)
	q = [u(k-tau) u(k-tau-1) y_m(k-1) y_m(k-2)];
	y_m(k) = q*w;
end

figure
	subplot(2, 1, 1)
		hold on
		stairs(1:sim_length, y)
		stairs(1:sim_length, y_m)
		title("Przebiegi dla zbioru weryfikuj¹cego")
		xlabel("k")
		ylabel("y")
		legend("Dane", "Model")
	subplot(2, 1, 2)
		hold on
		stairs(1:sim_length, u)
		xlabel("k")
		ylabel("u")
		
print(gcf,'./fig/2_8_przebieg_wer', '-dmeta')

figure
	scatter(y, y_m, '.');
	title("Relacja dla zbioru weryfikuj¹cego")
	xlabel("Dane")
	ylabel("Model")
print(gcf,'./fig/2_8_relacja_wer', '-dmeta')
	