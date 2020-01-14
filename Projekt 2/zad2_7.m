close all
clear
clc

tau = 4;

load("./data/2_6.mat")

best_ucz = min(E_ucz, [], 2);
best_wer = min(E_wer, [], 2);
[valI, i] = min(E_wer);

run("2_6_output/model_7_2.m")

load("./data/dane_ucz.mat")
sim_length= size(u, 2);
y_nn = zeros(size(u));
for k = tau+2:size(u, 2)
	q = [u(k-tau) u(k-tau-1) y_nn(k-1) y_nn(k-2)]';
	y_nn(k) = w20 + w2*tanh(w10 + w1*q);
end

figure
	subplot(2, 1, 1)
		hold on
		stairs(1:sim_length, y)
		stairs(1:sim_length, y_nn)
		title("Przebiegi dla zbioru ucz¹cego")
		xlabel("k")
		ylabel("y")
		legend("Dane", "Model")
	subplot(2, 1, 2)
		hold on
		stairs(1:sim_length, u)
		xlabel("k")
		ylabel("u")
print(gcf,'./fig/2_7_przebieg_ucz', '-dmeta')

figure
	scatter(y, y_nn, '.');
	title("Relacja dla zbioru ucz¹cego")
	xlabel("Dane")
	ylabel("Model")
print(gcf,'./fig/2_7_relacja_ucz', '-dmeta')


load("./data/dane_wer.mat")
sim_length= size(u, 2);
y_nn = zeros(size(u));
for k = tau+2:size(u, 2)
	q = [u(k-tau) u(k-tau-1) y_nn(k-1) y_nn(k-2)]';
	y_nn(k) = w20 + w2*tanh(w10 + w1*q);
end

figure
	subplot(2, 1, 1)
		hold on
		stairs(1:sim_length, y)
		stairs(1:sim_length, y_nn)
		title("Przebiegi dla zbioru weryfikuj¹cego")
		xlabel("k")
		ylabel("y")
		legend("Dane", "Model")
	subplot(2, 1, 2)
		hold on
		stairs(1:sim_length, u)
		xlabel("k")
		ylabel("u")
print(gcf,'./fig/2_7_przebieg_wer', '-dmeta')

figure
	scatter(y, y_nn, '.');
	title("Relacja dla zbioru weryfikuj¹cego")
	xlabel("Dane")
	ylabel("Model")
print(gcf,'./fig/2_7_relacja_wer', '-dmeta')
