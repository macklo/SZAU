close all
clear
clc

tau = 4;

load("./data/2_6.mat")

best_ucz = min(E_ucz, [], 2);
best_wer = min(E_wer, [], 2);
[valI, i] = min(E_wer);

w10(1,1)=-4.6198703823e-001; w1(1,1)=-3.5029875010e-001; w1(1,2)=-3.4199979359e-001; w1(1,3)=-2.2730766724e-003; w1(1,4)=-3.0706465060e-001; 
w10(2,1)=3.1820463853e-001; w1(2,1)=6.9097494306e-001; w1(2,2)=5.8780683327e-001; w1(2,3)=2.5070441985e-001; w1(2,4)=-9.7260211362e-001; 
w10(3,1)=1.9609413689e-001; w1(3,1)=6.6344282038e-001; w1(3,2)=6.7813363857e-001; w1(3,3)=-7.2342695999e-001; w1(3,4)=5.7175264093e-001; 
w10(4,1)=-1.0237351981e+000; w1(4,1)=5.8692840374e-001; w1(4,2)=6.1197512936e-001; w1(4,3)=-1.0922788685e+000; w1(4,4)=-7.5028153135e-001; 
w10(5,1)=7.6064771786e-001; w1(5,1)=-3.5319352731e-001; w1(5,2)=8.2915842497e-002; w1(5,3)=1.2669601768e+000; w1(5,4)=1.0523933415e+000; 
w10(6,1)=3.5274572645e-001; w1(6,1)=1.4220051502e-001; w1(6,2)=3.5484818020e-001; w1(6,3)=-6.6601109117e-001; w1(6,4)=-2.4651654314e-001; 
w10(7,1)=-1.0375293419e+000; w1(7,1)=-4.1205423860e-002; w1(7,2)=-2.6542124773e-002; w1(7,3)=5.1262090742e-001; w1(7,4)=-3.0077229963e-002; 
w20=1.1664227150e+000; w2(1)=-4.2618664672e-001; w2(2)=4.3193305647e-001; w2(3)=-4.7991956335e-001; w2(4)=1.5547681868e-001; w2(5)=2.2844205557e-001; w2(6)=-6.9613195381e-001; w2(7)=1.5186459277e+000; 

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
