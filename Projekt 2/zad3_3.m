close all
clear
clc

tau = 4;
alfa1 = -1.473409;
alfa2 = 0.525788;
beta1 = 0.026085;
beta2 = 0.021057;

load("./data/liczba_neuronow_wyniki.mat")

best_ucz = min(E_ucz, [], 2);
best_wer = min(E_wer, [], 2);

run("sieci_output/model_6_4.m")
sim_length = 2000;

y_zad = build_random_setpoints_array(struct("y", 0), sim_length, 100, 100, -0.4, 2.4);

u = zeros(1, sim_length);
y = zeros(1, sim_length);
y_m = zeros(1, sim_length);
x1 = zeros(1, sim_length);
x2 = zeros(1, sim_length);

load("data/lin.mat")

N = 20;
Nu = 20;
lambda = 250;

na = 2;
nb = 5;
b = [0 0 0 w(1) w(2)];
a = [-w(3) -w(4)];


b4=b(4);
b5=b(5);
a1=a(1);
a2=a(2);

s = zeros(1, N);
for j = 1:N
	sum1 = 0;
	sum2 = 0;
	for i = 1:min([j, size(b, 2)])
		sum1 = sum1+b(i);
	end
	for i = 1:min([j-1, size(a, 2)])
		sum2 = sum2+a(i)*s(j-i);
	end
	s(j) = sum1 - sum2;
end
figure
stairs(s)

M = zeros(N, Nu);
for i = 1:N
	for j = 1:Nu
		if (i >= j)
			M(i,j) = s(i-j+1);
		else
			M(i,j) = 0;
		end
	end
end

K =(M'*M + lambda*eye(Nu, Nu))\M';

for k = tau+2:sim_length
	g1 = (exp(6*u(k-4)) - 1)/(exp(6*u(k-4)) + 1);
	x1(k) = -alfa1*x1(k-1) + x2(k-1) + beta1*g1;
	x2(k) = -alfa2*x1(k-1) + beta2*g1;
	y(k) = -0.5*(1-exp(-2*x1(k)));
	
	q = [u(k-tau) u(k-tau-1) y(k-1) y(k-2)];
	y_m(k) = q*w;
	dk = y(k) - y_m(k);
	
	Y0 = zeros(N, 1);
	Y0(1) = b(4) * u(k-3) + b(5) * u(k-4) - a(1)*y(k)  - a(2)*y(k-1) + dk;
	Y0(2) = b(4) * u(k-2) + b(5) * u(k-3) - a(1)*Y0(1) - a(2)*y(k)   + dk;
	Y0(3) = b(4) * u(k-1) + b(5) * u(k-2) - a(1)*Y0(2) - a(2)*Y0(1)  + dk;
	Y0(4) = b(4) * u(k-1) + b(5) * u(k-1) - a(1)*Y0(3) - a(2)*Y0(2)  + dk;
	for i = 5:N
		Y0(i) = b(4) * u(k-1) + b(5) * u(k-1) - a(1)*Y0(i-1) - a(2)*Y0(i-2) + dk;
	end
	Y_zad  = ones(N, 1)*y_zad(k);
	deltaU = K*(Y_zad - Y0);
	u(k) = u(k-1) + deltaU(1);
	if(u(k) > 1)
		u(k) = 1;
	elseif (u(k) < -1)
		u(k) = -1;
	end
end

figure
	subplot(2, 1, 1)
		hold on
		stairs(1:sim_length, y)
		stairs(1:sim_length, y_zad)
		title("Przebieg")
		
	subplot(2, 1, 2)
		hold on
		stairs(1:sim_length, u)
