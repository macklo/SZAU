% close all
clear
clc

N = 20;
Nu = 20;
lambda = 200;

delta=1e-5;
tau=4;
sim_length = 2000;

mode = 1; % 1-NPL, 2-GPC

load("./data/setpoints.mat")

u = zeros(1, sim_length);
y_o = zeros(1, sim_length);
y_mlin = zeros(1, sim_length);
y_nn = zeros(1, sim_length);
y_m = zeros(1, sim_length);
x1_o = zeros(1, sim_length);
x2_o = zeros(1, sim_length);
b4 = zeros(1, sim_length);
b5 = zeros(1, sim_length);
a1 = zeros(1, sim_length);
a2 = zeros(1, sim_length);

if mode == 2
	load("data/lin.mat")
	b = [0 0 0 w(1) w(2)];
	a = [-w(3) -w(4)];
	
	s=zeros(1, N);
	s(4)=b(4);
	for i= 5:N
		s(i)=b(4)+b(5)-(a(1)*s(i-1)+a(2)*s(i-2));
	end
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
end

for k = tau+2:sim_length
	%1. Symulacja obiektu -> pomiar y(k)
	[x1_o(k), x2_o(k), y_o(k)] = object(u(k-4),x1_o(k-1),x2_o(k-1));

	% NPL
	if mode == 1
		%2. Linearyzacja modelu
		b4(k) = (nn(u(k-4)+delta, u(k-5), y_o(k-1), y_o(k-2)) -... 
			  nn(u(k-4), u(k-5), y_o(k-1), y_o(k-2)))/delta;
		b5(k) = (nn(u(k-4), u(k-5)+delta, y_o(k-1), y_o(k-2)) -... 
			  nn(u(k-4), u(k-5), y_o(k-1), y_o(k-2)))/delta;
		a1(k) = -(nn(u(k-4), u(k-5), y_o(k-1)+delta, y_o(k-2)) -... 
			  nn(u(k-4), u(k-5), y_o(k-1), y_o(k-2)))/delta;
		a2(k) = -(nn(u(k-4), u(k-5), y_o(k-1), y_o(k-2)+delta) -... 
			  nn(u(k-4), u(k-5), y_o(k-1), y_o(k-2)))/delta;

		y_mlin(k) = b4(k)*u(k-4) +b5(k)*u(k-5) -a1(k)*y_o(k-1)-a2(k)*y_o(k-2);
		y_nn(k) = nn(u(k-4), u(k-5), y_o(k-1), y_o(k-2));

		%3. Oblicz odp. skokow¹
		sc=zeros(1, N);
		sc(4)=b4(k);
		for i= 5:N
			sc(i)=b4(k)+b5(k)-(a1(k)*sc(i-1)+a2(k)*sc(i-2));
		end

		%4. Oblicz macierz dynamiczna
		M = zeros(N, Nu);
		for i = 1:N
			for j = 1:Nu
				if (i >= j)
					M(i,j) = sc(i-j+1);
				else
					M(i,j) = 0;
				end
			end
		end

		%5. Oblicz K
		K = ((M'*M + lambda*eye(Nu, Nu))^(-1))*M';

		%6. Oblicz d
		d = y_o(k)- nn(u(k-4), u(k-5), y_o(k-1), y_o(k-2));

		%7. Oblicz trajektorie swobodna
		Y0 = zeros(N, 1);
		Y0(1)=nn(u(k-3), u(k-4), y_o(k), y_o(k-1)) +d;
		Y0(2)=nn(u(k-2), u(k-3), Y0(1), y_o(k)) +d;
		Y0(3)=nn(u(k-1), u(k-2), Y0(2), Y0(1)) +d;
		for i =4:N
		Y0(i)=nn(u(k-1), u(k-1), Y0(i-1), Y0(i-2)) +d; 
		end
	end
	
	%GPC
	if mode == 2
		q = [u(k-tau) u(k-tau-1) y_o(k-1) y_o(k-2)];
		y_m(k) = q*w;
		dk = y_o(k) - y_m(k);

		Y0 = zeros(N, 1);
		Y0(1) = b(4) * u(k-3) + b(5) * u(k-4) - a(1)*y_o(k) - a(2)*y_o(k-1) + dk;
		Y0(2) = b(4) * u(k-2) + b(5) * u(k-3) - a(1)*Y0(1)  - a(2)*y_o(k)   + dk;
		Y0(3) = b(4) * u(k-1) + b(5) * u(k-2) - a(1)*Y0(2)  - a(2)*Y0(1)    + dk;
		Y0(4) = b(4) * u(k-1) + b(5) * u(k-1) - a(1)*Y0(3)  - a(2)*Y0(2)    + dk;
		for i = 5:N
			Y0(i) = b(4) * u(k-1) + b(5) * u(k-1) - a(1)*Y0(i-1) - a(2)*Y0(i-2) + dk;
		end
	end
	
	%8. Oblicz DELTAU
	Y_zad  = ones(N, 1)*y_zad(k);
	deltaU = K*(Y_zad - Y0);

	%9.Do sterowania pierwszy element
	u(k) = u(k-1) + deltaU(1);

	%10. Ewentualnie przytnij
	if(u(k) > 1)
		u(k) = 1;
	elseif (u(k) < -1)
		u(k) = -1;
	end
end

figure
	subplot(2, 1, 1)
		hold on
		stairs(1:sim_length, y_o)
		stairs(1:sim_length, y_zad)
		title("Wartoœæ wyjœciowa")
		ylabel("y")
		xlabel("k")
        legend("Wyjœcie obiektu","Wartoœæ zadana")
		
	subplot(2, 1, 2)
		stairs(1:sim_length, u)
		title("Sterowanie")
		ylabel("u")
		xlabel("k")
print(gcf,"./fig/zad3NPL" + " N= " + num2str(N) + " Nu= " + num2str(Nu) + " lambda= " + num2str(lambda) , '-dmeta')