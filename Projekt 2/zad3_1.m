close all
clear
clc

N = 200;
Nu = 200;
lambda = 100;


delta=1e-5;
tau=4;
sim_length = 200;

y_zad = build_random_setpoints_array(struct("y", 0), sim_length, 20, 20, 0, 8);

u = zeros(1, sim_length);
y_o = zeros(1, sim_length);
y_mlin = zeros(1, sim_length);
x1_o = zeros(1, sim_length);
x2_o = zeros(1, sim_length);
b4 = zeros(1, sim_length);
b5 = zeros(1, sim_length);
a1 = zeros(1, sim_length);
a2 = zeros(1, sim_length);


for k = tau+2:sim_length
  %1. Symulacja obiektu -> pomiar y(k)
  [x1_o(k), x2_o(k), y_o(k)] = object(u(k-4),x1_o(k-1),x2_o(k-1));
  
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
  
  %3. Oblicz odp. skokow¹
  sc=zeros(1, N);
  sc(4)=b4(k);
  for i= 5:N
      sc(N)=b4(k)+b5(k)-(a1(k)*sc(N-1)+a2(k)*sc(N-2));
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
  
  %8. Oblicz DELTAU
  Y_zad  = ones(N, 1)*y_zad(k);
  deltaU = K*(Y_zad - Y0);
  
  %9.Do sterowania pierwszy element
  u(k) = u(k-1) + deltaU(1);
  
  %10. Ewentualnie przytnij
% 	if(u(k) > 1)
% 		u(k) = 1;
% 	elseif (u(k) < -1)
% 		u(k) = -1;
%     end
    
end