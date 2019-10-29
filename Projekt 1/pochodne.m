close all
clear
clc

% syms F_in F_c F V k_0 E R 
% syms ro ro_c c_p c_pc T_in T_Cin T h a b E_R
% syms C_A T C_Ain F_c Ca0 T0

syms A1 C2 alfa1 alfa2
syms F1 FD V1 V2

dV1 = F1 + FD - alfa1 * (V1/A1)^(1/2);
dV2 = alfa1 * (V1/A1)^(1/2) - alfa2 * (V2/C2)^(1/4);
y = (V2/C2)^(1/2);

% dCadt=(F_in*C_Ain -F*C_A - V*k_0*exp(-E_R/T)*C_A)/V %(V przerzucamy na drug¹ stronê
% dTdt= (F_in*ro*c_p*T_in - F*ro*c_p*T +...
%     V*h*k_0*exp(-E_R/T)*C_A -((F_c^(b + 1)*a*(T-T_Cin))/(F_c + ((F_c^b*a)/(2*c_pc*ro)))))/(V*c_p*ro) %%V*c_p*ro przerzucamy na drug¹ stronê


%Wyliczamy pochodne cz¹stkowe
disp('dV1 po F1');
diff(dV1,F1) 
disp('dV1 po FD');
diff(dV1,FD) 
disp('dV1 po V1');
diff(dV1,V1)
disp('dV1 po V2');
diff(dV1,V2) 

%---------------------------------
disp('dV2 po F1');
diff(dV2,F1)
disp('dV2 po FD');
diff(dV2,FD)
disp('dV2 po V1');
diff(dV2,V1)
disp('dV2 po V2');
diff(dV2,V2)

%%
disp('y po F1');
diff(y,F1)
disp('y po FD');
diff(y,FD)
disp('y po V1');
diff(y,V1)
disp('y po V2');
diff(y,V2)