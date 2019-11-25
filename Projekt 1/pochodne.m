close all
clear
clc

syms A1 C2 alfa1 alfa2
syms F1 FD V1 V2

dV1 = F1 + FD - alfa1 * (V1/A1)^(1/2);
dV2 = alfa1 * (V1/A1)^(1/2) - alfa2 * (V2/C2)^(1/4);
y = (V2/C2)^(1/2);

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