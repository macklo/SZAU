%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.250370e+003; foe(n+1)=4.575397e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.483935e+003; foe(n+1)=2.743396e+003; krok(n+1)=4.526501e-004; ng(n+1)=6.001902e+003;
n=2; farx(n+1)=9.850708e+002; foe(n+1)=1.244963e+003; krok(n+1)=1.420248e-003; ng(n+1)=3.814401e+003;
n=3; farx(n+1)=1.259533e+003; foe(n+1)=7.865186e+002; krok(n+1)=1.476158e-004; ng(n+1)=1.179484e+004;
n=4; farx(n+1)=1.973186e+003; foe(n+1)=5.913651e+002; krok(n+1)=6.190058e-004; ng(n+1)=1.064958e+004;
n=5; farx(n+1)=1.851670e+003; foe(n+1)=5.461380e+002; krok(n+1)=4.625835e-003; ng(n+1)=1.128203e+003;
n=6; farx(n+1)=1.685603e+003; foe(n+1)=5.238628e+002; krok(n+1)=1.087346e-003; ng(n+1)=1.452306e+003;
n=7; farx(n+1)=1.584431e+003; foe(n+1)=5.184651e+002; krok(n+1)=3.970634e-004; ng(n+1)=1.272666e+003;
n=8; farx(n+1)=9.216613e+002; foe(n+1)=4.527827e+002; krok(n+1)=7.088910e-003; ng(n+1)=2.868791e+003;
n=9; farx(n+1)=5.743871e+002; foe(n+1)=4.181717e+002; krok(n+1)=6.295194e-004; ng(n+1)=2.781675e+003;
n=10; farx(n+1)=3.698502e+002; foe(n+1)=3.925899e+002; krok(n+1)=6.694872e-004; ng(n+1)=2.922986e+003;
n=11; farx(n+1)=2.816364e+002; foe(n+1)=3.734129e+002; krok(n+1)=1.040544e-003; ng(n+1)=3.978570e+003;
n=12; farx(n+1)=2.081822e+002; foe(n+1)=3.515110e+002; krok(n+1)=1.441974e-003; ng(n+1)=4.073792e+003;
n=13; farx(n+1)=1.671199e+002; foe(n+1)=3.301651e+002; krok(n+1)=8.368557e-003; ng(n+1)=3.617439e+003;
n=14; farx(n+1)=1.664734e+002; foe(n+1)=2.804022e+002; krok(n+1)=3.427028e-003; ng(n+1)=1.696404e+003;
n=15; farx(n+1)=1.673078e+002; foe(n+1)=2.727457e+002; krok(n+1)=4.966995e-004; ng(n+1)=1.323265e+003;
n=16; farx(n+1)=1.419975e+002; foe(n+1)=2.626636e+002; krok(n+1)=5.980899e-003; ng(n+1)=1.268557e+003;
n=17; farx(n+1)=9.356876e+001; foe(n+1)=2.247270e+002; krok(n+1)=1.505777e-003; ng(n+1)=1.448185e+003;
n=18; farx(n+1)=8.284185e+001; foe(n+1)=2.182244e+002; krok(n+1)=5.465324e-004; ng(n+1)=1.015000e+003;
n=19; farx(n+1)=8.444525e+001; foe(n+1)=2.146895e+002; krok(n+1)=9.481985e-004; ng(n+1)=1.896569e+003;
n=20; farx(n+1)=9.075969e+001; foe(n+1)=2.012341e+002; krok(n+1)=5.711090e-003; ng(n+1)=1.302567e+003;
n=21; farx(n+1)=9.356973e+001; foe(n+1)=1.934220e+002; krok(n+1)=3.973217e-003; ng(n+1)=1.140897e+003;
n=22; farx(n+1)=9.467109e+001; foe(n+1)=1.744071e+002; krok(n+1)=2.660052e-003; ng(n+1)=2.733795e+003;
n=23; farx(n+1)=1.138175e+002; foe(n+1)=1.520441e+002; krok(n+1)=9.053002e-004; ng(n+1)=2.089589e+003;
n=24; farx(n+1)=1.179893e+002; foe(n+1)=1.505729e+002; krok(n+1)=1.955288e-004; ng(n+1)=3.482445e+003;
n=25; farx(n+1)=1.236608e+002; foe(n+1)=1.489099e+002; krok(n+1)=1.284363e-003; ng(n+1)=3.326155e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=1.037942e+002; foe(n+1)=1.262714e+002; krok(n+1)=3.502676e-005; ng(n+1)=3.158162e+003;
n=27; farx(n+1)=6.138487e+001; foe(n+1)=8.656642e+001; krok(n+1)=4.461789e-005; ng(n+1)=4.416948e+003;
n=28; farx(n+1)=4.237537e+001; foe(n+1)=7.296479e+001; krok(n+1)=2.266194e-004; ng(n+1)=1.446813e+003;
n=29; farx(n+1)=3.684173e+001; foe(n+1)=5.299673e+001; krok(n+1)=5.098909e-005; ng(n+1)=3.085390e+003;
n=30; farx(n+1)=2.913738e+001; foe(n+1)=4.302318e+001; krok(n+1)=1.161194e-003; ng(n+1)=2.059076e+003;
n=31; farx(n+1)=2.550154e+001; foe(n+1)=3.795963e+001; krok(n+1)=1.293127e-003; ng(n+1)=1.391039e+003;
n=32; farx(n+1)=1.982557e+001; foe(n+1)=3.339828e+001; krok(n+1)=1.505777e-003; ng(n+1)=1.664181e+003;
n=33; farx(n+1)=1.677169e+001; foe(n+1)=2.909699e+001; krok(n+1)=2.043687e-003; ng(n+1)=2.582224e+003;
n=34; farx(n+1)=1.387983e+001; foe(n+1)=2.520607e+001; krok(n+1)=5.478088e-003; ng(n+1)=9.650274e+002;
n=35; farx(n+1)=1.120833e+001; foe(n+1)=2.263860e+001; krok(n+1)=4.323984e-003; ng(n+1)=1.258991e+003;
n=36; farx(n+1)=8.577846e+000; foe(n+1)=1.920577e+001; krok(n+1)=4.047652e-003; ng(n+1)=4.985344e+002;
n=37; farx(n+1)=8.357068e+000; foe(n+1)=1.773051e+001; krok(n+1)=2.349366e-003; ng(n+1)=7.802424e+002;
n=38; farx(n+1)=7.570653e+000; foe(n+1)=1.687598e+001; krok(n+1)=3.965892e-003; ng(n+1)=3.804657e+002;
n=39; farx(n+1)=6.985262e+000; foe(n+1)=1.619229e+001; krok(n+1)=3.721397e-003; ng(n+1)=7.234492e+002;
n=40; farx(n+1)=6.542448e+000; foe(n+1)=1.518954e+001; krok(n+1)=5.829556e-003; ng(n+1)=4.306223e+002;
n=41; farx(n+1)=6.255664e+000; foe(n+1)=1.403090e+001; krok(n+1)=5.107684e-003; ng(n+1)=3.476820e+002;
n=42; farx(n+1)=5.907409e+000; foe(n+1)=1.298975e+001; krok(n+1)=7.191676e-003; ng(n+1)=7.642274e+002;
n=43; farx(n+1)=5.470647e+000; foe(n+1)=1.253939e+001; krok(n+1)=5.242823e-003; ng(n+1)=7.226554e+002;
n=44; farx(n+1)=5.070353e+000; foe(n+1)=1.182498e+001; krok(n+1)=3.376031e-003; ng(n+1)=1.173676e+003;
n=45; farx(n+1)=4.833390e+000; foe(n+1)=1.058445e+001; krok(n+1)=1.715934e-002; ng(n+1)=7.904877e+002;
n=46; farx(n+1)=4.243908e+000; foe(n+1)=9.477002e+000; krok(n+1)=2.440235e-002; ng(n+1)=6.849680e+002;
n=47; farx(n+1)=4.167393e+000; foe(n+1)=9.315241e+000; krok(n+1)=4.473296e-003; ng(n+1)=2.702138e+002;
n=48; farx(n+1)=4.165903e+000; foe(n+1)=9.028375e+000; krok(n+1)=2.191235e-002; ng(n+1)=2.187094e+002;
n=49; farx(n+1)=4.202090e+000; foe(n+1)=8.918162e+000; krok(n+1)=8.580693e-003; ng(n+1)=1.531942e+002;
n=50; farx(n+1)=4.128849e+000; foe(n+1)=8.746960e+000; krok(n+1)=3.497807e-002; ng(n+1)=2.572809e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=4.138260e+000; foe(n+1)=8.725278e+000; krok(n+1)=1.595012e-005; ng(n+1)=1.671830e+002;
n=52; farx(n+1)=4.136130e+000; foe(n+1)=8.714789e+000; krok(n+1)=2.733737e-005; ng(n+1)=9.313175e+001;
n=53; farx(n+1)=4.105584e+000; foe(n+1)=8.671640e+000; krok(n+1)=2.431917e-004; ng(n+1)=6.478989e+001;
n=54; farx(n+1)=4.125400e+000; foe(n+1)=8.627255e+000; krok(n+1)=3.355212e-004; ng(n+1)=5.929872e+001;
n=55; farx(n+1)=4.239935e+000; foe(n+1)=8.347258e+000; krok(n+1)=2.419456e-003; ng(n+1)=5.850153e+001;
n=56; farx(n+1)=4.006266e+000; foe(n+1)=7.987062e+000; krok(n+1)=1.127471e-002; ng(n+1)=4.648731e+001;
n=57; farx(n+1)=3.820199e+000; foe(n+1)=7.835809e+000; krok(n+1)=2.929274e-003; ng(n+1)=1.957220e+002;
n=58; farx(n+1)=3.713644e+000; foe(n+1)=7.630619e+000; krok(n+1)=8.378898e-003; ng(n+1)=3.360347e+002;
n=59; farx(n+1)=3.790435e+000; foe(n+1)=7.415241e+000; krok(n+1)=7.680234e-003; ng(n+1)=4.227521e+002;
n=60; farx(n+1)=3.975446e+000; foe(n+1)=7.102720e+000; krok(n+1)=1.631060e-002; ng(n+1)=3.982643e+002;
n=61; farx(n+1)=3.931471e+000; foe(n+1)=7.036104e+000; krok(n+1)=4.874671e-003; ng(n+1)=1.774307e+002;
n=62; farx(n+1)=4.026355e+000; foe(n+1)=6.871652e+000; krok(n+1)=9.659422e-003; ng(n+1)=3.180495e+002;
n=63; farx(n+1)=4.130584e+000; foe(n+1)=6.434560e+000; krok(n+1)=5.335388e-002; ng(n+1)=2.031410e+002;
n=64; farx(n+1)=4.169369e+000; foe(n+1)=6.279880e+000; krok(n+1)=1.032555e-002; ng(n+1)=2.808796e+002;
n=65; farx(n+1)=4.283132e+000; foe(n+1)=6.145284e+000; krok(n+1)=2.426108e-002; ng(n+1)=1.229837e+002;
n=66; farx(n+1)=4.370452e+000; foe(n+1)=5.974294e+000; krok(n+1)=1.945154e-002; ng(n+1)=9.067692e+001;
n=67; farx(n+1)=4.198452e+000; foe(n+1)=5.837328e+000; krok(n+1)=1.725161e-002; ng(n+1)=8.335480e+001;
n=68; farx(n+1)=4.163924e+000; foe(n+1)=5.793086e+000; krok(n+1)=4.872032e-003; ng(n+1)=2.014420e+002;
n=69; farx(n+1)=4.059975e+000; foe(n+1)=5.539770e+000; krok(n+1)=5.962927e-002; ng(n+1)=1.031590e+002;
n=70; farx(n+1)=4.056689e+000; foe(n+1)=5.276823e+000; krok(n+1)=2.617232e-002; ng(n+1)=1.467907e+002;
n=71; farx(n+1)=4.023043e+000; foe(n+1)=5.138860e+000; krok(n+1)=3.285344e-002; ng(n+1)=1.792657e+002;
n=72; farx(n+1)=3.754543e+000; foe(n+1)=4.969777e+000; krok(n+1)=2.835564e-002; ng(n+1)=2.988322e+002;
n=73; farx(n+1)=3.673182e+000; foe(n+1)=4.856491e+000; krok(n+1)=7.242402e-003; ng(n+1)=2.179206e+002;
n=74; farx(n+1)=3.494981e+000; foe(n+1)=4.768860e+000; krok(n+1)=1.761298e-002; ng(n+1)=2.857780e+002;
n=75; farx(n+1)=3.012825e+000; foe(n+1)=4.429328e+000; krok(n+1)=6.736795e-002; ng(n+1)=8.247511e+001;
%odnowa zmiennej metryki
n=76; farx(n+1)=2.999726e+000; foe(n+1)=4.386485e+000; krok(n+1)=1.476242e-005; ng(n+1)=2.312379e+002;
n=77; farx(n+1)=2.988313e+000; foe(n+1)=4.377485e+000; krok(n+1)=1.663763e-005; ng(n+1)=9.903109e+001;
n=78; farx(n+1)=2.947929e+000; foe(n+1)=4.332041e+000; krok(n+1)=2.114007e-004; ng(n+1)=6.969002e+001;
n=79; farx(n+1)=2.887689e+000; foe(n+1)=4.257592e+000; krok(n+1)=4.175751e-004; ng(n+1)=5.837171e+001;
n=80; farx(n+1)=2.890530e+000; foe(n+1)=4.243596e+000; krok(n+1)=4.029366e-004; ng(n+1)=2.965652e+001;
n=81; farx(n+1)=2.944344e+000; foe(n+1)=4.192909e+000; krok(n+1)=4.534550e-003; ng(n+1)=2.129354e+001;
n=82; farx(n+1)=2.926614e+000; foe(n+1)=4.073184e+000; krok(n+1)=6.607631e-003; ng(n+1)=6.636907e+001;
n=83; farx(n+1)=2.913168e+000; foe(n+1)=4.007435e+000; krok(n+1)=5.000456e-003; ng(n+1)=2.528200e+002;
n=84; farx(n+1)=2.787380e+000; foe(n+1)=3.905737e+000; krok(n+1)=1.761298e-002; ng(n+1)=2.847580e+002;
n=85; farx(n+1)=2.787861e+000; foe(n+1)=3.834628e+000; krok(n+1)=9.630946e-003; ng(n+1)=2.426483e+002;
n=86; farx(n+1)=2.765702e+000; foe(n+1)=3.741380e+000; krok(n+1)=2.432476e-002; ng(n+1)=2.369230e+002;
n=87; farx(n+1)=2.623807e+000; foe(n+1)=3.628416e+000; krok(n+1)=1.640862e-002; ng(n+1)=2.754253e+002;
n=88; farx(n+1)=2.351148e+000; foe(n+1)=3.525687e+000; krok(n+1)=2.753897e-002; ng(n+1)=7.653651e+001;
n=89; farx(n+1)=2.226328e+000; foe(n+1)=3.455669e+000; krok(n+1)=2.251094e-002; ng(n+1)=1.651209e+002;
n=90; farx(n+1)=2.061602e+000; foe(n+1)=3.362450e+000; krok(n+1)=1.835014e-002; ng(n+1)=1.106143e+002;
n=91; farx(n+1)=1.966825e+000; foe(n+1)=3.259212e+000; krok(n+1)=2.421456e-002; ng(n+1)=9.638155e+001;
n=92; farx(n+1)=1.880000e+000; foe(n+1)=3.165032e+000; krok(n+1)=1.111520e-002; ng(n+1)=2.198994e+002;
n=93; farx(n+1)=1.794208e+000; foe(n+1)=3.090211e+000; krok(n+1)=1.564303e-002; ng(n+1)=1.614265e+002;
n=94; farx(n+1)=1.725470e+000; foe(n+1)=3.016356e+000; krok(n+1)=1.633292e-002; ng(n+1)=3.251107e+002;
n=95; farx(n+1)=1.697735e+000; foe(n+1)=2.954493e+000; krok(n+1)=7.801084e-003; ng(n+1)=4.160083e+002;
n=96; farx(n+1)=1.679545e+000; foe(n+1)=2.879486e+000; krok(n+1)=3.743558e-002; ng(n+1)=3.861292e+002;
n=97; farx(n+1)=1.642369e+000; foe(n+1)=2.699961e+000; krok(n+1)=4.799811e-002; ng(n+1)=3.993621e+002;
n=98; farx(n+1)=1.607143e+000; foe(n+1)=2.552860e+000; krok(n+1)=1.979909e-002; ng(n+1)=2.453151e+002;
n=99; farx(n+1)=1.545253e+000; foe(n+1)=2.443374e+000; krok(n+1)=3.039725e-002; ng(n+1)=3.497814e+002;
n=100; farx(n+1)=1.496893e+000; foe(n+1)=2.299590e+000; krok(n+1)=6.247411e-002; ng(n+1)=4.061880e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)