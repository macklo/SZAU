%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.032151e+003; foe(n+1)=4.054243e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.907187e+003; foe(n+1)=2.911152e+003; krok(n+1)=4.912716e-004; ng(n+1)=3.648084e+003;
n=2; farx(n+1)=1.205631e+003; foe(n+1)=1.175108e+003; krok(n+1)=2.607264e-003; ng(n+1)=1.794561e+003;
n=3; farx(n+1)=1.383885e+003; foe(n+1)=8.858372e+002; krok(n+1)=1.022356e-004; ng(n+1)=8.764306e+003;
n=4; farx(n+1)=2.347252e+003; foe(n+1)=5.793787e+002; krok(n+1)=6.791933e-004; ng(n+1)=1.027402e+004;
n=5; farx(n+1)=1.424023e+003; foe(n+1)=4.514409e+002; krok(n+1)=3.576264e-003; ng(n+1)=9.540359e+002;
n=6; farx(n+1)=7.265146e+002; foe(n+1)=3.721595e+002; krok(n+1)=1.970836e-003; ng(n+1)=8.243518e+002;
n=7; farx(n+1)=6.550094e+002; foe(n+1)=3.546364e+002; krok(n+1)=6.052270e-004; ng(n+1)=1.295737e+003;
n=8; farx(n+1)=5.954361e+002; foe(n+1)=3.366728e+002; krok(n+1)=3.127663e-003; ng(n+1)=8.053979e+002;
n=9; farx(n+1)=4.475841e+002; foe(n+1)=3.157621e+002; krok(n+1)=9.328627e-003; ng(n+1)=5.333581e+002;
n=10; farx(n+1)=3.479962e+002; foe(n+1)=2.993153e+002; krok(n+1)=4.186106e-003; ng(n+1)=7.592549e+002;
n=11; farx(n+1)=1.671118e+002; foe(n+1)=2.595378e+002; krok(n+1)=8.637826e-003; ng(n+1)=9.428736e+002;
n=12; farx(n+1)=1.621298e+002; foe(n+1)=2.581715e+002; krok(n+1)=1.032284e-004; ng(n+1)=1.251742e+003;
n=13; farx(n+1)=1.488863e+002; foe(n+1)=2.323050e+002; krok(n+1)=5.980537e-003; ng(n+1)=1.707383e+003;
n=14; farx(n+1)=1.489267e+002; foe(n+1)=2.321373e+002; krok(n+1)=6.909928e-007; ng(n+1)=1.655288e+003;
n=15; farx(n+1)=1.486135e+002; foe(n+1)=2.319524e+002; krok(n+1)=2.780862e-005; ng(n+1)=9.969834e+003;
n=16; farx(n+1)=1.229180e+002; foe(n+1)=2.049525e+002; krok(n+1)=9.064776e-004; ng(n+1)=2.998306e+005;
n=17; farx(n+1)=1.215607e+002; foe(n+1)=2.045280e+002; krok(n+1)=4.999771e-006; ng(n+1)=1.853872e+003;
n=18; farx(n+1)=1.154573e+002; foe(n+1)=2.007743e+002; krok(n+1)=9.516502e-004; ng(n+1)=3.108804e+003;
n=19; farx(n+1)=1.014346e+002; foe(n+1)=1.683261e+002; krok(n+1)=1.086709e-002; ng(n+1)=2.306331e+003;
n=20; farx(n+1)=1.108426e+002; foe(n+1)=1.591490e+002; krok(n+1)=8.713055e-004; ng(n+1)=2.824406e+003;
n=21; farx(n+1)=1.148165e+002; foe(n+1)=1.551658e+002; krok(n+1)=3.660956e-003; ng(n+1)=3.311593e+003;
n=22; farx(n+1)=9.949251e+001; foe(n+1)=1.411153e+002; krok(n+1)=2.388829e-003; ng(n+1)=3.083892e+003;
n=23; farx(n+1)=8.571706e+001; foe(n+1)=1.358323e+002; krok(n+1)=1.965086e-003; ng(n+1)=3.033132e+003;
n=24; farx(n+1)=6.321832e+001; foe(n+1)=1.287852e+002; krok(n+1)=6.475792e-003; ng(n+1)=3.158413e+003;
n=25; farx(n+1)=6.732288e+001; foe(n+1)=1.190194e+002; krok(n+1)=2.107682e-002; ng(n+1)=3.295428e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=8.314523e+001; foe(n+1)=8.821697e+001; krok(n+1)=3.022327e-005; ng(n+1)=3.772443e+003;
n=27; farx(n+1)=7.564572e+001; foe(n+1)=8.470668e+001; krok(n+1)=2.832743e-005; ng(n+1)=1.555804e+003;
n=28; farx(n+1)=5.131365e+001; foe(n+1)=6.578033e+001; krok(n+1)=2.288097e-004; ng(n+1)=2.009932e+003;
n=29; farx(n+1)=3.713510e+001; foe(n+1)=4.411387e+001; krok(n+1)=4.449379e-004; ng(n+1)=1.095767e+003;
n=30; farx(n+1)=2.960406e+001; foe(n+1)=3.824544e+001; krok(n+1)=1.415702e-003; ng(n+1)=6.751058e+002;
n=31; farx(n+1)=2.230679e+001; foe(n+1)=3.299194e+001; krok(n+1)=8.861138e-004; ng(n+1)=5.983958e+002;
n=32; farx(n+1)=1.936890e+001; foe(n+1)=3.081685e+001; krok(n+1)=1.265711e-003; ng(n+1)=1.887393e+003;
n=33; farx(n+1)=1.592916e+001; foe(n+1)=2.644821e+001; krok(n+1)=4.683883e-003; ng(n+1)=1.439647e+003;
n=34; farx(n+1)=1.558947e+001; foe(n+1)=2.589592e+001; krok(n+1)=3.240245e-003; ng(n+1)=5.272002e+002;
n=35; farx(n+1)=1.442711e+001; foe(n+1)=2.498937e+001; krok(n+1)=8.799720e-003; ng(n+1)=1.003209e+003;
n=36; farx(n+1)=1.245491e+001; foe(n+1)=2.264103e+001; krok(n+1)=1.567500e-002; ng(n+1)=1.397449e+003;
n=37; farx(n+1)=1.085003e+001; foe(n+1)=2.170576e+001; krok(n+1)=8.424126e-003; ng(n+1)=2.869499e+002;
n=38; farx(n+1)=1.018736e+001; foe(n+1)=2.136610e+001; krok(n+1)=3.426502e-003; ng(n+1)=3.966289e+002;
n=39; farx(n+1)=9.278995e+000; foe(n+1)=2.057103e+001; krok(n+1)=4.049572e-003; ng(n+1)=3.896507e+002;
n=40; farx(n+1)=8.878284e+000; foe(n+1)=1.991378e+001; krok(n+1)=9.317111e-003; ng(n+1)=6.030014e+002;
n=41; farx(n+1)=7.784396e+000; foe(n+1)=1.896662e+001; krok(n+1)=1.370811e-002; ng(n+1)=2.664582e+002;
n=42; farx(n+1)=6.502107e+000; foe(n+1)=1.784209e+001; krok(n+1)=1.547432e-002; ng(n+1)=3.724978e+002;
n=43; farx(n+1)=6.268193e+000; foe(n+1)=1.758632e+001; krok(n+1)=5.767896e-003; ng(n+1)=7.419775e+002;
n=44; farx(n+1)=5.747164e+000; foe(n+1)=1.679204e+001; krok(n+1)=3.863769e-002; ng(n+1)=3.133875e+002;
n=45; farx(n+1)=5.297490e+000; foe(n+1)=1.620879e+001; krok(n+1)=1.552269e-002; ng(n+1)=3.083515e+002;
n=46; farx(n+1)=5.101969e+000; foe(n+1)=1.543855e+001; krok(n+1)=3.535818e-002; ng(n+1)=6.065388e+002;
n=47; farx(n+1)=4.887541e+000; foe(n+1)=1.486034e+001; krok(n+1)=1.663330e-001; ng(n+1)=3.697041e+002;
n=48; farx(n+1)=4.796296e+000; foe(n+1)=1.434561e+001; krok(n+1)=7.853583e-002; ng(n+1)=2.991345e+002;
n=49; farx(n+1)=4.774475e+000; foe(n+1)=1.349209e+001; krok(n+1)=1.158784e-001; ng(n+1)=4.195952e+002;
n=50; farx(n+1)=4.527943e+000; foe(n+1)=1.277020e+001; krok(n+1)=2.522670e-001; ng(n+1)=1.702039e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=4.493824e+000; foe(n+1)=1.261334e+001; krok(n+1)=4.279756e-005; ng(n+1)=2.496596e+002;
n=52; farx(n+1)=4.495321e+000; foe(n+1)=1.260064e+001; krok(n+1)=3.517612e-005; ng(n+1)=8.069263e+001;
n=53; farx(n+1)=4.560806e+000; foe(n+1)=1.222274e+001; krok(n+1)=2.890523e-004; ng(n+1)=1.457753e+002;
n=54; farx(n+1)=4.517405e+000; foe(n+1)=1.205777e+001; krok(n+1)=8.489916e-005; ng(n+1)=1.538364e+002;
n=55; farx(n+1)=4.384083e+000; foe(n+1)=1.152393e+001; krok(n+1)=2.947309e-003; ng(n+1)=5.747883e+001;
n=56; farx(n+1)=4.319807e+000; foe(n+1)=1.142310e+001; krok(n+1)=1.812955e-003; ng(n+1)=4.664913e+001;
n=57; farx(n+1)=4.122495e+000; foe(n+1)=1.094530e+001; krok(n+1)=4.866267e-003; ng(n+1)=5.351951e+001;
n=58; farx(n+1)=3.775852e+000; foe(n+1)=1.063665e+001; krok(n+1)=1.333847e-002; ng(n+1)=1.518134e+002;
n=59; farx(n+1)=3.687837e+000; foe(n+1)=1.041793e+001; krok(n+1)=7.874995e-003; ng(n+1)=1.703071e+002;
n=60; farx(n+1)=3.300284e+000; foe(n+1)=9.880535e+000; krok(n+1)=2.414071e-002; ng(n+1)=1.727874e+002;
n=61; farx(n+1)=3.150462e+000; foe(n+1)=9.638985e+000; krok(n+1)=3.303815e-003; ng(n+1)=3.343243e+002;
n=62; farx(n+1)=2.929556e+000; foe(n+1)=9.253743e+000; krok(n+1)=1.806042e-002; ng(n+1)=1.747998e+002;
n=63; farx(n+1)=2.791264e+000; foe(n+1)=9.122554e+000; krok(n+1)=6.227174e-003; ng(n+1)=4.352660e+002;
n=64; farx(n+1)=2.432078e+000; foe(n+1)=8.472981e+000; krok(n+1)=9.506102e-003; ng(n+1)=2.486991e+002;
n=65; farx(n+1)=2.374258e+000; foe(n+1)=8.330802e+000; krok(n+1)=5.819038e-003; ng(n+1)=1.730574e+002;
n=66; farx(n+1)=2.216291e+000; foe(n+1)=7.944842e+000; krok(n+1)=3.850210e-002; ng(n+1)=2.412544e+002;
n=67; farx(n+1)=2.098067e+000; foe(n+1)=7.703008e+000; krok(n+1)=8.174750e-003; ng(n+1)=6.236123e+002;
n=68; farx(n+1)=2.016469e+000; foe(n+1)=7.437726e+000; krok(n+1)=2.896961e-002; ng(n+1)=6.571863e+002;
n=69; farx(n+1)=1.989018e+000; foe(n+1)=7.017222e+000; krok(n+1)=2.617232e-002; ng(n+1)=6.016781e+002;
n=70; farx(n+1)=1.949698e+000; foe(n+1)=5.994226e+000; krok(n+1)=6.786752e-002; ng(n+1)=3.772578e+002;
n=71; farx(n+1)=1.940616e+000; foe(n+1)=5.902859e+000; krok(n+1)=1.184956e-002; ng(n+1)=4.418539e+002;
n=72; farx(n+1)=1.905809e+000; foe(n+1)=5.737987e+000; krok(n+1)=5.160128e-002; ng(n+1)=3.017709e+002;
n=73; farx(n+1)=1.879362e+000; foe(n+1)=5.572226e+000; krok(n+1)=5.241643e-003; ng(n+1)=7.295291e+002;
n=74; farx(n+1)=1.862292e+000; foe(n+1)=5.358746e+000; krok(n+1)=1.280543e-001; ng(n+1)=5.012005e+002;
n=75; farx(n+1)=1.807933e+000; foe(n+1)=5.119604e+000; krok(n+1)=1.668649e-001; ng(n+1)=6.564238e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.807528e+000; foe(n+1)=5.052527e+000; krok(n+1)=1.021552e-006; ng(n+1)=1.009758e+003;
n=77; farx(n+1)=1.804392e+000; foe(n+1)=4.913528e+000; krok(n+1)=1.929081e-005; ng(n+1)=3.226068e+002;
n=78; farx(n+1)=1.796275e+000; foe(n+1)=4.897696e+000; krok(n+1)=4.621507e-005; ng(n+1)=7.886007e+001;
n=79; farx(n+1)=1.787862e+000; foe(n+1)=4.829202e+000; krok(n+1)=6.694872e-004; ng(n+1)=4.905841e+001;
n=80; farx(n+1)=1.745906e+000; foe(n+1)=4.728011e+000; krok(n+1)=1.104943e-003; ng(n+1)=5.592445e+001;
n=81; farx(n+1)=1.714911e+000; foe(n+1)=4.603264e+000; krok(n+1)=6.118871e-003; ng(n+1)=2.190109e+002;
n=82; farx(n+1)=1.597248e+000; foe(n+1)=4.220155e+000; krok(n+1)=6.427616e-003; ng(n+1)=3.309564e+002;
n=83; farx(n+1)=1.540183e+000; foe(n+1)=3.965313e+000; krok(n+1)=6.303106e-003; ng(n+1)=2.924802e+002;
n=84; farx(n+1)=1.515224e+000; foe(n+1)=3.802655e+000; krok(n+1)=2.778799e-003; ng(n+1)=2.816339e+002;
n=85; farx(n+1)=1.514671e+000; foe(n+1)=3.671950e+000; krok(n+1)=3.816075e-003; ng(n+1)=6.933513e+002;
n=86; farx(n+1)=1.511354e+000; foe(n+1)=3.579047e+000; krok(n+1)=4.258262e-003; ng(n+1)=4.055821e+002;
n=87; farx(n+1)=1.483350e+000; foe(n+1)=3.375017e+000; krok(n+1)=7.442795e-003; ng(n+1)=5.159915e+002;
n=88; farx(n+1)=1.466701e+000; foe(n+1)=3.296887e+000; krok(n+1)=3.660956e-003; ng(n+1)=3.683356e+002;
n=89; farx(n+1)=1.405478e+000; foe(n+1)=3.099075e+000; krok(n+1)=1.953985e-002; ng(n+1)=6.193407e+002;
n=90; farx(n+1)=1.371520e+000; foe(n+1)=3.001680e+000; krok(n+1)=1.669856e-002; ng(n+1)=3.442991e+002;
n=91; farx(n+1)=1.351424e+000; foe(n+1)=2.932066e+000; krok(n+1)=9.913504e-003; ng(n+1)=3.633527e+002;
n=92; farx(n+1)=1.344169e+000; foe(n+1)=2.803565e+000; krok(n+1)=6.255374e-002; ng(n+1)=2.918068e+002;
n=93; farx(n+1)=1.328708e+000; foe(n+1)=2.573657e+000; krok(n+1)=5.121364e-002; ng(n+1)=5.961037e+002;
n=94; farx(n+1)=1.393731e+000; foe(n+1)=2.459740e+000; krok(n+1)=9.137743e-002; ng(n+1)=9.838725e+001;
n=95; farx(n+1)=1.377965e+000; foe(n+1)=2.396846e+000; krok(n+1)=6.046686e-002; ng(n+1)=1.854326e+002;
n=96; farx(n+1)=1.308640e+000; foe(n+1)=2.266138e+000; krok(n+1)=7.213047e-002; ng(n+1)=6.571565e+001;
n=97; farx(n+1)=1.346958e+000; foe(n+1)=2.188689e+000; krok(n+1)=5.670481e-002; ng(n+1)=2.064892e+002;
n=98; farx(n+1)=1.293641e+000; foe(n+1)=2.112774e+000; krok(n+1)=1.865214e-001; ng(n+1)=3.065485e+002;
n=99; farx(n+1)=1.301354e+000; foe(n+1)=2.050500e+000; krok(n+1)=1.241815e-001; ng(n+1)=3.987074e+002;
n=100; farx(n+1)=1.255661e+000; foe(n+1)=1.994782e+000; krok(n+1)=2.064051e-001; ng(n+1)=1.544308e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)