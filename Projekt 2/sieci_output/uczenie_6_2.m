%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.895105e+003; foe(n+1)=4.808262e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.840313e+003; foe(n+1)=3.633338e+003; krok(n+1)=3.738062e-004; ng(n+1)=4.482049e+003;
n=2; farx(n+1)=1.224186e+003; foe(n+1)=9.521918e+002; krok(n+1)=1.420511e-003; ng(n+1)=3.479817e+003;
n=3; farx(n+1)=1.393936e+003; foe(n+1)=7.999008e+002; krok(n+1)=2.405530e-004; ng(n+1)=3.776949e+003;
n=4; farx(n+1)=4.478667e+002; foe(n+1)=6.354531e+002; krok(n+1)=8.443488e-003; ng(n+1)=7.018922e+002;
n=5; farx(n+1)=3.863914e+002; foe(n+1)=5.704260e+002; krok(n+1)=1.701807e-004; ng(n+1)=5.745084e+003;
n=6; farx(n+1)=3.370926e+002; foe(n+1)=4.933565e+002; krok(n+1)=4.283785e-004; ng(n+1)=3.776864e+003;
n=7; farx(n+1)=2.446706e+002; foe(n+1)=4.043808e+002; krok(n+1)=4.011171e-003; ng(n+1)=1.960897e+003;
n=8; farx(n+1)=2.225217e+002; foe(n+1)=3.873490e+002; krok(n+1)=4.205124e-004; ng(n+1)=1.112673e+003;
n=9; farx(n+1)=2.063252e+002; foe(n+1)=3.629296e+002; krok(n+1)=4.916421e-003; ng(n+1)=1.190365e+003;
n=10; farx(n+1)=2.053910e+002; foe(n+1)=3.445614e+002; krok(n+1)=2.739044e-003; ng(n+1)=7.113146e+002;
n=11; farx(n+1)=2.136767e+002; foe(n+1)=3.220950e+002; krok(n+1)=2.547833e-003; ng(n+1)=8.683217e+002;
n=12; farx(n+1)=2.131212e+002; foe(n+1)=3.216936e+002; krok(n+1)=1.022356e-004; ng(n+1)=6.919514e+002;
n=13; farx(n+1)=2.133035e+002; foe(n+1)=3.143042e+002; krok(n+1)=4.108016e-004; ng(n+1)=1.861529e+003;
n=14; farx(n+1)=2.028420e+002; foe(n+1)=3.097773e+002; krok(n+1)=1.232480e-003; ng(n+1)=7.849222e+002;
n=15; farx(n+1)=1.951933e+002; foe(n+1)=3.082403e+002; krok(n+1)=4.375567e-004; ng(n+1)=3.885835e+003;
n=16; farx(n+1)=1.925228e+002; foe(n+1)=3.071200e+002; krok(n+1)=5.534736e-004; ng(n+1)=3.137887e+003;
n=17; farx(n+1)=1.843131e+002; foe(n+1)=3.016691e+002; krok(n+1)=1.298462e-003; ng(n+1)=8.622081e+003;
n=18; farx(n+1)=1.828980e+002; foe(n+1)=3.011018e+002; krok(n+1)=9.144347e-005; ng(n+1)=1.168298e+004;
n=19; farx(n+1)=1.740681e+002; foe(n+1)=2.961416e+002; krok(n+1)=1.413789e-003; ng(n+1)=8.422273e+002;
n=20; farx(n+1)=1.160309e+002; foe(n+1)=2.697770e+002; krok(n+1)=2.090702e-003; ng(n+1)=2.654865e+004;
n=21; farx(n+1)=1.197987e+002; foe(n+1)=2.539079e+002; krok(n+1)=1.021844e-003; ng(n+1)=6.112404e+003;
n=22; farx(n+1)=1.108236e+002; foe(n+1)=2.210598e+002; krok(n+1)=1.506159e-002; ng(n+1)=6.204397e+002;
n=23; farx(n+1)=1.072289e+002; foe(n+1)=2.190722e+002; krok(n+1)=4.927090e-004; ng(n+1)=2.547925e+003;
n=24; farx(n+1)=9.730319e+001; foe(n+1)=2.091458e+002; krok(n+1)=9.876893e-003; ng(n+1)=1.358239e+003;
n=25; farx(n+1)=9.182940e+001; foe(n+1)=2.044824e+002; krok(n+1)=2.303305e-003; ng(n+1)=3.747295e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=9.233290e+001; foe(n+1)=1.781582e+002; krok(n+1)=2.034152e-005; ng(n+1)=4.151113e+003;
n=27; farx(n+1)=8.052070e+001; foe(n+1)=1.706246e+002; krok(n+1)=2.793956e-005; ng(n+1)=2.501653e+003;
n=28; farx(n+1)=5.849602e+001; foe(n+1)=1.501112e+002; krok(n+1)=2.141893e-004; ng(n+1)=1.314421e+003;
n=29; farx(n+1)=4.070497e+001; foe(n+1)=1.260905e+002; krok(n+1)=3.488003e-005; ng(n+1)=2.887592e+003;
n=30; farx(n+1)=1.719937e+001; foe(n+1)=9.573075e+001; krok(n+1)=6.790684e-004; ng(n+1)=1.903315e+003;
n=31; farx(n+1)=1.709277e+001; foe(n+1)=9.533336e+001; krok(n+1)=2.369246e-005; ng(n+1)=1.218020e+003;
n=32; farx(n+1)=1.611627e+001; foe(n+1)=7.772427e+001; krok(n+1)=1.596151e-004; ng(n+1)=3.760836e+003;
n=33; farx(n+1)=1.534101e+001; foe(n+1)=5.763150e+001; krok(n+1)=3.880097e-003; ng(n+1)=1.978575e+003;
n=34; farx(n+1)=1.379818e+001; foe(n+1)=5.654718e+001; krok(n+1)=1.460638e-004; ng(n+1)=2.395934e+003;
n=35; farx(n+1)=1.120901e+001; foe(n+1)=4.959377e+001; krok(n+1)=7.839040e-004; ng(n+1)=1.966009e+003;
n=36; farx(n+1)=1.154991e+001; foe(n+1)=2.977854e+001; krok(n+1)=3.928638e-003; ng(n+1)=1.818175e+003;
n=37; farx(n+1)=1.159023e+001; foe(n+1)=2.930323e+001; krok(n+1)=1.018217e-004; ng(n+1)=1.940512e+003;
n=38; farx(n+1)=1.126887e+001; foe(n+1)=2.624496e+001; krok(n+1)=4.829711e-003; ng(n+1)=1.555343e+003;
n=39; farx(n+1)=1.132565e+001; foe(n+1)=2.513554e+001; krok(n+1)=9.208067e-004; ng(n+1)=9.936761e+002;
n=40; farx(n+1)=1.162663e+001; foe(n+1)=2.378623e+001; krok(n+1)=6.065270e-003; ng(n+1)=1.445346e+003;
n=41; farx(n+1)=1.241806e+001; foe(n+1)=2.234442e+001; krok(n+1)=1.220118e-002; ng(n+1)=3.288365e+002;
n=42; farx(n+1)=1.239113e+001; foe(n+1)=1.918039e+001; krok(n+1)=1.321526e-002; ng(n+1)=7.008673e+002;
n=43; farx(n+1)=1.171693e+001; foe(n+1)=1.859870e+001; krok(n+1)=4.611188e-003; ng(n+1)=5.258464e+002;
n=44; farx(n+1)=1.173634e+001; foe(n+1)=1.776088e+001; krok(n+1)=1.532510e-002; ng(n+1)=6.874284e+002;
n=45; farx(n+1)=1.077306e+001; foe(n+1)=1.699248e+001; krok(n+1)=1.142218e-002; ng(n+1)=4.739366e+002;
n=46; farx(n+1)=1.112434e+001; foe(n+1)=1.631070e+001; krok(n+1)=2.876670e-002; ng(n+1)=1.111196e+003;
n=47; farx(n+1)=9.788583e+000; foe(n+1)=1.533693e+001; krok(n+1)=5.359529e-002; ng(n+1)=6.949295e+002;
n=48; farx(n+1)=8.948710e+000; foe(n+1)=1.473228e+001; krok(n+1)=3.904156e-002; ng(n+1)=5.843543e+002;
n=49; farx(n+1)=8.828870e+000; foe(n+1)=1.448408e+001; krok(n+1)=4.663470e-002; ng(n+1)=4.848828e+002;
n=50; farx(n+1)=8.246344e+000; foe(n+1)=1.403986e+001; krok(n+1)=4.069805e-002; ng(n+1)=5.718246e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=8.209022e+000; foe(n+1)=1.386480e+001; krok(n+1)=8.010217e-006; ng(n+1)=5.707485e+002;
n=52; farx(n+1)=8.152789e+000; foe(n+1)=1.373221e+001; krok(n+1)=5.385244e-005; ng(n+1)=1.851596e+002;
n=53; farx(n+1)=8.127707e+000; foe(n+1)=1.357810e+001; krok(n+1)=7.295488e-005; ng(n+1)=1.712406e+002;
n=54; farx(n+1)=7.980995e+000; foe(n+1)=1.332499e+001; krok(n+1)=3.970634e-004; ng(n+1)=8.908820e+001;
n=55; farx(n+1)=8.102849e+000; foe(n+1)=1.307203e+001; krok(n+1)=8.398512e-004; ng(n+1)=7.580865e+001;
n=56; farx(n+1)=7.840604e+000; foe(n+1)=1.288559e+001; krok(n+1)=1.889294e-003; ng(n+1)=9.636403e+001;
n=57; farx(n+1)=7.417143e+000; foe(n+1)=1.229542e+001; krok(n+1)=6.100588e-003; ng(n+1)=2.029808e+002;
n=58; farx(n+1)=6.668203e+000; foe(n+1)=1.199923e+001; krok(n+1)=8.533243e-003; ng(n+1)=2.574554e+002;
n=59; farx(n+1)=5.615090e+000; foe(n+1)=1.100698e+001; krok(n+1)=7.870080e-003; ng(n+1)=3.713945e+002;
n=60; farx(n+1)=5.253199e+000; foe(n+1)=1.058214e+001; krok(n+1)=3.050294e-003; ng(n+1)=3.100197e+002;
n=61; farx(n+1)=4.768336e+000; foe(n+1)=1.022477e+001; krok(n+1)=9.939766e-003; ng(n+1)=2.524314e+002;
n=62; farx(n+1)=4.174029e+000; foe(n+1)=9.710050e+000; krok(n+1)=6.271232e-003; ng(n+1)=5.650976e+002;
n=63; farx(n+1)=3.219334e+000; foe(n+1)=8.684158e+000; krok(n+1)=1.428969e-002; ng(n+1)=4.125449e+002;
n=64; farx(n+1)=2.999507e+000; foe(n+1)=8.341303e+000; krok(n+1)=6.606617e-003; ng(n+1)=8.675868e+001;
n=65; farx(n+1)=3.016868e+000; foe(n+1)=8.224529e+000; krok(n+1)=6.425309e-003; ng(n+1)=5.618059e+002;
n=66; farx(n+1)=2.968321e+000; foe(n+1)=7.897738e+000; krok(n+1)=1.647425e-002; ng(n+1)=4.146300e+002;
n=67; farx(n+1)=2.989637e+000; foe(n+1)=7.382231e+000; krok(n+1)=2.972664e-002; ng(n+1)=5.956740e+002;
n=68; farx(n+1)=2.978730e+000; foe(n+1)=7.152762e+000; krok(n+1)=1.438322e-002; ng(n+1)=6.011985e+002;
n=69; farx(n+1)=2.938709e+000; foe(n+1)=6.823454e+000; krok(n+1)=2.128042e-002; ng(n+1)=4.788493e+002;
n=70; farx(n+1)=2.799987e+000; foe(n+1)=6.421529e+000; krok(n+1)=5.945329e-002; ng(n+1)=1.500224e+002;
n=71; farx(n+1)=2.848661e+000; foe(n+1)=6.330194e+000; krok(n+1)=2.426108e-002; ng(n+1)=3.064957e+002;
n=72; farx(n+1)=2.771825e+000; foe(n+1)=6.156110e+000; krok(n+1)=1.072012e-001; ng(n+1)=1.037231e+002;
n=73; farx(n+1)=2.424812e+000; foe(n+1)=5.903555e+000; krok(n+1)=3.202323e-001; ng(n+1)=1.035194e+002;
n=74; farx(n+1)=2.112692e+000; foe(n+1)=5.722894e+000; krok(n+1)=1.497077e-001; ng(n+1)=1.782200e+002;
n=75; farx(n+1)=1.998454e+000; foe(n+1)=5.616600e+000; krok(n+1)=1.134226e-001; ng(n+1)=1.803734e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.993668e+000; foe(n+1)=5.579591e+000; krok(n+1)=1.351734e-005; ng(n+1)=1.809061e+002;
n=77; farx(n+1)=1.995557e+000; foe(n+1)=5.570815e+000; krok(n+1)=1.506057e-005; ng(n+1)=1.045069e+002;
n=78; farx(n+1)=1.997417e+000; foe(n+1)=5.562661e+000; krok(n+1)=1.915690e-005; ng(n+1)=7.164854e+001;
n=79; farx(n+1)=1.999822e+000; foe(n+1)=5.529637e+000; krok(n+1)=4.794822e-004; ng(n+1)=3.497466e+001;
n=80; farx(n+1)=1.977771e+000; foe(n+1)=5.504628e+000; krok(n+1)=8.651099e-004; ng(n+1)=2.367725e+001;
n=81; farx(n+1)=1.978075e+000; foe(n+1)=5.488868e+000; krok(n+1)=3.157424e-003; ng(n+1)=2.103315e+001;
n=82; farx(n+1)=1.982871e+000; foe(n+1)=5.442313e+000; krok(n+1)=1.254988e-002; ng(n+1)=3.025506e+001;
n=83; farx(n+1)=1.936563e+000; foe(n+1)=5.411995e+000; krok(n+1)=4.964322e-003; ng(n+1)=2.827338e+001;
n=84; farx(n+1)=1.889150e+000; foe(n+1)=5.354266e+000; krok(n+1)=2.339865e-002; ng(n+1)=1.315830e+002;
n=85; farx(n+1)=1.866141e+000; foe(n+1)=5.323690e+000; krok(n+1)=3.473809e-003; ng(n+1)=2.410392e+002;
n=86; farx(n+1)=1.844484e+000; foe(n+1)=5.287738e+000; krok(n+1)=1.238474e-002; ng(n+1)=1.732856e+002;
n=87; farx(n+1)=1.846352e+000; foe(n+1)=5.261743e+000; krok(n+1)=1.684825e-002; ng(n+1)=2.512164e+002;
n=88; farx(n+1)=1.792458e+000; foe(n+1)=5.167133e+000; krok(n+1)=5.205007e-002; ng(n+1)=3.216778e+002;
n=89; farx(n+1)=1.783791e+000; foe(n+1)=5.107959e+000; krok(n+1)=1.868912e-002; ng(n+1)=3.410599e+002;
n=90; farx(n+1)=1.755189e+000; foe(n+1)=5.009337e+000; krok(n+1)=3.243320e-002; ng(n+1)=2.500094e+002;
n=91; farx(n+1)=1.679141e+000; foe(n+1)=4.927154e+000; krok(n+1)=4.000365e-002; ng(n+1)=3.591969e+002;
n=92; farx(n+1)=1.601613e+000; foe(n+1)=4.842730e+000; krok(n+1)=3.076926e-002; ng(n+1)=8.130860e+001;
n=93; farx(n+1)=1.503821e+000; foe(n+1)=4.710404e+000; krok(n+1)=6.109779e-002; ng(n+1)=1.670501e+002;
n=94; farx(n+1)=1.279906e+000; foe(n+1)=4.592603e+000; krok(n+1)=1.104861e-001; ng(n+1)=2.725600e+002;
n=95; farx(n+1)=1.144212e+000; foe(n+1)=4.452860e+000; krok(n+1)=1.227474e-001; ng(n+1)=5.098919e+002;
n=96; farx(n+1)=1.022208e+000; foe(n+1)=4.317851e+000; krok(n+1)=6.680629e-002; ng(n+1)=3.891137e+001;
n=97; farx(n+1)=1.019054e+000; foe(n+1)=4.297826e+000; krok(n+1)=1.257761e-002; ng(n+1)=2.810780e+002;
n=98; farx(n+1)=9.471325e-001; foe(n+1)=4.168480e+000; krok(n+1)=1.919924e-001; ng(n+1)=1.463628e+002;
n=99; farx(n+1)=9.268457e-001; foe(n+1)=3.955581e+000; krok(n+1)=2.301336e-001; ng(n+1)=1.709759e+002;
n=100; farx(n+1)=9.097274e-001; foe(n+1)=3.779542e+000; krok(n+1)=6.384007e-002; ng(n+1)=3.072323e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
