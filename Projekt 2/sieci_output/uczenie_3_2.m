%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.070929e+003; foe(n+1)=4.054548e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.916690e+003; foe(n+1)=2.997119e+003; krok(n+1)=4.927090e-004; ng(n+1)=2.666040e+003;
n=2; farx(n+1)=8.497573e+002; foe(n+1)=9.779125e+002; krok(n+1)=5.548179e-003; ng(n+1)=7.988982e+002;
n=3; farx(n+1)=9.412580e+002; foe(n+1)=8.354610e+002; krok(n+1)=5.561724e-005; ng(n+1)=6.245231e+003;
n=4; farx(n+1)=1.391186e+003; foe(n+1)=6.176821e+002; krok(n+1)=2.834094e-004; ng(n+1)=7.536942e+003;
n=5; farx(n+1)=1.520209e+003; foe(n+1)=5.414610e+002; krok(n+1)=2.581388e-003; ng(n+1)=1.148301e+003;
n=6; farx(n+1)=1.112271e+003; foe(n+1)=4.951815e+002; krok(n+1)=2.841023e-003; ng(n+1)=1.404323e+003;
n=7; farx(n+1)=5.061696e+002; foe(n+1)=3.893430e+002; krok(n+1)=1.111520e-002; ng(n+1)=2.385296e+003;
n=8; farx(n+1)=4.239624e+002; foe(n+1)=3.781417e+002; krok(n+1)=1.244425e-005; ng(n+1)=3.055152e+003;
n=9; farx(n+1)=3.884613e+002; foe(n+1)=3.672579e+002; krok(n+1)=9.824964e-004; ng(n+1)=2.499467e+003;
n=10; farx(n+1)=3.773598e+002; foe(n+1)=3.615337e+002; krok(n+1)=8.408970e-005; ng(n+1)=2.814196e+003;
n=11; farx(n+1)=3.574228e+002; foe(n+1)=3.578647e+002; krok(n+1)=2.435962e-003; ng(n+1)=1.448932e+003;
n=12; farx(n+1)=3.358838e+002; foe(n+1)=3.516705e+002; krok(n+1)=2.722891e-003; ng(n+1)=1.204799e+003;
n=13; farx(n+1)=3.145064e+002; foe(n+1)=3.486595e+002; krok(n+1)=2.100547e-003; ng(n+1)=8.789164e+002;
n=14; farx(n+1)=3.018061e+002; foe(n+1)=3.366688e+002; krok(n+1)=1.887308e-002; ng(n+1)=8.363965e+002;
n=15; farx(n+1)=3.049698e+002; foe(n+1)=3.302907e+002; krok(n+1)=1.464382e-002; ng(n+1)=1.145767e+003;
n=16; farx(n+1)=2.140349e+002; foe(n+1)=3.213929e+002; krok(n+1)=2.097129e-002; ng(n+1)=1.625159e+003;
n=17; farx(n+1)=1.145583e+002; foe(n+1)=2.722631e+002; krok(n+1)=1.392689e-001; ng(n+1)=1.653651e+003;
n=18; farx(n+1)=1.005552e+002; foe(n+1)=2.555940e+002; krok(n+1)=3.108824e-002; ng(n+1)=2.943997e+003;
n=19; farx(n+1)=1.148146e+002; foe(n+1)=2.404780e+002; krok(n+1)=1.394132e-002; ng(n+1)=6.600246e+003;
n=20; farx(n+1)=4.990304e+001; foe(n+1)=1.737355e+002; krok(n+1)=9.372047e-001; ng(n+1)=3.806111e+003;
n=21; farx(n+1)=4.422196e+001; foe(n+1)=1.711784e+002; krok(n+1)=1.142218e-002; ng(n+1)=4.838227e+003;
n=22; farx(n+1)=3.887311e+001; foe(n+1)=1.523096e+002; krok(n+1)=1.158784e-001; ng(n+1)=5.897558e+003;
n=23; farx(n+1)=3.370088e+001; foe(n+1)=1.282167e+002; krok(n+1)=7.406706e-002; ng(n+1)=4.013947e+003;
n=24; farx(n+1)=2.942519e+001; foe(n+1)=1.165989e+002; krok(n+1)=2.896961e-002; ng(n+1)=4.876154e+003;
n=25; farx(n+1)=2.400186e+001; foe(n+1)=1.081372e+002; krok(n+1)=6.647460e-002; ng(n+1)=6.563120e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=1.960620e+001; foe(n+1)=8.097144e+001; krok(n+1)=3.525905e-006; ng(n+1)=9.370298e+003;
n=27; farx(n+1)=1.982038e+001; foe(n+1)=7.665801e+001; krok(n+1)=5.161420e-005; ng(n+1)=8.412384e+002;
n=28; farx(n+1)=1.822797e+001; foe(n+1)=7.139126e+001; krok(n+1)=1.523335e-004; ng(n+1)=5.966663e+002;
n=29; farx(n+1)=1.797701e+001; foe(n+1)=6.582155e+001; krok(n+1)=2.927427e-004; ng(n+1)=5.674546e+002;
n=30; farx(n+1)=1.515824e+001; foe(n+1)=5.597801e+001; krok(n+1)=1.499941e-003; ng(n+1)=3.554039e+002;
n=31; farx(n+1)=1.486826e+001; foe(n+1)=5.169857e+001; krok(n+1)=2.476023e-003; ng(n+1)=2.752237e+003;
n=32; farx(n+1)=1.248753e+001; foe(n+1)=4.374327e+001; krok(n+1)=7.737158e-003; ng(n+1)=1.179380e+003;
n=33; farx(n+1)=8.656506e+000; foe(n+1)=3.885033e+001; krok(n+1)=9.240905e-003; ng(n+1)=2.663115e+002;
n=34; farx(n+1)=8.380695e+000; foe(n+1)=3.823079e+001; krok(n+1)=1.586441e-002; ng(n+1)=1.233772e+003;
n=35; farx(n+1)=7.890971e+000; foe(n+1)=3.622396e+001; krok(n+1)=2.194150e-002; ng(n+1)=7.003841e+002;
n=36; farx(n+1)=7.816310e+000; foe(n+1)=3.596651e+001; krok(n+1)=7.824004e-003; ng(n+1)=1.540120e+003;
n=37; farx(n+1)=8.201043e+000; foe(n+1)=3.530721e+001; krok(n+1)=8.252113e-003; ng(n+1)=2.043119e+003;
n=38; farx(n+1)=8.503029e+000; foe(n+1)=3.478796e+001; krok(n+1)=6.023109e-003; ng(n+1)=1.662779e+003;
n=39; farx(n+1)=8.391893e+000; foe(n+1)=3.424135e+001; krok(n+1)=3.571281e-002; ng(n+1)=1.511004e+003;
n=40; farx(n+1)=7.653803e+000; foe(n+1)=3.268080e+001; krok(n+1)=5.131523e-002; ng(n+1)=1.150017e+003;
n=41; farx(n+1)=6.724810e+000; foe(n+1)=2.771104e+001; krok(n+1)=1.717868e-001; ng(n+1)=4.434705e+002;
n=42; farx(n+1)=5.964494e+000; foe(n+1)=2.552441e+001; krok(n+1)=3.457398e-001; ng(n+1)=1.379379e+003;
n=43; farx(n+1)=5.475375e+000; foe(n+1)=2.097477e+001; krok(n+1)=8.923142e-001; ng(n+1)=4.330537e+002;
n=44; farx(n+1)=4.868832e+000; foe(n+1)=1.898468e+001; krok(n+1)=1.817918e-001; ng(n+1)=1.082528e+003;
n=45; farx(n+1)=3.699040e+000; foe(n+1)=1.618423e+001; krok(n+1)=5.403454e-001; ng(n+1)=4.689187e+002;
n=46; farx(n+1)=3.466195e+000; foe(n+1)=1.426186e+001; krok(n+1)=1.486329e-001; ng(n+1)=4.729895e+002;
n=47; farx(n+1)=3.438545e+000; foe(n+1)=1.114706e+001; krok(n+1)=2.893371e-001; ng(n+1)=1.440658e+003;
n=48; farx(n+1)=3.126857e+000; foe(n+1)=9.346793e+000; krok(n+1)=2.725909e-001; ng(n+1)=5.286382e+002;
n=49; farx(n+1)=2.484325e+000; foe(n+1)=7.247394e+000; krok(n+1)=5.962886e-001; ng(n+1)=7.113800e+002;
n=50; farx(n+1)=2.431557e+000; foe(n+1)=6.986259e+000; krok(n+1)=5.183452e-001; ng(n+1)=3.134076e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=2.420514e+000; foe(n+1)=6.913419e+000; krok(n+1)=4.273013e-006; ng(n+1)=3.834533e+002;
n=52; farx(n+1)=2.410400e+000; foe(n+1)=6.878059e+000; krok(n+1)=3.273007e-005; ng(n+1)=1.068978e+002;
n=53; farx(n+1)=2.401958e+000; foe(n+1)=6.862654e+000; krok(n+1)=6.549403e-005; ng(n+1)=5.473360e+001;
n=54; farx(n+1)=2.371622e+000; foe(n+1)=6.709588e+000; krok(n+1)=5.188110e-004; ng(n+1)=6.070314e+001;
n=55; farx(n+1)=2.386163e+000; foe(n+1)=6.535580e+000; krok(n+1)=5.201755e-003; ng(n+1)=2.263318e+001;
n=56; farx(n+1)=2.320457e+000; foe(n+1)=6.440576e+000; krok(n+1)=2.254077e-003; ng(n+1)=1.451391e+002;
n=57; farx(n+1)=2.348010e+000; foe(n+1)=6.272259e+000; krok(n+1)=1.920518e-002; ng(n+1)=1.980046e+002;
n=58; farx(n+1)=2.314574e+000; foe(n+1)=6.248187e+000; krok(n+1)=1.165911e-002; ng(n+1)=1.475096e+002;
n=59; farx(n+1)=2.255482e+000; foe(n+1)=6.190762e+000; krok(n+1)=1.352965e-002; ng(n+1)=1.911511e+002;
n=60; farx(n+1)=2.252392e+000; foe(n+1)=6.120779e+000; krok(n+1)=1.920145e-002; ng(n+1)=2.039389e+002;
n=61; farx(n+1)=2.241838e+000; foe(n+1)=6.056959e+000; krok(n+1)=1.031518e-001; ng(n+1)=2.884265e+002;
n=62; farx(n+1)=2.316677e+000; foe(n+1)=5.993478e+000; krok(n+1)=7.173481e-002; ng(n+1)=4.551130e+002;
n=63; farx(n+1)=2.401329e+000; foe(n+1)=5.724398e+000; krok(n+1)=9.636975e-002; ng(n+1)=3.627742e+002;
n=64; farx(n+1)=2.247021e+000; foe(n+1)=5.650657e+000; krok(n+1)=4.256083e-002; ng(n+1)=1.413526e+002;
n=65; farx(n+1)=2.252648e+000; foe(n+1)=5.615310e+000; krok(n+1)=9.407061e-002; ng(n+1)=1.820200e+002;
n=66; farx(n+1)=1.977108e+000; foe(n+1)=5.265808e+000; krok(n+1)=6.614848e-001; ng(n+1)=2.468839e+002;
n=67; farx(n+1)=1.859307e+000; foe(n+1)=4.984518e+000; krok(n+1)=6.480438e-001; ng(n+1)=3.460225e+002;
n=68; farx(n+1)=1.629291e+000; foe(n+1)=4.567527e+000; krok(n+1)=7.113726e-001; ng(n+1)=2.773657e+002;
n=69; farx(n+1)=1.575737e+000; foe(n+1)=4.438119e+000; krok(n+1)=2.418297e-001; ng(n+1)=2.113187e+002;
n=70; farx(n+1)=1.519324e+000; foe(n+1)=4.185493e+000; krok(n+1)=8.536620e-001; ng(n+1)=1.525912e+002;
n=71; farx(n+1)=1.527844e+000; foe(n+1)=4.111767e+000; krok(n+1)=3.435737e-001; ng(n+1)=9.423362e+001;
n=72; farx(n+1)=1.535215e+000; foe(n+1)=4.051683e+000; krok(n+1)=2.507212e-001; ng(n+1)=9.264568e+001;
n=73; farx(n+1)=1.525018e+000; foe(n+1)=4.005365e+000; krok(n+1)=5.327344e-001; ng(n+1)=9.266667e+001;
n=74; farx(n+1)=1.545700e+000; foe(n+1)=3.970287e+000; krok(n+1)=8.340965e-001; ng(n+1)=5.490541e+001;
n=75; farx(n+1)=1.558475e+000; foe(n+1)=3.942700e+000; krok(n+1)=1.492366e+000; ng(n+1)=8.904371e+001;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.559203e+000; foe(n+1)=3.941055e+000; krok(n+1)=4.568717e-006; ng(n+1)=4.824225e+001;
n=77; farx(n+1)=1.559021e+000; foe(n+1)=3.940298e+000; krok(n+1)=1.099167e-005; ng(n+1)=2.414270e+001;
n=78; farx(n+1)=1.557372e+000; foe(n+1)=3.938815e+000; krok(n+1)=2.596664e-004; ng(n+1)=7.248188e+000;
n=79; farx(n+1)=1.554275e+000; foe(n+1)=3.922703e+000; krok(n+1)=4.465528e-004; ng(n+1)=1.774528e+001;
n=80; farx(n+1)=1.557160e+000; foe(n+1)=3.907445e+000; krok(n+1)=3.255599e-003; ng(n+1)=7.344148e+000;
n=81; farx(n+1)=1.555564e+000; foe(n+1)=3.868204e+000; krok(n+1)=9.564855e-003; ng(n+1)=9.261698e+000;
n=82; farx(n+1)=1.542273e+000; foe(n+1)=3.848572e+000; krok(n+1)=3.364099e-003; ng(n+1)=6.342248e+001;
n=83; farx(n+1)=1.533607e+000; foe(n+1)=3.836310e+000; krok(n+1)=6.543079e-003; ng(n+1)=1.152739e+002;
n=84; farx(n+1)=1.551186e+000; foe(n+1)=3.822733e+000; krok(n+1)=3.192743e-002; ng(n+1)=1.502077e+002;
n=85; farx(n+1)=1.554025e+000; foe(n+1)=3.811013e+000; krok(n+1)=1.726831e-002; ng(n+1)=1.648616e+002;
n=86; farx(n+1)=1.564093e+000; foe(n+1)=3.797270e+000; krok(n+1)=1.925105e-002; ng(n+1)=1.570084e+002;
n=87; farx(n+1)=1.561169e+000; foe(n+1)=3.790289e+000; krok(n+1)=3.034235e-002; ng(n+1)=1.467170e+002;
n=88; farx(n+1)=1.550598e+000; foe(n+1)=3.766208e+000; krok(n+1)=1.409038e-001; ng(n+1)=1.345684e+002;
n=89; farx(n+1)=1.535361e+000; foe(n+1)=3.712349e+000; krok(n+1)=6.754317e-002; ng(n+1)=1.818650e+002;
n=90; farx(n+1)=1.516534e+000; foe(n+1)=3.680939e+000; krok(n+1)=9.118713e-002; ng(n+1)=9.124902e+001;
n=91; farx(n+1)=1.497077e+000; foe(n+1)=3.608804e+000; krok(n+1)=8.536620e-001; ng(n+1)=5.887979e+001;
n=92; farx(n+1)=1.495160e+000; foe(n+1)=3.573387e+000; krok(n+1)=7.807902e-001; ng(n+1)=1.135342e+002;
n=93; farx(n+1)=1.475406e+000; foe(n+1)=3.545341e+000; krok(n+1)=6.228822e-001; ng(n+1)=1.115295e+002;
n=94; farx(n+1)=1.463005e+000; foe(n+1)=3.517808e+000; krok(n+1)=6.101385e-001; ng(n+1)=8.235557e+001;
n=95; farx(n+1)=1.438697e+000; foe(n+1)=3.470319e+000; krok(n+1)=5.022690e-001; ng(n+1)=1.257038e+002;
n=96; farx(n+1)=1.422214e+000; foe(n+1)=3.440397e+000; krok(n+1)=4.536903e-001; ng(n+1)=1.362872e+002;
n=97; farx(n+1)=1.380832e+000; foe(n+1)=3.386817e+000; krok(n+1)=9.205262e-001; ng(n+1)=1.235160e+002;
n=98; farx(n+1)=1.349395e+000; foe(n+1)=3.306566e+000; krok(n+1)=7.725027e-001; ng(n+1)=1.282878e+002;
n=99; farx(n+1)=1.329200e+000; foe(n+1)=3.212325e+000; krok(n+1)=8.179879e-001; ng(n+1)=2.719565e+002;
n=100; farx(n+1)=1.327858e+000; foe(n+1)=3.165929e+000; krok(n+1)=1.181740e-001; ng(n+1)=1.427207e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
