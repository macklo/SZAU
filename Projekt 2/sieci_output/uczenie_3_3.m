%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.920393e+003; foe(n+1)=4.933380e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.925086e+003; foe(n+1)=3.925700e+003; krok(n+1)=3.383827e-004; ng(n+1)=3.113277e+003;
n=2; farx(n+1)=2.238249e+003; foe(n+1)=1.529184e+003; krok(n+1)=1.680696e-002; ng(n+1)=3.704503e+002;
n=3; farx(n+1)=1.965308e+003; foe(n+1)=1.449308e+003; krok(n+1)=3.586316e-006; ng(n+1)=1.105317e+004;
n=4; farx(n+1)=3.109800e+003; foe(n+1)=1.035803e+003; krok(n+1)=1.367364e-004; ng(n+1)=7.929398e+003;
n=5; farx(n+1)=2.560119e+003; foe(n+1)=9.526869e+002; krok(n+1)=1.280208e-004; ng(n+1)=3.163569e+003;
n=6; farx(n+1)=3.707620e+003; foe(n+1)=9.097896e+002; krok(n+1)=8.471113e-004; ng(n+1)=1.204481e+003;
n=7; farx(n+1)=3.676602e+003; foe(n+1)=8.857513e+002; krok(n+1)=4.749570e-004; ng(n+1)=1.588193e+003;
n=8; farx(n+1)=3.485969e+003; foe(n+1)=8.557099e+002; krok(n+1)=3.604935e-004; ng(n+1)=2.642433e+003;
n=9; farx(n+1)=3.400666e+003; foe(n+1)=8.234572e+002; krok(n+1)=2.213894e-003; ng(n+1)=2.620497e+003;
n=10; farx(n+1)=3.179574e+003; foe(n+1)=8.049255e+002; krok(n+1)=4.575114e-003; ng(n+1)=3.065843e+003;
n=11; farx(n+1)=2.541085e+003; foe(n+1)=7.810610e+002; krok(n+1)=6.113282e-003; ng(n+1)=3.551570e+003;
n=12; farx(n+1)=3.232800e+003; foe(n+1)=7.441276e+002; krok(n+1)=1.767909e-002; ng(n+1)=4.992040e+003;
n=13; farx(n+1)=3.491022e+003; foe(n+1)=7.387927e+002; krok(n+1)=1.161194e-003; ng(n+1)=1.956827e+003;
n=14; farx(n+1)=3.299575e+003; foe(n+1)=7.233697e+002; krok(n+1)=2.409244e-002; ng(n+1)=1.251896e+003;
n=15; farx(n+1)=3.003740e+003; foe(n+1)=7.097561e+002; krok(n+1)=5.647473e-003; ng(n+1)=1.356713e+003;
n=16; farx(n+1)=2.438367e+003; foe(n+1)=6.846995e+002; krok(n+1)=1.576669e-002; ng(n+1)=2.131950e+003;
n=17; farx(n+1)=2.150358e+003; foe(n+1)=6.770676e+002; krok(n+1)=1.207035e-002; ng(n+1)=3.829646e+003;
n=18; farx(n+1)=1.929351e+003; foe(n+1)=6.107457e+002; krok(n+1)=1.809047e-002; ng(n+1)=4.821312e+003;
n=19; farx(n+1)=1.197185e+003; foe(n+1)=5.471998e+002; krok(n+1)=2.847182e-002; ng(n+1)=3.363501e+003;
n=20; farx(n+1)=9.361613e+002; foe(n+1)=5.301522e+002; krok(n+1)=8.589342e-002; ng(n+1)=4.077922e+003;
n=21; farx(n+1)=7.463607e+002; foe(n+1)=4.933699e+002; krok(n+1)=2.114117e-001; ng(n+1)=4.637960e+003;
n=22; farx(n+1)=6.010907e+002; foe(n+1)=4.743369e+002; krok(n+1)=2.954196e-001; ng(n+1)=1.426214e+003;
n=23; farx(n+1)=5.684942e+002; foe(n+1)=4.672924e+002; krok(n+1)=3.920841e-002; ng(n+1)=8.638295e+002;
n=24; farx(n+1)=4.951338e+002; foe(n+1)=4.407848e+002; krok(n+1)=3.307424e-001; ng(n+1)=1.490383e+003;
n=25; farx(n+1)=2.329769e+002; foe(n+1)=4.065360e+002; krok(n+1)=7.491049e-001; ng(n+1)=1.176603e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=2.193088e+002; foe(n+1)=3.874263e+002; krok(n+1)=6.370154e-006; ng(n+1)=4.619797e+003;
n=27; farx(n+1)=2.169257e+002; foe(n+1)=3.830857e+002; krok(n+1)=3.006731e-006; ng(n+1)=3.510790e+003;
n=28; farx(n+1)=1.158745e+002; foe(n+1)=3.286717e+002; krok(n+1)=2.775430e-005; ng(n+1)=2.665197e+003;
n=29; farx(n+1)=8.715014e+001; foe(n+1)=3.078329e+002; krok(n+1)=5.676782e-004; ng(n+1)=4.255456e+003;
n=30; farx(n+1)=8.373256e+001; foe(n+1)=3.054500e+002; krok(n+1)=2.045351e-004; ng(n+1)=8.542672e+003;
n=31; farx(n+1)=4.940640e+001; foe(n+1)=2.732433e+002; krok(n+1)=6.847610e-004; ng(n+1)=8.574699e+003;
n=32; farx(n+1)=2.717410e+001; foe(n+1)=2.436678e+002; krok(n+1)=5.631478e-004; ng(n+1)=1.138616e+004;
n=33; farx(n+1)=1.484126e+001; foe(n+1)=1.646975e+002; krok(n+1)=1.347687e-003; ng(n+1)=1.655579e+004;
n=34; farx(n+1)=1.484453e+001; foe(n+1)=1.528534e+002; krok(n+1)=9.237209e-006; ng(n+1)=1.039225e+004;
n=35; farx(n+1)=1.513050e+001; foe(n+1)=1.491864e+002; krok(n+1)=1.858436e-004; ng(n+1)=1.204982e+004;
n=36; farx(n+1)=1.406577e+001; foe(n+1)=1.434403e+002; krok(n+1)=8.631728e-004; ng(n+1)=9.882687e+003;
n=37; farx(n+1)=1.337016e+001; foe(n+1)=1.416006e+002; krok(n+1)=2.614439e-003; ng(n+1)=7.537501e+003;
n=38; farx(n+1)=1.402912e+001; foe(n+1)=1.388397e+002; krok(n+1)=2.774090e-003; ng(n+1)=8.538334e+003;
n=39; farx(n+1)=1.162273e+001; foe(n+1)=1.303579e+002; krok(n+1)=6.480489e-003; ng(n+1)=8.766803e+003;
n=40; farx(n+1)=1.056988e+001; foe(n+1)=1.264426e+002; krok(n+1)=1.688016e-003; ng(n+1)=1.265008e+004;
n=41; farx(n+1)=1.040849e+001; foe(n+1)=1.217619e+002; krok(n+1)=1.488820e-002; ng(n+1)=1.268221e+004;
n=42; farx(n+1)=1.002129e+001; foe(n+1)=1.189966e+002; krok(n+1)=1.340015e-002; ng(n+1)=1.081413e+004;
n=43; farx(n+1)=8.623140e+000; foe(n+1)=1.084382e+002; krok(n+1)=1.238474e-002; ng(n+1)=9.584931e+003;
n=44; farx(n+1)=7.329093e+000; foe(n+1)=8.437757e+001; krok(n+1)=1.011873e-001; ng(n+1)=9.463812e+003;
n=45; farx(n+1)=5.820429e+000; foe(n+1)=6.264635e+001; krok(n+1)=1.468011e-001; ng(n+1)=6.474701e+003;
n=46; farx(n+1)=5.734238e+000; foe(n+1)=6.142643e+001; krok(n+1)=4.547265e-001; ng(n+1)=2.864794e+003;
n=47; farx(n+1)=5.909743e+000; foe(n+1)=5.995772e+001; krok(n+1)=5.657309e-001; ng(n+1)=1.517774e+003;
n=48; farx(n+1)=5.690180e+000; foe(n+1)=5.908768e+001; krok(n+1)=7.178655e-001; ng(n+1)=1.457728e+003;
n=49; farx(n+1)=5.570263e+000; foe(n+1)=5.850912e+001; krok(n+1)=1.104103e+000; ng(n+1)=9.690673e+002;
n=50; farx(n+1)=5.595421e+000; foe(n+1)=5.823045e+001; krok(n+1)=1.814761e+000; ng(n+1)=6.597760e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=5.585718e+000; foe(n+1)=5.820614e+001; krok(n+1)=3.323679e-007; ng(n+1)=6.606085e+002;
n=52; farx(n+1)=5.575868e+000; foe(n+1)=5.818854e+001; krok(n+1)=3.787925e-006; ng(n+1)=1.577833e+002;
n=53; farx(n+1)=5.248482e+000; foe(n+1)=5.784948e+001; krok(n+1)=2.979192e-004; ng(n+1)=1.199200e+002;
n=54; farx(n+1)=5.168084e+000; foe(n+1)=5.779542e+001; krok(n+1)=5.538211e-005; ng(n+1)=4.481472e+002;
n=55; farx(n+1)=5.997153e+000; foe(n+1)=5.447690e+001; krok(n+1)=1.088038e-003; ng(n+1)=4.984684e+002;
n=56; farx(n+1)=6.088911e+000; foe(n+1)=5.427705e+001; krok(n+1)=1.528320e-003; ng(n+1)=1.315257e+003;
n=57; farx(n+1)=6.084858e+000; foe(n+1)=5.418010e+001; krok(n+1)=2.288181e-003; ng(n+1)=9.013707e+002;
n=58; farx(n+1)=6.178007e+000; foe(n+1)=5.381170e+001; krok(n+1)=3.067118e-002; ng(n+1)=8.832746e+002;
n=59; farx(n+1)=6.159180e+000; foe(n+1)=5.353956e+001; krok(n+1)=3.778589e-003; ng(n+1)=1.758206e+003;
n=60; farx(n+1)=6.288548e+000; foe(n+1)=5.264142e+001; krok(n+1)=1.605804e-002; ng(n+1)=1.485073e+003;
n=61; farx(n+1)=6.429923e+000; foe(n+1)=5.123783e+001; krok(n+1)=1.296098e-002; ng(n+1)=3.564421e+003;
n=62; farx(n+1)=6.491046e+000; foe(n+1)=5.104230e+001; krok(n+1)=4.312903e-003; ng(n+1)=5.126005e+003;
n=63; farx(n+1)=6.547920e+000; foe(n+1)=5.071323e+001; krok(n+1)=1.121678e-002; ng(n+1)=4.999797e+003;
n=64; farx(n+1)=6.334088e+000; foe(n+1)=5.035235e+001; krok(n+1)=6.774452e-002; ng(n+1)=5.230773e+003;
n=65; farx(n+1)=6.865068e+000; foe(n+1)=4.893336e+001; krok(n+1)=6.392790e-002; ng(n+1)=4.865991e+003;
n=66; farx(n+1)=7.003819e+000; foe(n+1)=4.843584e+001; krok(n+1)=4.895097e-002; ng(n+1)=1.344849e+003;
n=67; farx(n+1)=7.072708e+000; foe(n+1)=4.759387e+001; krok(n+1)=1.096649e-001; ng(n+1)=1.816477e+003;
n=68; farx(n+1)=7.780482e+000; foe(n+1)=4.593089e+001; krok(n+1)=1.261335e-001; ng(n+1)=2.931956e+003;
n=69; farx(n+1)=6.115942e+000; foe(n+1)=4.233561e+001; krok(n+1)=7.807902e-001; ng(n+1)=2.174410e+003;
n=70; farx(n+1)=4.405520e+000; foe(n+1)=3.775079e+001; krok(n+1)=2.349459e-001; ng(n+1)=2.935186e+003;
n=71; farx(n+1)=4.264530e+000; foe(n+1)=3.547991e+001; krok(n+1)=3.556863e-001; ng(n+1)=2.739110e+003;
n=72; farx(n+1)=4.055061e+000; foe(n+1)=3.332737e+001; krok(n+1)=4.406062e-001; ng(n+1)=2.354070e+003;
n=73; farx(n+1)=3.969740e+000; foe(n+1)=3.191822e+001; krok(n+1)=1.666377e-001; ng(n+1)=1.584877e+003;
n=74; farx(n+1)=4.089994e+000; foe(n+1)=2.920520e+001; krok(n+1)=1.300904e-001; ng(n+1)=5.409349e+002;
n=75; farx(n+1)=3.710207e+000; foe(n+1)=2.675508e+001; krok(n+1)=3.963116e-001; ng(n+1)=2.106040e+003;
%odnowa zmiennej metryki
n=76; farx(n+1)=3.607619e+000; foe(n+1)=2.583594e+001; krok(n+1)=7.343206e-007; ng(n+1)=2.659227e+003;
n=77; farx(n+1)=3.607733e+000; foe(n+1)=2.549571e+001; krok(n+1)=2.684775e-006; ng(n+1)=9.887784e+002;
n=78; farx(n+1)=3.517978e+000; foe(n+1)=2.514714e+001; krok(n+1)=5.807097e-006; ng(n+1)=6.738561e+002;
n=79; farx(n+1)=3.139782e+000; foe(n+1)=2.282015e+001; krok(n+1)=1.210454e-003; ng(n+1)=1.734188e+002;
n=80; farx(n+1)=2.548945e+000; foe(n+1)=2.073421e+001; krok(n+1)=9.701682e-004; ng(n+1)=5.523361e+002;
n=81; farx(n+1)=2.089428e+000; foe(n+1)=1.854961e+001; krok(n+1)=1.737584e-003; ng(n+1)=7.937982e+002;
n=82; farx(n+1)=2.110326e+000; foe(n+1)=1.813977e+001; krok(n+1)=9.143281e-004; ng(n+1)=2.080536e+003;
n=83; farx(n+1)=2.225686e+000; foe(n+1)=1.623366e+001; krok(n+1)=1.054870e-002; ng(n+1)=1.735184e+003;
n=84; farx(n+1)=2.064103e+000; foe(n+1)=1.426444e+001; krok(n+1)=7.314625e-003; ng(n+1)=3.261655e+003;
n=85; farx(n+1)=2.252782e+000; foe(n+1)=1.312436e+001; krok(n+1)=6.446985e-003; ng(n+1)=1.766479e+003;
n=86; farx(n+1)=2.131048e+000; foe(n+1)=1.216018e+001; krok(n+1)=8.029022e-003; ng(n+1)=1.982983e+003;
n=87; farx(n+1)=2.055278e+000; foe(n+1)=1.212465e+001; krok(n+1)=3.096185e-003; ng(n+1)=3.914678e+002;
n=88; farx(n+1)=1.967700e+000; foe(n+1)=1.208639e+001; krok(n+1)=4.761479e-003; ng(n+1)=2.615599e+002;
n=89; farx(n+1)=1.807732e+000; foe(n+1)=1.172979e+001; krok(n+1)=5.576355e-002; ng(n+1)=5.494877e+002;
n=90; farx(n+1)=1.497353e+000; foe(n+1)=1.111507e+001; krok(n+1)=1.142218e-002; ng(n+1)=1.661038e+002;
n=91; farx(n+1)=1.555515e+000; foe(n+1)=9.937888e+000; krok(n+1)=1.076361e-001; ng(n+1)=1.840437e+003;
n=92; farx(n+1)=1.872867e+000; foe(n+1)=8.381110e+000; krok(n+1)=9.205345e-001; ng(n+1)=6.471071e+002;
n=93; farx(n+1)=1.919378e+000; foe(n+1)=7.722572e+000; krok(n+1)=9.632832e-001; ng(n+1)=6.106569e+002;
n=94; farx(n+1)=2.070471e+000; foe(n+1)=7.496668e+000; krok(n+1)=7.382907e-001; ng(n+1)=2.984302e+002;
n=95; farx(n+1)=2.056050e+000; foe(n+1)=7.304042e+000; krok(n+1)=8.008859e-001; ng(n+1)=3.356531e+002;
n=96; farx(n+1)=1.969242e+000; foe(n+1)=7.141087e+000; krok(n+1)=4.763389e-001; ng(n+1)=1.004344e+003;
n=97; farx(n+1)=1.687458e+000; foe(n+1)=7.019505e+000; krok(n+1)=6.614848e-001; ng(n+1)=2.781451e+002;
n=98; farx(n+1)=1.631229e+000; foe(n+1)=6.868810e+000; krok(n+1)=8.055989e-001; ng(n+1)=3.215432e+002;
n=99; farx(n+1)=1.546917e+000; foe(n+1)=6.759343e+000; krok(n+1)=5.350602e-001; ng(n+1)=1.032137e+003;
n=100; farx(n+1)=1.584092e+000; foe(n+1)=6.693880e+000; krok(n+1)=4.536903e-001; ng(n+1)=5.320354e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
