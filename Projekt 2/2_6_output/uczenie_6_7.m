%uczenie predyktora arx
clear all;
n=0; farx(n+1)=4.323627e+003; foe(n+1)=4.235298e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.194607e+003; foe(n+1)=3.114106e+003; krok(n+1)=4.216476e-004; ng(n+1)=2.982120e+003;
n=2; farx(n+1)=3.044537e+003; foe(n+1)=3.045592e+003; krok(n+1)=2.523757e-004; ng(n+1)=1.671833e+003;
n=3; farx(n+1)=1.104747e+003; foe(n+1)=7.217208e+003; krok(n+1)=1.669364e-002; ng(n+1)=2.603617e+002;
n=4; farx(n+1)=1.092987e+003; foe(n+1)=7.722653e+003; krok(n+1)=1.729605e-005; ng(n+1)=2.190369e+003;
n=5; farx(n+1)=1.077527e+003; foe(n+1)=8.080297e+003; krok(n+1)=1.949868e-002; ng(n+1)=2.199821e+003;
n=6; farx(n+1)=4.206428e+002; foe(n+1)=7.886474e+003; krok(n+1)=2.240289e+000; ng(n+1)=2.390724e+003;
n=7; farx(n+1)=2.272866e+002; foe(n+1)=6.644324e+003; krok(n+1)=8.316025e-001; ng(n+1)=7.177957e+002;
n=8; farx(n+1)=1.325916e+002; foe(n+1)=6.861807e+003; krok(n+1)=3.108763e+000; ng(n+1)=5.040789e+002;
n=9; farx(n+1)=6.615472e+001; foe(n+1)=1.103436e+004; krok(n+1)=1.445753e+000; ng(n+1)=3.777809e+002;
n=10; farx(n+1)=5.832859e+001; foe(n+1)=1.088114e+004; krok(n+1)=1.009068e+000; ng(n+1)=1.772838e+002;
n=11; farx(n+1)=5.631373e+001; foe(n+1)=9.173707e+003; krok(n+1)=9.073805e-001; ng(n+1)=5.743710e+001;
n=12; farx(n+1)=5.116719e+001; foe(n+1)=5.245518e+003; krok(n+1)=3.982210e+000; ng(n+1)=8.332996e+001;
n=13; farx(n+1)=4.764786e+001; foe(n+1)=2.941725e+002; krok(n+1)=9.769152e-001; ng(n+1)=1.161710e+002;
n=14; farx(n+1)=4.667321e+001; foe(n+1)=1.012633e+003; krok(n+1)=7.382907e-001; ng(n+1)=9.073568e+001;
n=15; farx(n+1)=4.428436e+001; foe(n+1)=2.234375e+002; krok(n+1)=4.650420e+000; ng(n+1)=7.044440e+001;
n=16; farx(n+1)=3.328426e+001; foe(n+1)=1.685775e+002; krok(n+1)=3.944619e+000; ng(n+1)=1.186197e+002;
n=17; farx(n+1)=2.930880e+001; foe(n+1)=7.720290e+002; krok(n+1)=4.536903e-001; ng(n+1)=1.577481e+002;
n=18; farx(n+1)=2.508931e+001; foe(n+1)=1.086169e+003; krok(n+1)=8.892951e-001; ng(n+1)=1.774243e+002;
n=19; farx(n+1)=2.349986e+001; foe(n+1)=3.226735e+002; krok(n+1)=4.641165e-001; ng(n+1)=3.173257e+001;
n=20; farx(n+1)=2.282399e+001; foe(n+1)=2.516470e+002; krok(n+1)=4.268310e-001; ng(n+1)=1.072178e+002;
n=21; farx(n+1)=2.233539e+001; foe(n+1)=7.520847e+002; krok(n+1)=1.242513e+000; ng(n+1)=3.696864e+001;
n=22; farx(n+1)=2.205444e+001; foe(n+1)=7.582129e+002; krok(n+1)=7.487569e-001; ng(n+1)=2.851566e+001;
n=23; farx(n+1)=2.197649e+001; foe(n+1)=3.749520e+002; krok(n+1)=7.610808e-001; ng(n+1)=9.765991e+000;
n=24; farx(n+1)=2.193795e+001; foe(n+1)=4.289848e+002; krok(n+1)=1.462039e+000; ng(n+1)=2.391962e+001;
n=25; farx(n+1)=2.191411e+001; foe(n+1)=4.163089e+002; krok(n+1)=1.902505e+000; ng(n+1)=8.017406e+000;
%odnowa zmiennej metryki
n=26; farx(n+1)=2.191200e+001; foe(n+1)=4.230382e+002; krok(n+1)=4.031350e-004; ng(n+1)=5.826855e+000;
n=27; farx(n+1)=2.191111e+001; foe(n+1)=4.024856e+002; krok(n+1)=3.534472e-004; ng(n+1)=3.786884e+000;
n=28; farx(n+1)=2.190531e+001; foe(n+1)=3.464352e+002; krok(n+1)=2.659087e-003; ng(n+1)=3.258322e+000;
n=29; farx(n+1)=2.190072e+001; foe(n+1)=3.409147e+002; krok(n+1)=1.142218e-002; ng(n+1)=1.509540e+000;
n=30; farx(n+1)=2.189886e+001; foe(n+1)=3.480183e+002; krok(n+1)=1.109636e-002; ng(n+1)=9.201110e-001;
n=31; farx(n+1)=2.189286e+001; foe(n+1)=3.177661e+002; krok(n+1)=7.238598e-001; ng(n+1)=7.584464e-001;
n=32; farx(n+1)=7.961528e+000; foe(n+1)=2.272225e+002; krok(n+1)=6.474872e+001; ng(n+1)=4.467029e-001;
n=33; farx(n+1)=7.955854e+000; foe(n+1)=2.266000e+002; krok(n+1)=7.481826e-003; ng(n+1)=6.411923e+001;
n=34; farx(n+1)=6.110860e+000; foe(n+1)=1.830143e+002; krok(n+1)=2.658374e-001; ng(n+1)=6.368464e+001;
n=35; farx(n+1)=5.483779e+000; foe(n+1)=2.220986e+002; krok(n+1)=1.454604e+000; ng(n+1)=1.008212e+002;
n=36; farx(n+1)=4.805508e+000; foe(n+1)=1.674267e+002; krok(n+1)=8.008859e-001; ng(n+1)=9.298507e+001;
n=37; farx(n+1)=4.509488e+000; foe(n+1)=1.408443e+002; krok(n+1)=3.349286e-001; ng(n+1)=5.192593e+001;
n=38; farx(n+1)=4.220342e+000; foe(n+1)=1.444519e+002; krok(n+1)=6.900158e-001; ng(n+1)=7.342329e+001;
n=39; farx(n+1)=3.947673e+000; foe(n+1)=1.615671e+002; krok(n+1)=1.242513e+000; ng(n+1)=8.375986e+001;
n=40; farx(n+1)=3.632094e+000; foe(n+1)=1.474123e+002; krok(n+1)=1.585247e+000; ng(n+1)=6.134983e+001;
n=41; farx(n+1)=3.431479e+000; foe(n+1)=1.249042e+002; krok(n+1)=2.876382e-001; ng(n+1)=6.937866e+001;
n=42; farx(n+1)=3.334953e+000; foe(n+1)=1.295279e+002; krok(n+1)=4.884576e-001; ng(n+1)=5.372575e+001;
n=43; farx(n+1)=3.247667e+000; foe(n+1)=1.553170e+002; krok(n+1)=1.045307e+000; ng(n+1)=8.118519e+001;
n=44; farx(n+1)=3.129240e+000; foe(n+1)=1.327950e+002; krok(n+1)=1.247213e+000; ng(n+1)=2.400388e+001;
n=45; farx(n+1)=3.025564e+000; foe(n+1)=1.242331e+002; krok(n+1)=6.167589e-001; ng(n+1)=4.597762e+001;
n=46; farx(n+1)=2.966220e+000; foe(n+1)=1.255173e+002; krok(n+1)=5.945315e-001; ng(n+1)=7.914322e+001;
n=47; farx(n+1)=2.878972e+000; foe(n+1)=1.280433e+002; krok(n+1)=8.377184e-001; ng(n+1)=1.949098e+001;
n=48; farx(n+1)=2.845823e+000; foe(n+1)=1.358238e+002; krok(n+1)=7.461831e-001; ng(n+1)=3.260195e+001;
n=49; farx(n+1)=2.784863e+000; foe(n+1)=1.268102e+002; krok(n+1)=1.893841e+000; ng(n+1)=2.589427e+001;
n=50; farx(n+1)=2.764841e+000; foe(n+1)=1.294391e+002; krok(n+1)=4.004429e-001; ng(n+1)=6.405755e+001;
%odnowa zmiennej metryki
n=51; farx(n+1)=2.761458e+000; foe(n+1)=1.311744e+002; krok(n+1)=2.769106e-005; ng(n+1)=2.360220e+001;
n=52; farx(n+1)=2.758797e+000; foe(n+1)=1.297095e+002; krok(n+1)=1.042068e-004; ng(n+1)=1.246741e+001;
n=53; farx(n+1)=2.758045e+000; foe(n+1)=1.268815e+002; krok(n+1)=7.692316e-003; ng(n+1)=6.502573e-001;
n=54; farx(n+1)=2.742247e+000; foe(n+1)=1.309519e+002; krok(n+1)=2.160874e-002; ng(n+1)=1.853849e+000;
n=55; farx(n+1)=2.729005e+000; foe(n+1)=1.289379e+002; krok(n+1)=4.858193e-002; ng(n+1)=7.552377e+000;
n=56; farx(n+1)=2.724172e+000; foe(n+1)=1.226889e+002; krok(n+1)=8.005809e-002; ng(n+1)=2.825583e+001;
n=57; farx(n+1)=2.718129e+000; foe(n+1)=1.200550e+002; krok(n+1)=2.621075e-001; ng(n+1)=3.966330e+001;
n=58; farx(n+1)=2.700899e+000; foe(n+1)=1.221999e+002; krok(n+1)=3.332754e-001; ng(n+1)=3.612237e+001;
n=59; farx(n+1)=2.680938e+000; foe(n+1)=1.206550e+002; krok(n+1)=1.386016e+000; ng(n+1)=2.340632e+001;
n=60; farx(n+1)=2.670944e+000; foe(n+1)=1.226786e+002; krok(n+1)=8.008859e-001; ng(n+1)=3.001254e+001;
n=61; farx(n+1)=2.661747e+000; foe(n+1)=1.159615e+002; krok(n+1)=7.101669e-001; ng(n+1)=1.365434e+001;
n=62; farx(n+1)=2.652728e+000; foe(n+1)=1.150546e+002; krok(n+1)=1.307567e+000; ng(n+1)=2.171329e+001;
n=63; farx(n+1)=2.649568e+000; foe(n+1)=1.177602e+002; krok(n+1)=6.338619e-001; ng(n+1)=2.518684e+001;
n=64; farx(n+1)=2.642934e+000; foe(n+1)=1.183350e+002; krok(n+1)=9.574267e-001; ng(n+1)=7.951874e+000;
n=65; farx(n+1)=2.640676e+000; foe(n+1)=1.155710e+002; krok(n+1)=3.839848e-001; ng(n+1)=8.760184e+000;
n=66; farx(n+1)=2.637333e+000; foe(n+1)=1.136068e+002; krok(n+1)=1.352447e+000; ng(n+1)=1.828416e+001;
n=67; farx(n+1)=2.634436e+000; foe(n+1)=1.116543e+002; krok(n+1)=8.413217e-001; ng(n+1)=1.040875e+001;
n=68; farx(n+1)=2.632622e+000; foe(n+1)=1.157151e+002; krok(n+1)=1.065517e+000; ng(n+1)=1.413765e+001;
n=69; farx(n+1)=2.631320e+000; foe(n+1)=1.140709e+002; krok(n+1)=1.280117e+000; ng(n+1)=4.332710e+000;
n=70; farx(n+1)=2.630392e+000; foe(n+1)=1.129592e+002; krok(n+1)=6.614848e-001; ng(n+1)=1.244743e+001;
n=71; farx(n+1)=2.629802e+000; foe(n+1)=1.111087e+002; krok(n+1)=7.655551e-001; ng(n+1)=6.587291e+000;
n=72; farx(n+1)=2.629257e+000; foe(n+1)=1.128751e+002; krok(n+1)=1.569476e+000; ng(n+1)=5.314198e+000;
n=73; farx(n+1)=2.629082e+000; foe(n+1)=1.131514e+002; krok(n+1)=9.397835e-001; ng(n+1)=4.652451e+000;
n=74; farx(n+1)=2.628919e+000; foe(n+1)=1.119787e+002; krok(n+1)=7.847380e-001; ng(n+1)=1.241018e+000;
n=75; farx(n+1)=2.628860e+000; foe(n+1)=1.116562e+002; krok(n+1)=1.200452e+000; ng(n+1)=3.394337e+000;
%odnowa zmiennej metryki
n=76; farx(n+1)=2.628850e+000; foe(n+1)=1.117489e+002; krok(n+1)=1.171829e-005; ng(n+1)=1.777768e+000;
n=77; farx(n+1)=2.628847e+000; foe(n+1)=1.118138e+002; krok(n+1)=8.233132e-005; ng(n+1)=2.791480e-001;
n=78; farx(n+1)=2.628845e+000; foe(n+1)=1.119975e+002; krok(n+1)=5.352819e-004; ng(n+1)=1.250694e-001;
n=79; farx(n+1)=2.628817e+000; foe(n+1)=1.119598e+002; krok(n+1)=5.695205e-002; ng(n+1)=4.292040e-002;
n=80; farx(n+1)=2.628797e+000; foe(n+1)=1.119235e+002; krok(n+1)=3.805310e-002; ng(n+1)=4.648575e-002;
n=81; farx(n+1)=2.628787e+000; foe(n+1)=1.118285e+002; krok(n+1)=6.045742e-002; ng(n+1)=2.585008e-002;
n=82; farx(n+1)=2.628781e+000; foe(n+1)=1.117441e+002; krok(n+1)=8.137925e+000; ng(n+1)=2.194907e-002;
n=83; farx(n+1)=2.628781e+000; foe(n+1)=1.117338e+002; krok(n+1)=1.009068e+000; ng(n+1)=1.219983e-001;
n=84; farx(n+1)=2.628781e+000; foe(n+1)=1.117340e+002; krok(n+1)=9.834271e-001; ng(n+1)=3.582110e-003;
n=85; farx(n+1)=2.628781e+000; foe(n+1)=1.117340e+002; krok(n+1)=3.620407e-006; ng(n+1)=3.915756e-004;
n=86; farx(n+1)=2.628781e+000; foe(n+1)=1.117340e+002; krok(n+1)=5.586175e-008; ng(n+1)=3.915741e-004;
 % z�y kierunek w metodzie zm - odnowa 
n=87; farx(n+1)=2.628781e+000; foe(n+1)=1.117340e+002; krok(n+1)=1.203661e-005; ng(n+1)=3.915741e-004;
n=88; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=6.836820e-005; ng(n+1)=6.937432e-005;
n=89; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=3.575152e-006; ng(n+1)=1.618682e-005;
n=90; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=2.190752e-005; ng(n+1)=1.618912e-005;
n=91; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.412899e-006; ng(n+1)=1.614656e-005;
n=92; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=4.283100e-006; ng(n+1)=1.614651e-005;
n=93; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=8.040168e-008; ng(n+1)=1.614639e-005;
n=94; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=2.140943e-008; ng(n+1)=1.614639e-005;
n=95; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=5.996128e-008; ng(n+1)=1.614639e-005;
n=96; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.639991e-005; ng(n+1)=1.614639e-005;
n=97; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=9.952471e-006; ng(n+1)=1.614612e-005;
n=98; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=2.681339e-009; ng(n+1)=1.614595e-005;
n=99; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=7.086453e-009; ng(n+1)=1.614595e-005;
n=100; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.523371e-006; ng(n+1)=1.614595e-005;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora ARX');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
