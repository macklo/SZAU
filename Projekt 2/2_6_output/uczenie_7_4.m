%uczenie predyktora arx
clear all;
n=0; farx(n+1)=4.699213e+003; foe(n+1)=4.904655e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.222836e+003; foe(n+1)=4.080967e+003; krok(n+1)=4.233525e-004; ng(n+1)=9.559610e+003;
n=2; farx(n+1)=9.656271e+002; foe(n+1)=1.224006e+004; krok(n+1)=1.582138e-004; ng(n+1)=1.410208e+004;
n=3; farx(n+1)=3.167117e+002; foe(n+1)=7.718392e+003; krok(n+1)=2.556540e-004; ng(n+1)=6.899378e+003;
n=4; farx(n+1)=1.550110e+002; foe(n+1)=1.479747e+004; krok(n+1)=3.660956e-003; ng(n+1)=5.182368e+003;
n=5; farx(n+1)=1.075137e+002; foe(n+1)=1.492047e+004; krok(n+1)=4.651747e-004; ng(n+1)=3.538211e+003;
n=6; farx(n+1)=7.723355e+001; foe(n+1)=2.196838e+004; krok(n+1)=1.166078e-003; ng(n+1)=1.281428e+003;
n=7; farx(n+1)=6.394005e+001; foe(n+1)=1.465132e+004; krok(n+1)=8.410247e-004; ng(n+1)=2.365319e+003;
n=8; farx(n+1)=3.713572e+001; foe(n+1)=8.457099e+003; krok(n+1)=4.782427e-003; ng(n+1)=2.091894e+003;
n=9; farx(n+1)=8.416732e+000; foe(n+1)=2.481396e+003; krok(n+1)=3.937497e-003; ng(n+1)=2.138025e+003;
n=10; farx(n+1)=5.786168e+000; foe(n+1)=7.765409e+002; krok(n+1)=1.111857e-003; ng(n+1)=4.583022e+002;
n=11; farx(n+1)=3.089128e+000; foe(n+1)=1.855486e+002; krok(n+1)=7.029269e-003; ng(n+1)=4.313472e+002;
n=12; farx(n+1)=2.646664e+000; foe(n+1)=1.868128e+002; krok(n+1)=8.349282e-003; ng(n+1)=1.195581e+002;
n=13; farx(n+1)=2.047647e+000; foe(n+1)=9.609911e+001; krok(n+1)=2.398987e-002; ng(n+1)=1.499398e+002;
n=14; farx(n+1)=1.821321e+000; foe(n+1)=8.421423e+001; krok(n+1)=2.721013e-002; ng(n+1)=1.110010e+002;
n=15; farx(n+1)=1.712905e+000; foe(n+1)=9.528886e+001; krok(n+1)=2.741202e-002; ng(n+1)=8.262428e+001;
n=16; farx(n+1)=1.429047e+000; foe(n+1)=1.116135e+002; krok(n+1)=1.467726e-001; ng(n+1)=6.680286e+001;
n=17; farx(n+1)=1.343599e+000; foe(n+1)=9.601689e+001; krok(n+1)=5.793922e-002; ng(n+1)=7.919153e+001;
n=18; farx(n+1)=1.174725e+000; foe(n+1)=9.835423e+001; krok(n+1)=1.580662e-001; ng(n+1)=5.935361e+001;
n=19; farx(n+1)=1.120215e+000; foe(n+1)=1.224572e+002; krok(n+1)=1.019950e-001; ng(n+1)=4.182709e+001;
n=20; farx(n+1)=1.034023e+000; foe(n+1)=1.172899e+002; krok(n+1)=1.827549e-001; ng(n+1)=3.824098e+001;
n=21; farx(n+1)=8.920273e-001; foe(n+1)=1.264877e+002; krok(n+1)=2.775080e-001; ng(n+1)=2.476007e+001;
n=22; farx(n+1)=8.756877e-001; foe(n+1)=1.306425e+002; krok(n+1)=4.852216e-002; ng(n+1)=3.447602e+001;
n=23; farx(n+1)=8.517571e-001; foe(n+1)=1.657119e+002; krok(n+1)=7.356694e-002; ng(n+1)=1.468851e+001;
n=24; farx(n+1)=8.261369e-001; foe(n+1)=1.534830e+002; krok(n+1)=2.286351e-001; ng(n+1)=1.612564e+001;
n=25; farx(n+1)=7.485881e-001; foe(n+1)=1.710389e+002; krok(n+1)=4.602673e-001; ng(n+1)=4.429250e+001;
%odnowa zmiennej metryki
n=26; farx(n+1)=7.353731e-001; foe(n+1)=1.529110e+002; krok(n+1)=1.307587e-004; ng(n+1)=5.167625e+001;
n=27; farx(n+1)=7.337518e-001; foe(n+1)=1.665397e+002; krok(n+1)=1.112345e-004; ng(n+1)=1.722013e+001;
n=28; farx(n+1)=7.224180e-001; foe(n+1)=1.555060e+002; krok(n+1)=3.364099e-003; ng(n+1)=8.952249e+000;
n=29; farx(n+1)=7.019274e-001; foe(n+1)=1.610520e+002; krok(n+1)=1.631060e-002; ng(n+1)=5.802868e+000;
n=30; farx(n+1)=6.715673e-001; foe(n+1)=2.160800e+002; krok(n+1)=2.384691e-002; ng(n+1)=8.270606e+000;
n=31; farx(n+1)=6.054611e-001; foe(n+1)=1.052777e+002; krok(n+1)=4.438543e-002; ng(n+1)=3.328071e+001;
n=32; farx(n+1)=5.736288e-001; foe(n+1)=8.397161e+001; krok(n+1)=3.022871e-002; ng(n+1)=3.346368e+001;
n=33; farx(n+1)=5.564496e-001; foe(n+1)=7.086406e+001; krok(n+1)=2.803148e-002; ng(n+1)=5.409134e+001;
n=34; farx(n+1)=5.397710e-001; foe(n+1)=6.646102e+001; krok(n+1)=4.614317e-002; ng(n+1)=4.830373e+001;
n=35; farx(n+1)=5.106961e-001; foe(n+1)=8.387835e+001; krok(n+1)=1.134226e-001; ng(n+1)=2.214864e+001;
n=36; farx(n+1)=4.664173e-001; foe(n+1)=5.355693e+001; krok(n+1)=2.152723e-001; ng(n+1)=4.999906e+001;
n=37; farx(n+1)=4.463239e-001; foe(n+1)=4.857784e+001; krok(n+1)=1.487872e-001; ng(n+1)=2.408038e+001;
n=38; farx(n+1)=4.379147e-001; foe(n+1)=5.122064e+001; krok(n+1)=3.103068e-002; ng(n+1)=3.181273e+001;
n=39; farx(n+1)=4.280157e-001; foe(n+1)=4.393119e+001; krok(n+1)=1.728699e-001; ng(n+1)=1.045639e+001;
n=40; farx(n+1)=4.227462e-001; foe(n+1)=4.200296e+001; krok(n+1)=7.317146e-002; ng(n+1)=2.108959e+001;
n=41; farx(n+1)=4.162062e-001; foe(n+1)=4.019422e+001; krok(n+1)=2.601807e-001; ng(n+1)=4.507587e+000;
n=42; farx(n+1)=4.107534e-001; foe(n+1)=4.437126e+001; krok(n+1)=1.277097e-001; ng(n+1)=1.275710e+001;
n=43; farx(n+1)=4.042206e-001; foe(n+1)=3.527925e+001; krok(n+1)=1.176675e-001; ng(n+1)=3.040054e+001;
n=44; farx(n+1)=3.952908e-001; foe(n+1)=3.358683e+001; krok(n+1)=1.338969e-001; ng(n+1)=4.097558e+001;
n=45; farx(n+1)=3.905354e-001; foe(n+1)=3.092654e+001; krok(n+1)=3.482091e-002; ng(n+1)=2.160693e+001;
n=46; farx(n+1)=3.780026e-001; foe(n+1)=2.249011e+001; krok(n+1)=1.871606e-001; ng(n+1)=1.154574e+001;
n=47; farx(n+1)=3.689558e-001; foe(n+1)=2.125078e+001; krok(n+1)=1.134226e-001; ng(n+1)=3.070248e+001;
n=48; farx(n+1)=3.605543e-001; foe(n+1)=1.916245e+001; krok(n+1)=4.128102e-001; ng(n+1)=5.336793e+000;
n=49; farx(n+1)=3.555066e-001; foe(n+1)=1.809155e+001; krok(n+1)=2.745494e-001; ng(n+1)=9.467329e+000;
n=50; farx(n+1)=3.439952e-001; foe(n+1)=1.415476e+001; krok(n+1)=5.316748e-001; ng(n+1)=2.427294e+001;
%odnowa zmiennej metryki
n=51; farx(n+1)=3.410985e-001; foe(n+1)=1.335294e+001; krok(n+1)=7.813213e-005; ng(n+1)=2.584300e+001;
n=52; farx(n+1)=3.406635e-001; foe(n+1)=1.321294e+001; krok(n+1)=9.242591e-005; ng(n+1)=1.073790e+001;
n=53; farx(n+1)=3.390363e-001; foe(n+1)=1.383812e+001; krok(n+1)=1.109636e-002; ng(n+1)=1.612991e+000;
n=54; farx(n+1)=3.379390e-001; foe(n+1)=1.365036e+001; krok(n+1)=4.696827e-003; ng(n+1)=2.374114e+000;
n=55; farx(n+1)=3.354782e-001; foe(n+1)=1.456059e+001; krok(n+1)=9.659422e-003; ng(n+1)=2.522139e+000;
n=56; farx(n+1)=3.332221e-001; foe(n+1)=1.371294e+001; krok(n+1)=1.031518e-001; ng(n+1)=3.067275e+000;
n=57; farx(n+1)=3.312383e-001; foe(n+1)=1.254727e+001; krok(n+1)=1.297775e-001; ng(n+1)=5.958680e+000;
n=58; farx(n+1)=3.288967e-001; foe(n+1)=1.197720e+001; krok(n+1)=5.857529e-002; ng(n+1)=1.115250e+001;
n=59; farx(n+1)=3.266300e-001; foe(n+1)=1.055759e+001; krok(n+1)=9.599621e-002; ng(n+1)=2.237220e+001;
n=60; farx(n+1)=3.231395e-001; foe(n+1)=9.245918e+000; krok(n+1)=4.063898e-001; ng(n+1)=1.948195e+001;
n=61; farx(n+1)=3.204519e-001; foe(n+1)=7.784124e+000; krok(n+1)=1.238533e-001; ng(n+1)=1.218282e+001;
n=62; farx(n+1)=3.141019e-001; foe(n+1)=7.431221e+000; krok(n+1)=6.428475e-001; ng(n+1)=7.007722e+000;
n=63; farx(n+1)=3.121322e-001; foe(n+1)=7.370601e+000; krok(n+1)=8.268560e-002; ng(n+1)=2.160540e+001;
n=64; farx(n+1)=3.083618e-001; foe(n+1)=7.232279e+000; krok(n+1)=2.223238e-001; ng(n+1)=1.216326e+001;
n=65; farx(n+1)=3.054298e-001; foe(n+1)=7.289024e+000; krok(n+1)=2.203031e-001; ng(n+1)=1.062324e+001;
n=66; farx(n+1)=3.020478e-001; foe(n+1)=7.147968e+000; krok(n+1)=1.919924e-001; ng(n+1)=2.095020e+001;
n=67; farx(n+1)=2.974623e-001; foe(n+1)=6.751535e+000; krok(n+1)=3.240219e-001; ng(n+1)=9.041102e+000;
n=68; farx(n+1)=2.945805e-001; foe(n+1)=6.936621e+000; krok(n+1)=1.042621e-001; ng(n+1)=1.803818e+001;
n=69; farx(n+1)=2.927469e-001; foe(n+1)=6.850491e+000; krok(n+1)=1.221144e-001; ng(n+1)=7.412466e+000;
n=70; farx(n+1)=2.895269e-001; foe(n+1)=6.243938e+000; krok(n+1)=4.461224e-001; ng(n+1)=7.974169e+000;
n=71; farx(n+1)=2.870931e-001; foe(n+1)=6.277264e+000; krok(n+1)=2.846355e-001; ng(n+1)=6.933669e+000;
n=72; farx(n+1)=2.840554e-001; foe(n+1)=5.822634e+000; krok(n+1)=4.063898e-001; ng(n+1)=8.854590e+000;
n=73; farx(n+1)=2.825261e-001; foe(n+1)=5.199864e+000; krok(n+1)=4.004429e-001; ng(n+1)=2.432204e+000;
n=74; farx(n+1)=2.803651e-001; foe(n+1)=4.749675e+000; krok(n+1)=8.838889e-001; ng(n+1)=2.570485e+000;
n=75; farx(n+1)=2.790596e-001; foe(n+1)=4.322998e+000; krok(n+1)=3.076393e-001; ng(n+1)=7.413552e+000;
%odnowa zmiennej metryki
n=76; farx(n+1)=2.783661e-001; foe(n+1)=4.246338e+000; krok(n+1)=7.131570e-005; ng(n+1)=1.576301e+001;
n=77; farx(n+1)=2.783333e-001; foe(n+1)=4.239462e+000; krok(n+1)=1.419195e-004; ng(n+1)=2.343834e+000;
n=78; farx(n+1)=2.782095e-001; foe(n+1)=4.306045e+000; krok(n+1)=1.726346e-003; ng(n+1)=1.340854e+000;
n=79; farx(n+1)=2.776127e-001; foe(n+1)=4.440712e+000; krok(n+1)=4.000365e-002; ng(n+1)=6.427584e-001;
n=80; farx(n+1)=2.773039e-001; foe(n+1)=4.501207e+000; krok(n+1)=2.788178e-002; ng(n+1)=5.722896e-001;
n=81; farx(n+1)=2.760651e-001; foe(n+1)=4.416558e+000; krok(n+1)=3.452721e-002; ng(n+1)=1.458928e+000;
n=82; farx(n+1)=2.756259e-001; foe(n+1)=4.486226e+000; krok(n+1)=7.799474e-002; ng(n+1)=1.076706e+001;
n=83; farx(n+1)=2.741502e-001; foe(n+1)=4.373348e+000; krok(n+1)=2.320583e-001; ng(n+1)=1.580279e+001;
n=84; farx(n+1)=2.733245e-001; foe(n+1)=4.640503e+000; krok(n+1)=6.224674e-002; ng(n+1)=1.319155e+001;
n=85; farx(n+1)=2.728116e-001; foe(n+1)=4.588279e+000; krok(n+1)=2.114117e-001; ng(n+1)=9.250433e+000;
n=86; farx(n+1)=2.706620e-001; foe(n+1)=4.621172e+000; krok(n+1)=4.406235e-001; ng(n+1)=8.701722e+000;
n=87; farx(n+1)=2.687127e-001; foe(n+1)=4.477853e+000; krok(n+1)=8.029948e-001; ng(n+1)=3.200545e+000;
n=88; farx(n+1)=2.670004e-001; foe(n+1)=4.495184e+000; krok(n+1)=3.691454e-001; ng(n+1)=8.551262e+000;
n=89; farx(n+1)=2.660755e-001; foe(n+1)=4.324594e+000; krok(n+1)=2.114117e-001; ng(n+1)=4.696505e+000;
n=90; farx(n+1)=2.640069e-001; foe(n+1)=3.680910e+000; krok(n+1)=2.286351e-001; ng(n+1)=9.613455e+000;
n=91; farx(n+1)=2.610691e-001; foe(n+1)=3.519542e+000; krok(n+1)=5.956610e-001; ng(n+1)=9.993829e+000;
n=92; farx(n+1)=2.587278e-001; foe(n+1)=3.553811e+000; krok(n+1)=2.048545e-001; ng(n+1)=7.839376e+000;
n=93; farx(n+1)=2.565115e-001; foe(n+1)=3.540309e+000; krok(n+1)=4.686024e-001; ng(n+1)=6.186798e+000;
n=94; farx(n+1)=2.544932e-001; foe(n+1)=3.804211e+000; krok(n+1)=5.352197e-001; ng(n+1)=6.025591e+000;
n=95; farx(n+1)=2.528805e-001; foe(n+1)=3.840820e+000; krok(n+1)=1.653712e-001; ng(n+1)=1.201417e+001;
n=96; farx(n+1)=2.509942e-001; foe(n+1)=3.670708e+000; krok(n+1)=3.806077e-001; ng(n+1)=1.387347e+001;
n=97; farx(n+1)=2.485167e-001; foe(n+1)=3.734858e+000; krok(n+1)=9.934522e-001; ng(n+1)=6.611355e+000;
n=98; farx(n+1)=2.470211e-001; foe(n+1)=3.330943e+000; krok(n+1)=2.863419e-001; ng(n+1)=6.278982e+000;
n=99; farx(n+1)=2.462501e-001; foe(n+1)=3.234869e+000; krok(n+1)=2.418674e-001; ng(n+1)=6.743065e+000;
n=100; farx(n+1)=2.451061e-001; foe(n+1)=3.181849e+000; krok(n+1)=1.738735e-001; ng(n+1)=1.147346e+001;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora ARX');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
