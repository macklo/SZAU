%uczenie predyktora oe
clear all;
n=0; farx(n+1)=5.006360e+003; foe(n+1)=4.886956e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.922632e+003; foe(n+1)=3.876772e+003; krok(n+1)=3.473499e-004; ng(n+1)=2.597530e+003;
n=2; farx(n+1)=8.103088e+002; foe(n+1)=7.963431e+002; krok(n+1)=3.831275e-003; ng(n+1)=1.109854e+003;
n=3; farx(n+1)=8.060322e+002; foe(n+1)=7.783936e+002; krok(n+1)=1.813010e-004; ng(n+1)=8.720779e+002;
n=4; farx(n+1)=7.428056e+002; foe(n+1)=7.738482e+002; krok(n+1)=2.954608e-004; ng(n+1)=7.321185e+002;
n=5; farx(n+1)=4.054082e+002; foe(n+1)=5.410540e+002; krok(n+1)=2.320282e-003; ng(n+1)=1.349843e+003;
n=6; farx(n+1)=3.006169e+002; foe(n+1)=4.759773e+002; krok(n+1)=3.952630e-004; ng(n+1)=2.843720e+003;
n=7; farx(n+1)=2.286766e+002; foe(n+1)=4.339089e+002; krok(n+1)=5.604282e-004; ng(n+1)=3.796636e+003;
n=8; farx(n+1)=2.203325e+002; foe(n+1)=4.299944e+002; krok(n+1)=9.532169e-005; ng(n+1)=6.241471e+003;
n=9; farx(n+1)=2.006890e+002; foe(n+1)=3.794515e+002; krok(n+1)=8.442897e-003; ng(n+1)=6.301503e+003;
n=10; farx(n+1)=2.008507e+002; foe(n+1)=3.782934e+002; krok(n+1)=1.221743e-005; ng(n+1)=3.734412e+003;
n=11; farx(n+1)=1.438642e+002; foe(n+1)=3.396437e+002; krok(n+1)=1.048565e-002; ng(n+1)=4.340166e+003;
n=12; farx(n+1)=1.085478e+002; foe(n+1)=2.743056e+002; krok(n+1)=6.643424e-003; ng(n+1)=5.459452e+003;
n=13; farx(n+1)=1.068127e+002; foe(n+1)=2.737372e+002; krok(n+1)=7.910691e-005; ng(n+1)=3.111397e+003;
n=14; farx(n+1)=1.107097e+002; foe(n+1)=2.729564e+002; krok(n+1)=8.058731e-004; ng(n+1)=2.604118e+003;
n=15; farx(n+1)=9.842372e+001; foe(n+1)=2.632571e+002; krok(n+1)=4.614317e-002; ng(n+1)=2.161191e+003;
n=16; farx(n+1)=8.966268e+001; foe(n+1)=2.532150e+002; krok(n+1)=2.556212e-002; ng(n+1)=2.417076e+003;
n=17; farx(n+1)=6.163754e+001; foe(n+1)=2.212241e+002; krok(n+1)=2.069473e-002; ng(n+1)=3.632404e+003;
n=18; farx(n+1)=5.678385e+001; foe(n+1)=2.070078e+002; krok(n+1)=1.045039e-003; ng(n+1)=1.961978e+004;
n=19; farx(n+1)=4.243872e+001; foe(n+1)=1.869111e+002; krok(n+1)=4.438543e-002; ng(n+1)=2.894824e+004;
n=20; farx(n+1)=3.721955e+001; foe(n+1)=1.625888e+002; krok(n+1)=2.798246e-001; ng(n+1)=2.453774e+004;
n=21; farx(n+1)=3.473866e+001; foe(n+1)=1.509781e+002; krok(n+1)=6.813219e-002; ng(n+1)=1.343517e+004;
n=22; farx(n+1)=3.276265e+001; foe(n+1)=1.445724e+002; krok(n+1)=1.755320e-001; ng(n+1)=3.320332e+003;
n=23; farx(n+1)=3.615816e+001; foe(n+1)=1.197110e+002; krok(n+1)=8.191626e-002; ng(n+1)=1.031320e+004;
n=24; farx(n+1)=3.398372e+001; foe(n+1)=1.127883e+002; krok(n+1)=2.056099e-001; ng(n+1)=3.827712e+003;
n=25; farx(n+1)=3.229589e+001; foe(n+1)=1.053607e+002; krok(n+1)=3.007188e-001; ng(n+1)=2.009831e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=3.196088e+001; foe(n+1)=1.048499e+002; krok(n+1)=8.187860e-007; ng(n+1)=2.874014e+003;
n=27; farx(n+1)=3.089973e+001; foe(n+1)=1.038949e+002; krok(n+1)=4.765565e-005; ng(n+1)=4.888840e+002;
n=28; farx(n+1)=2.719876e+001; foe(n+1)=1.004123e+002; krok(n+1)=2.247378e-004; ng(n+1)=4.515591e+002;
n=29; farx(n+1)=2.455371e+001; foe(n+1)=9.886595e+001; krok(n+1)=1.614953e-004; ng(n+1)=5.408880e+002;
n=30; farx(n+1)=2.290488e+001; foe(n+1)=9.109619e+001; krok(n+1)=1.342085e-003; ng(n+1)=1.363924e+003;
n=31; farx(n+1)=1.859979e+001; foe(n+1)=7.894997e+001; krok(n+1)=1.225263e-003; ng(n+1)=6.450937e+003;
n=32; farx(n+1)=1.482099e+001; foe(n+1)=7.507352e+001; krok(n+1)=9.915259e-004; ng(n+1)=1.525378e+003;
n=33; farx(n+1)=1.312078e+001; foe(n+1)=7.052005e+001; krok(n+1)=1.369522e-003; ng(n+1)=2.293177e+003;
n=34; farx(n+1)=1.401529e+001; foe(n+1)=6.793673e+001; krok(n+1)=2.804195e-003; ng(n+1)=3.119278e+003;
n=35; farx(n+1)=1.460236e+001; foe(n+1)=6.680967e+001; krok(n+1)=8.786586e-004; ng(n+1)=1.896883e+003;
n=36; farx(n+1)=1.663510e+001; foe(n+1)=6.496358e+001; krok(n+1)=1.033570e-002; ng(n+1)=4.930320e+002;
n=37; farx(n+1)=1.766489e+001; foe(n+1)=6.396346e+001; krok(n+1)=4.051549e-003; ng(n+1)=1.117527e+003;
n=38; farx(n+1)=1.843459e+001; foe(n+1)=6.312665e+001; krok(n+1)=4.736926e-003; ng(n+1)=2.514249e+003;
n=39; farx(n+1)=1.506521e+001; foe(n+1)=5.886082e+001; krok(n+1)=7.000904e-002; ng(n+1)=1.038431e+003;
n=40; farx(n+1)=1.456421e+001; foe(n+1)=5.844925e+001; krok(n+1)=1.901096e-002; ng(n+1)=9.287245e+002;
n=41; farx(n+1)=1.347948e+001; foe(n+1)=5.665958e+001; krok(n+1)=1.253365e-001; ng(n+1)=9.408954e+002;
n=42; farx(n+1)=1.297448e+001; foe(n+1)=5.413067e+001; krok(n+1)=2.155114e-001; ng(n+1)=1.303851e+003;
n=43; farx(n+1)=1.187093e+001; foe(n+1)=5.196142e+001; krok(n+1)=3.044248e-001; ng(n+1)=2.088134e+003;
n=44; farx(n+1)=1.260965e+001; foe(n+1)=5.088530e+001; krok(n+1)=5.114232e-001; ng(n+1)=6.977424e+002;
n=45; farx(n+1)=1.286386e+001; foe(n+1)=5.053577e+001; krok(n+1)=6.109364e-001; ng(n+1)=9.679835e+002;
n=46; farx(n+1)=1.256236e+001; foe(n+1)=5.013987e+001; krok(n+1)=4.128102e-001; ng(n+1)=2.032774e+002;
n=47; farx(n+1)=1.143320e+001; foe(n+1)=4.916486e+001; krok(n+1)=1.150553e+000; ng(n+1)=6.902770e+002;
n=48; farx(n+1)=9.061998e+000; foe(n+1)=4.757148e+001; krok(n+1)=9.495372e-001; ng(n+1)=9.694060e+002;
n=49; farx(n+1)=7.331076e+000; foe(n+1)=4.648556e+001; krok(n+1)=2.192961e-001; ng(n+1)=9.749053e+002;
n=50; farx(n+1)=5.448301e+000; foe(n+1)=4.439356e+001; krok(n+1)=2.896177e-001; ng(n+1)=9.678033e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=5.415957e+000; foe(n+1)=4.437524e+001; krok(n+1)=1.099167e-005; ng(n+1)=1.352012e+002;
n=52; farx(n+1)=5.308744e+000; foe(n+1)=4.428093e+001; krok(n+1)=1.197415e-005; ng(n+1)=3.350305e+002;
n=53; farx(n+1)=5.185818e+000; foe(n+1)=4.417997e+001; krok(n+1)=1.304908e-005; ng(n+1)=3.372821e+002;
n=54; farx(n+1)=5.128736e+000; foe(n+1)=4.317695e+001; krok(n+1)=2.084136e-004; ng(n+1)=2.556449e+002;
n=55; farx(n+1)=5.137387e+000; foe(n+1)=4.230112e+001; krok(n+1)=1.524981e-003; ng(n+1)=1.258973e+002;
n=56; farx(n+1)=4.960922e+000; foe(n+1)=4.174262e+001; krok(n+1)=1.043353e-003; ng(n+1)=1.501101e+002;
n=57; farx(n+1)=4.895867e+000; foe(n+1)=4.141765e+001; krok(n+1)=3.174920e-003; ng(n+1)=3.419842e+002;
n=58; farx(n+1)=4.574161e+000; foe(n+1)=4.059219e+001; krok(n+1)=1.617035e-002; ng(n+1)=1.123539e+003;
n=59; farx(n+1)=4.466081e+000; foe(n+1)=4.035289e+001; krok(n+1)=1.504792e-003; ng(n+1)=9.107452e+002;
n=60; farx(n+1)=4.448547e+000; foe(n+1)=4.001457e+001; krok(n+1)=1.321323e-002; ng(n+1)=1.000926e+003;
n=61; farx(n+1)=4.296741e+000; foe(n+1)=3.981271e+001; krok(n+1)=1.923658e-002; ng(n+1)=6.028271e+002;
n=62; farx(n+1)=4.174911e+000; foe(n+1)=3.973619e+001; krok(n+1)=1.617035e-002; ng(n+1)=4.809879e+002;
n=63; farx(n+1)=4.143766e+000; foe(n+1)=3.968000e+001; krok(n+1)=6.669234e-003; ng(n+1)=2.468264e+002;
n=64; farx(n+1)=3.594897e+000; foe(n+1)=3.934368e+001; krok(n+1)=1.138873e-001; ng(n+1)=7.008249e+002;
n=65; farx(n+1)=3.671726e+000; foe(n+1)=3.875546e+001; krok(n+1)=3.222183e-002; ng(n+1)=4.861245e+002;
n=66; farx(n+1)=3.935593e+000; foe(n+1)=3.811281e+001; krok(n+1)=3.897625e-002; ng(n+1)=8.516721e+002;
n=67; farx(n+1)=3.889994e+000; foe(n+1)=3.726201e+001; krok(n+1)=3.427774e-001; ng(n+1)=5.367638e+002;
n=68; farx(n+1)=3.913899e+000; foe(n+1)=3.674129e+001; krok(n+1)=1.196783e-001; ng(n+1)=1.123153e+003;
n=69; farx(n+1)=4.598638e+000; foe(n+1)=3.532437e+001; krok(n+1)=3.161324e-001; ng(n+1)=3.353247e+002;
n=70; farx(n+1)=5.284294e+000; foe(n+1)=3.460215e+001; krok(n+1)=4.236083e-002; ng(n+1)=1.336179e+003;
n=71; farx(n+1)=5.976930e+000; foe(n+1)=3.148941e+001; krok(n+1)=4.027995e-001; ng(n+1)=2.128590e+003;
n=72; farx(n+1)=5.860552e+000; foe(n+1)=2.994345e+001; krok(n+1)=3.879553e-002; ng(n+1)=2.285734e+003;
n=73; farx(n+1)=5.325549e+000; foe(n+1)=2.581615e+001; krok(n+1)=1.134226e-001; ng(n+1)=2.645824e+003;
n=74; farx(n+1)=4.874733e+000; foe(n+1)=2.088467e+001; krok(n+1)=1.016482e-001; ng(n+1)=4.060121e+003;
n=75; farx(n+1)=5.150642e+000; foe(n+1)=1.650003e+001; krok(n+1)=1.607119e-001; ng(n+1)=1.252456e+003;
%odnowa zmiennej metryki
n=76; farx(n+1)=5.133195e+000; foe(n+1)=1.553741e+001; krok(n+1)=7.466456e-007; ng(n+1)=3.902030e+003;
n=77; farx(n+1)=5.121543e+000; foe(n+1)=1.540181e+001; krok(n+1)=5.846756e-006; ng(n+1)=5.514118e+002;
n=78; farx(n+1)=4.897094e+000; foe(n+1)=1.359526e+001; krok(n+1)=3.754715e-005; ng(n+1)=7.593051e+002;
n=79; farx(n+1)=5.076604e+000; foe(n+1)=1.233188e+001; krok(n+1)=3.923630e-005; ng(n+1)=6.537269e+002;
n=80; farx(n+1)=5.072542e+000; foe(n+1)=9.851372e+000; krok(n+1)=6.847610e-004; ng(n+1)=2.575720e+002;
n=81; farx(n+1)=5.052454e+000; foe(n+1)=9.043109e+000; krok(n+1)=1.298034e-004; ng(n+1)=4.979871e+002;
n=82; farx(n+1)=4.617965e+000; foe(n+1)=8.568814e+000; krok(n+1)=2.569527e-003; ng(n+1)=1.052145e+003;
n=83; farx(n+1)=4.284806e+000; foe(n+1)=8.323014e+000; krok(n+1)=7.242402e-003; ng(n+1)=7.464286e+002;
n=84; farx(n+1)=4.272688e+000; foe(n+1)=7.852227e+000; krok(n+1)=7.739087e-003; ng(n+1)=5.070319e+002;
n=85; farx(n+1)=4.132343e+000; foe(n+1)=7.622860e+000; krok(n+1)=5.062842e-003; ng(n+1)=9.925293e+002;
n=86; farx(n+1)=3.953737e+000; foe(n+1)=7.579936e+000; krok(n+1)=1.023953e-002; ng(n+1)=4.228528e+002;
n=87; farx(n+1)=3.912000e+000; foe(n+1)=7.474676e+000; krok(n+1)=1.948770e-002; ng(n+1)=7.192574e+002;
n=88; farx(n+1)=4.237488e+000; foe(n+1)=7.342511e+000; krok(n+1)=3.497807e-002; ng(n+1)=8.923157e+002;
n=89; farx(n+1)=4.333805e+000; foe(n+1)=7.147231e+000; krok(n+1)=9.175071e-003; ng(n+1)=9.339201e+002;
n=90; farx(n+1)=4.007956e+000; foe(n+1)=7.052627e+000; krok(n+1)=3.450322e-002; ng(n+1)=1.626083e+002;
n=91; farx(n+1)=3.354533e+000; foe(n+1)=6.632317e+000; krok(n+1)=1.291716e-001; ng(n+1)=1.398245e+002;
n=92; farx(n+1)=3.055647e+000; foe(n+1)=6.528251e+000; krok(n+1)=1.316760e-001; ng(n+1)=4.083630e+002;
n=93; farx(n+1)=2.907138e+000; foe(n+1)=6.262158e+000; krok(n+1)=2.114442e-001; ng(n+1)=4.964361e+002;
n=94; farx(n+1)=2.739554e+000; foe(n+1)=6.032608e+000; krok(n+1)=3.388866e-001; ng(n+1)=5.933051e+002;
n=95; farx(n+1)=2.288645e+000; foe(n+1)=5.691914e+000; krok(n+1)=2.357265e-001; ng(n+1)=1.072808e+003;
n=96; farx(n+1)=2.074620e+000; foe(n+1)=5.439479e+000; krok(n+1)=2.442288e-001; ng(n+1)=6.814470e+002;
n=97; farx(n+1)=1.930424e+000; foe(n+1)=4.987581e+000; krok(n+1)=4.364225e-001; ng(n+1)=3.561843e+002;
n=98; farx(n+1)=1.954558e+000; foe(n+1)=4.727271e+000; krok(n+1)=1.158784e-001; ng(n+1)=1.042404e+003;
n=99; farx(n+1)=2.114939e+000; foe(n+1)=4.500157e+000; krok(n+1)=3.041754e-001; ng(n+1)=7.920382e+002;
n=100; farx(n+1)=2.050528e+000; foe(n+1)=4.126325e+000; krok(n+1)=2.628275e-001; ng(n+1)=6.139026e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
