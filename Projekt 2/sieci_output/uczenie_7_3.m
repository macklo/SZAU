%uczenie predyktora oe
clear all;
n=0; farx(n+1)=3.923223e+003; foe(n+1)=3.984579e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.629732e+003; foe(n+1)=2.739779e+003; krok(n+1)=5.188110e-004; ng(n+1)=4.350435e+003;
n=2; farx(n+1)=1.017573e+003; foe(n+1)=1.127627e+003; krok(n+1)=1.454759e-003; ng(n+1)=3.272105e+003;
n=3; farx(n+1)=1.283623e+003; foe(n+1)=6.587955e+002; krok(n+1)=1.650796e-004; ng(n+1)=9.477435e+003;
n=4; farx(n+1)=1.834092e+003; foe(n+1)=5.596694e+002; krok(n+1)=5.747417e-004; ng(n+1)=6.766179e+003;
n=5; farx(n+1)=1.593633e+003; foe(n+1)=5.035143e+002; krok(n+1)=5.557598e-003; ng(n+1)=9.129568e+002;
n=6; farx(n+1)=1.539719e+003; foe(n+1)=4.981444e+002; krok(n+1)=5.688063e-004; ng(n+1)=9.824566e+002;
n=7; farx(n+1)=7.025078e+002; foe(n+1)=4.018327e+002; krok(n+1)=1.291963e-003; ng(n+1)=1.398122e+003;
n=8; farx(n+1)=4.901011e+002; foe(n+1)=3.763648e+002; krok(n+1)=2.911585e-005; ng(n+1)=1.343095e+004;
n=9; farx(n+1)=4.707003e+002; foe(n+1)=3.615278e+002; krok(n+1)=2.044712e-004; ng(n+1)=1.506577e+004;
n=10; farx(n+1)=4.239506e+002; foe(n+1)=3.460822e+002; krok(n+1)=2.241713e-003; ng(n+1)=6.243382e+003;
n=11; farx(n+1)=2.960096e+002; foe(n+1)=3.134295e+002; krok(n+1)=1.046652e-002; ng(n+1)=1.370046e+003;
n=12; farx(n+1)=2.958280e+002; foe(n+1)=3.131833e+002; krok(n+1)=2.555890e-005; ng(n+1)=1.759263e+003;
n=13; farx(n+1)=2.748268e+002; foe(n+1)=3.058403e+002; krok(n+1)=2.647020e-003; ng(n+1)=1.563415e+003;
n=14; farx(n+1)=2.832979e+002; foe(n+1)=2.798038e+002; krok(n+1)=7.680234e-003; ng(n+1)=3.162138e+003;
n=15; farx(n+1)=2.963500e+002; foe(n+1)=2.768019e+002; krok(n+1)=1.357061e-004; ng(n+1)=2.738391e+003;
n=16; farx(n+1)=2.923303e+002; foe(n+1)=2.759273e+002; krok(n+1)=3.229906e-004; ng(n+1)=5.306995e+002;
n=17; farx(n+1)=2.822574e+002; foe(n+1)=2.678822e+002; krok(n+1)=1.830478e-003; ng(n+1)=1.602910e+003;
n=18; farx(n+1)=2.813052e+002; foe(n+1)=2.618341e+002; krok(n+1)=5.608389e-003; ng(n+1)=1.995909e+003;
n=19; farx(n+1)=2.503273e+002; foe(n+1)=2.455367e+002; krok(n+1)=5.355898e-003; ng(n+1)=1.493973e+003;
n=20; farx(n+1)=2.362536e+002; foe(n+1)=2.429066e+002; krok(n+1)=7.624904e-004; ng(n+1)=1.472997e+003;
n=21; farx(n+1)=2.171613e+002; foe(n+1)=2.371399e+002; krok(n+1)=1.081189e-003; ng(n+1)=5.094245e+003;
n=22; farx(n+1)=2.149637e+002; foe(n+1)=2.365485e+002; krok(n+1)=1.180809e-004; ng(n+1)=3.923561e+003;
n=23; farx(n+1)=1.675520e+002; foe(n+1)=2.132620e+002; krok(n+1)=1.713076e-002; ng(n+1)=3.630413e+003;
n=24; farx(n+1)=1.637599e+002; foe(n+1)=2.117447e+002; krok(n+1)=5.491616e-005; ng(n+1)=3.532113e+003;
n=25; farx(n+1)=1.528319e+002; foe(n+1)=1.978305e+002; krok(n+1)=1.489680e-003; ng(n+1)=3.775350e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=1.527870e+002; foe(n+1)=1.977767e+002; krok(n+1)=4.863692e-008; ng(n+1)=4.599613e+003;
n=27; farx(n+1)=1.186554e+002; foe(n+1)=1.490185e+002; krok(n+1)=4.934594e-005; ng(n+1)=4.162310e+003;
n=28; farx(n+1)=1.229831e+002; foe(n+1)=1.310188e+002; krok(n+1)=9.441524e-005; ng(n+1)=1.977489e+003;
n=29; farx(n+1)=9.824531e+001; foe(n+1)=1.112323e+002; krok(n+1)=2.471070e-004; ng(n+1)=1.855495e+003;
n=30; farx(n+1)=9.278352e+001; foe(n+1)=1.077358e+002; krok(n+1)=2.160041e-004; ng(n+1)=4.465556e+003;
n=31; farx(n+1)=4.117079e+001; foe(n+1)=6.837914e+001; krok(n+1)=5.836445e-003; ng(n+1)=4.474383e+003;
n=32; farx(n+1)=3.870479e+001; foe(n+1)=6.569311e+001; krok(n+1)=2.266194e-004; ng(n+1)=9.570971e+002;
n=33; farx(n+1)=3.061567e+001; foe(n+1)=5.913308e+001; krok(n+1)=1.420248e-003; ng(n+1)=7.981068e+002;
n=34; farx(n+1)=2.621130e+001; foe(n+1)=5.332330e+001; krok(n+1)=1.265711e-003; ng(n+1)=1.551057e+003;
n=35; farx(n+1)=2.109945e+001; foe(n+1)=4.425306e+001; krok(n+1)=2.707062e-003; ng(n+1)=8.362102e+002;
n=36; farx(n+1)=7.786492e+000; foe(n+1)=2.079515e+001; krok(n+1)=9.041785e-003; ng(n+1)=2.051243e+003;
n=37; farx(n+1)=5.264528e+000; foe(n+1)=1.263088e+001; krok(n+1)=8.017326e-004; ng(n+1)=1.663445e+003;
n=38; farx(n+1)=4.869591e+000; foe(n+1)=1.164247e+001; krok(n+1)=8.790130e-004; ng(n+1)=4.283520e+002;
n=39; farx(n+1)=4.493878e+000; foe(n+1)=1.020207e+001; krok(n+1)=3.135616e-003; ng(n+1)=4.818249e+002;
n=40; farx(n+1)=4.202667e+000; foe(n+1)=9.893938e+000; krok(n+1)=4.473296e-003; ng(n+1)=3.876897e+002;
n=41; farx(n+1)=3.946917e+000; foe(n+1)=8.814526e+000; krok(n+1)=2.399905e-002; ng(n+1)=4.320008e+002;
n=42; farx(n+1)=3.827893e+000; foe(n+1)=8.314250e+000; krok(n+1)=1.004449e-002; ng(n+1)=4.523853e+002;
n=43; farx(n+1)=3.672635e+000; foe(n+1)=7.947739e+000; krok(n+1)=2.476948e-002; ng(n+1)=1.980522e+002;
n=44; farx(n+1)=3.375019e+000; foe(n+1)=7.081064e+000; krok(n+1)=1.945154e-002; ng(n+1)=1.079147e+002;
n=45; farx(n+1)=3.230470e+000; foe(n+1)=6.735055e+000; krok(n+1)=1.547432e-002; ng(n+1)=4.568971e+002;
n=46; farx(n+1)=3.126035e+000; foe(n+1)=6.226101e+000; krok(n+1)=2.757245e-002; ng(n+1)=2.115857e+002;
n=47; farx(n+1)=3.164725e+000; foe(n+1)=5.957677e+000; krok(n+1)=3.661089e-002; ng(n+1)=1.797341e+002;
n=48; farx(n+1)=3.118589e+000; foe(n+1)=5.879199e+000; krok(n+1)=3.134999e-002; ng(n+1)=1.920844e+002;
n=49; farx(n+1)=2.920485e+000; foe(n+1)=5.527407e+000; krok(n+1)=1.372911e-001; ng(n+1)=2.630044e+002;
n=50; farx(n+1)=2.788021e+000; foe(n+1)=5.354914e+000; krok(n+1)=1.241631e-001; ng(n+1)=2.020965e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=2.784604e+000; foe(n+1)=5.337804e+000; krok(n+1)=6.169752e-006; ng(n+1)=1.980658e+002;
n=52; farx(n+1)=2.783451e+000; foe(n+1)=5.334283e+000; krok(n+1)=4.285449e-005; ng(n+1)=4.775611e+001;
n=53; farx(n+1)=2.792365e+000; foe(n+1)=5.322977e+000; krok(n+1)=2.130969e-004; ng(n+1)=4.045822e+001;
n=54; farx(n+1)=2.798997e+000; foe(n+1)=5.289327e+000; krok(n+1)=8.132823e-004; ng(n+1)=3.215323e+001;
n=55; farx(n+1)=2.823159e+000; foe(n+1)=5.234124e+000; krok(n+1)=3.739948e-003; ng(n+1)=1.867996e+001;
n=56; farx(n+1)=2.818212e+000; foe(n+1)=5.134804e+000; krok(n+1)=1.561071e-003; ng(n+1)=5.270438e+001;
n=57; farx(n+1)=2.799715e+000; foe(n+1)=5.074548e+000; krok(n+1)=1.023953e-002; ng(n+1)=2.045721e+002;
n=58; farx(n+1)=2.686729e+000; foe(n+1)=4.851066e+000; krok(n+1)=1.210728e-002; ng(n+1)=3.232801e+002;
n=59; farx(n+1)=2.550063e+000; foe(n+1)=4.622447e+000; krok(n+1)=1.321526e-002; ng(n+1)=2.133940e+002;
n=60; farx(n+1)=2.508408e+000; foe(n+1)=4.522300e+000; krok(n+1)=8.237127e-003; ng(n+1)=3.806357e+002;
n=61; farx(n+1)=2.477255e+000; foe(n+1)=4.419811e+000; krok(n+1)=2.223039e-002; ng(n+1)=1.961250e+002;
n=62; farx(n+1)=2.427790e+000; foe(n+1)=4.306721e+000; krok(n+1)=1.614645e-002; ng(n+1)=1.441838e+002;
n=63; farx(n+1)=2.459303e+000; foe(n+1)=4.134088e+000; krok(n+1)=2.080702e-002; ng(n+1)=2.263504e+002;
n=64; farx(n+1)=2.439672e+000; foe(n+1)=4.075666e+000; krok(n+1)=7.242402e-003; ng(n+1)=3.130541e+002;
n=65; farx(n+1)=2.278563e+000; foe(n+1)=3.804462e+000; krok(n+1)=4.320574e-002; ng(n+1)=3.566907e+002;
n=66; farx(n+1)=2.243514e+000; foe(n+1)=3.728460e+000; krok(n+1)=3.127439e-003; ng(n+1)=2.980162e+002;
n=67; farx(n+1)=2.198888e+000; foe(n+1)=3.643623e+000; krok(n+1)=2.065111e-002; ng(n+1)=1.044559e+002;
n=68; farx(n+1)=2.108658e+000; foe(n+1)=3.462785e+000; krok(n+1)=3.264474e-002; ng(n+1)=1.590723e+002;
n=69; farx(n+1)=2.058781e+000; foe(n+1)=3.280282e+000; krok(n+1)=4.130221e-002; ng(n+1)=1.447331e+002;
n=70; farx(n+1)=2.012468e+000; foe(n+1)=3.142983e+000; krok(n+1)=2.525939e-002; ng(n+1)=2.945116e+002;
n=71; farx(n+1)=1.976395e+000; foe(n+1)=3.073560e+000; krok(n+1)=3.468850e-002; ng(n+1)=8.952543e+001;
n=72; farx(n+1)=1.908281e+000; foe(n+1)=2.974067e+000; krok(n+1)=3.300845e-002; ng(n+1)=9.236718e+001;
n=73; farx(n+1)=1.736225e+000; foe(n+1)=2.782822e+000; krok(n+1)=1.047148e-001; ng(n+1)=7.888625e+001;
n=74; farx(n+1)=1.653050e+000; foe(n+1)=2.668387e+000; krok(n+1)=1.088038e-001; ng(n+1)=1.984653e+002;
n=75; farx(n+1)=1.621519e+000; foe(n+1)=2.547225e+000; krok(n+1)=1.778431e-001; ng(n+1)=1.358848e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.620325e+000; foe(n+1)=2.537705e+000; krok(n+1)=6.971762e-006; ng(n+1)=1.569464e+002;
n=77; farx(n+1)=1.620695e+000; foe(n+1)=2.536694e+000; krok(n+1)=1.259177e-005; ng(n+1)=3.980462e+001;
n=78; farx(n+1)=1.621677e+000; foe(n+1)=2.531579e+000; krok(n+1)=1.973390e-004; ng(n+1)=2.719469e+001;
n=79; farx(n+1)=1.612493e+000; foe(n+1)=2.517215e+000; krok(n+1)=8.605927e-004; ng(n+1)=2.020979e+001;
n=80; farx(n+1)=1.605520e+000; foe(n+1)=2.479798e+000; krok(n+1)=1.081070e-003; ng(n+1)=2.713883e+001;
n=81; farx(n+1)=1.593965e+000; foe(n+1)=2.462703e+000; krok(n+1)=1.854565e-003; ng(n+1)=2.132201e+001;
n=82; farx(n+1)=1.536660e+000; foe(n+1)=2.379742e+000; krok(n+1)=1.957537e-002; ng(n+1)=1.993920e+001;
n=83; farx(n+1)=1.526457e+000; foe(n+1)=2.362354e+000; krok(n+1)=2.716773e-003; ng(n+1)=9.585302e+001;
n=84; farx(n+1)=1.483330e+000; foe(n+1)=2.328435e+000; krok(n+1)=1.464382e-002; ng(n+1)=1.176167e+002;
n=85; farx(n+1)=1.459413e+000; foe(n+1)=2.290785e+000; krok(n+1)=2.000182e-002; ng(n+1)=1.166240e+002;
n=86; farx(n+1)=1.444332e+000; foe(n+1)=2.268523e+000; krok(n+1)=1.064021e-002; ng(n+1)=1.789135e+002;
n=87; farx(n+1)=1.414231e+000; foe(n+1)=2.215312e+000; krok(n+1)=3.012318e-002; ng(n+1)=8.866144e+001;
n=88; farx(n+1)=1.384568e+000; foe(n+1)=2.160980e+000; krok(n+1)=1.064021e-002; ng(n+1)=2.712335e+002;
n=89; farx(n+1)=1.362424e+000; foe(n+1)=2.100428e+000; krok(n+1)=1.974372e-002; ng(n+1)=2.593735e+002;
n=90; farx(n+1)=1.343847e+000; foe(n+1)=2.080292e+000; krok(n+1)=1.889482e-002; ng(n+1)=1.114407e+002;
n=91; farx(n+1)=1.331568e+000; foe(n+1)=2.057185e+000; krok(n+1)=5.433546e-003; ng(n+1)=1.955353e+002;
n=92; farx(n+1)=1.309173e+000; foe(n+1)=2.017911e+000; krok(n+1)=6.222343e-002; ng(n+1)=8.890522e+001;
n=93; farx(n+1)=1.264996e+000; foe(n+1)=1.965347e+000; krok(n+1)=3.923333e-002; ng(n+1)=1.551296e+002;
n=94; farx(n+1)=1.265132e+000; foe(n+1)=1.943700e+000; krok(n+1)=3.127687e-002; ng(n+1)=3.127679e+002;
n=95; farx(n+1)=1.285265e+000; foe(n+1)=1.879648e+000; krok(n+1)=1.016482e-001; ng(n+1)=2.281566e+002;
n=96; farx(n+1)=1.264949e+000; foe(n+1)=1.833473e+000; krok(n+1)=8.589342e-002; ng(n+1)=1.121545e+002;
n=97; farx(n+1)=1.272490e+000; foe(n+1)=1.799357e+000; krok(n+1)=1.242132e-001; ng(n+1)=1.294350e+002;
n=98; farx(n+1)=1.278333e+000; foe(n+1)=1.769138e+000; krok(n+1)=4.331299e-002; ng(n+1)=1.161353e+002;
n=99; farx(n+1)=1.282559e+000; foe(n+1)=1.724170e+000; krok(n+1)=5.411858e-002; ng(n+1)=3.989303e+002;
n=100; farx(n+1)=1.256490e+000; foe(n+1)=1.678361e+000; krok(n+1)=1.027628e-001; ng(n+1)=1.604198e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)