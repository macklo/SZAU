%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.897050e+003; foe(n+1)=4.828634e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.928939e+003; foe(n+1)=3.822986e+003; krok(n+1)=3.505243e-004; ng(n+1)=3.735712e+003;
n=2; farx(n+1)=1.079798e+003; foe(n+1)=8.721602e+002; krok(n+1)=2.641497e-003; ng(n+1)=2.294952e+003;
n=3; farx(n+1)=1.171740e+003; foe(n+1)=8.167219e+002; krok(n+1)=2.463545e-004; ng(n+1)=2.323875e+003;
n=4; farx(n+1)=4.353980e+002; foe(n+1)=6.766639e+002; krok(n+1)=7.557177e-003; ng(n+1)=7.241290e+002;
n=5; farx(n+1)=3.759910e+002; foe(n+1)=6.515699e+002; krok(n+1)=5.371265e-005; ng(n+1)=6.925769e+003;
n=6; farx(n+1)=2.998680e+002; foe(n+1)=4.816623e+002; krok(n+1)=9.064776e-004; ng(n+1)=6.507044e+003;
n=7; farx(n+1)=3.023797e+002; foe(n+1)=4.354862e+002; krok(n+1)=5.823131e-004; ng(n+1)=5.024900e+003;
n=8; farx(n+1)=2.119209e+002; foe(n+1)=3.705426e+002; krok(n+1)=2.999882e-003; ng(n+1)=4.653075e+003;
n=9; farx(n+1)=2.005071e+002; foe(n+1)=3.510765e+002; krok(n+1)=5.277180e-004; ng(n+1)=1.765896e+003;
n=10; farx(n+1)=2.136837e+002; foe(n+1)=3.247133e+002; krok(n+1)=7.032104e-003; ng(n+1)=1.104008e+003;
n=11; farx(n+1)=1.876172e+002; foe(n+1)=3.106222e+002; krok(n+1)=4.038148e-003; ng(n+1)=6.682385e+002;
n=12; farx(n+1)=1.836858e+002; foe(n+1)=3.078236e+002; krok(n+1)=7.988163e-004; ng(n+1)=8.396610e+002;
n=13; farx(n+1)=1.750572e+002; foe(n+1)=3.026686e+002; krok(n+1)=2.501815e-003; ng(n+1)=9.813011e+002;
n=14; farx(n+1)=1.527296e+002; foe(n+1)=2.886070e+002; krok(n+1)=5.557598e-003; ng(n+1)=6.958989e+002;
n=15; farx(n+1)=1.386202e+002; foe(n+1)=2.843908e+002; krok(n+1)=1.554197e-003; ng(n+1)=8.414886e+002;
n=16; farx(n+1)=9.384304e+001; foe(n+1)=2.700867e+002; krok(n+1)=6.681129e-003; ng(n+1)=7.040025e+002;
n=17; farx(n+1)=6.201271e+001; foe(n+1)=2.567970e+002; krok(n+1)=3.127269e-004; ng(n+1)=2.377167e+003;
n=18; farx(n+1)=5.338083e+001; foe(n+1)=2.500743e+002; krok(n+1)=7.543972e-004; ng(n+1)=3.886970e+003;
n=19; farx(n+1)=4.730835e+001; foe(n+1)=2.456100e+002; krok(n+1)=1.244450e-003; ng(n+1)=3.230554e+003;
n=20; farx(n+1)=4.466089e+001; foe(n+1)=2.434316e+002; krok(n+1)=6.935224e-004; ng(n+1)=2.389299e+003;
n=21; farx(n+1)=4.091586e+001; foe(n+1)=2.385546e+002; krok(n+1)=5.139053e-003; ng(n+1)=2.716621e+003;
n=22; farx(n+1)=3.902975e+001; foe(n+1)=2.350397e+002; krok(n+1)=2.332157e-003; ng(n+1)=2.120272e+003;
n=23; farx(n+1)=3.967161e+001; foe(n+1)=2.246361e+002; krok(n+1)=7.321912e-003; ng(n+1)=2.528792e+003;
n=24; farx(n+1)=3.750631e+001; foe(n+1)=2.163253e+002; krok(n+1)=2.569527e-003; ng(n+1)=1.767990e+003;
n=25; farx(n+1)=3.954313e+001; foe(n+1)=2.073650e+002; krok(n+1)=4.483426e-003; ng(n+1)=3.328657e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=3.558541e+001; foe(n+1)=1.962497e+002; krok(n+1)=1.381179e-004; ng(n+1)=1.356173e+003;
n=27; farx(n+1)=3.370561e+001; foe(n+1)=1.900932e+002; krok(n+1)=5.482924e-006; ng(n+1)=3.820807e+003;
n=28; farx(n+1)=3.162582e+001; foe(n+1)=1.834967e+002; krok(n+1)=1.856539e-005; ng(n+1)=1.788647e+003;
n=29; farx(n+1)=3.248610e+001; foe(n+1)=1.707180e+002; krok(n+1)=9.921626e-005; ng(n+1)=2.083438e+003;
n=30; farx(n+1)=3.406288e+001; foe(n+1)=1.501531e+002; krok(n+1)=9.700243e-004; ng(n+1)=2.378207e+003;
n=31; farx(n+1)=3.269145e+001; foe(n+1)=1.406804e+002; krok(n+1)=6.369581e-004; ng(n+1)=1.013967e+003;
n=32; farx(n+1)=2.562492e+001; foe(n+1)=1.293183e+002; krok(n+1)=7.821151e-004; ng(n+1)=8.409089e+002;
n=33; farx(n+1)=2.204715e+001; foe(n+1)=1.184324e+002; krok(n+1)=8.018005e-004; ng(n+1)=2.539939e+003;
n=34; farx(n+1)=1.769151e+001; foe(n+1)=8.584828e+001; krok(n+1)=1.861026e-003; ng(n+1)=1.784727e+003;
n=35; farx(n+1)=1.829817e+001; foe(n+1)=7.262157e+001; krok(n+1)=1.758026e-003; ng(n+1)=3.786114e+003;
n=36; farx(n+1)=1.454929e+001; foe(n+1)=5.339423e+001; krok(n+1)=3.334617e-003; ng(n+1)=2.305925e+003;
n=37; farx(n+1)=1.409373e+001; foe(n+1)=4.487107e+001; krok(n+1)=7.821151e-004; ng(n+1)=2.403319e+003;
n=38; farx(n+1)=1.021203e+001; foe(n+1)=3.740125e+001; krok(n+1)=4.747676e-003; ng(n+1)=4.801808e+002;
n=39; farx(n+1)=9.279286e+000; foe(n+1)=3.359032e+001; krok(n+1)=1.635770e-003; ng(n+1)=1.514444e+003;
n=40; farx(n+1)=6.458174e+000; foe(n+1)=2.569636e+001; krok(n+1)=8.534063e-003; ng(n+1)=1.306660e+003;
n=41; farx(n+1)=5.587755e+000; foe(n+1)=2.120043e+001; krok(n+1)=6.555717e-003; ng(n+1)=6.041306e+002;
n=42; farx(n+1)=5.202904e+000; foe(n+1)=1.947236e+001; krok(n+1)=9.157950e-003; ng(n+1)=7.768623e+002;
n=43; farx(n+1)=4.552775e+000; foe(n+1)=1.645002e+001; krok(n+1)=6.734730e-003; ng(n+1)=7.582818e+002;
n=44; farx(n+1)=4.824524e+000; foe(n+1)=1.530122e+001; krok(n+1)=5.119766e-003; ng(n+1)=1.074726e+003;
n=45; farx(n+1)=4.737737e+000; foe(n+1)=1.315496e+001; krok(n+1)=1.034736e-002; ng(n+1)=1.122482e+003;
n=46; farx(n+1)=4.369738e+000; foe(n+1)=1.178977e+001; krok(n+1)=4.672280e-003; ng(n+1)=1.278251e+003;
n=47; farx(n+1)=4.470106e+000; foe(n+1)=1.103873e+001; krok(n+1)=1.064021e-002; ng(n+1)=1.505762e+003;
n=48; farx(n+1)=4.360663e+000; foe(n+1)=1.040542e+001; krok(n+1)=1.881038e-002; ng(n+1)=8.143100e+002;
n=49; farx(n+1)=4.280992e+000; foe(n+1)=9.612223e+000; krok(n+1)=2.727640e-002; ng(n+1)=9.149028e+002;
n=50; farx(n+1)=3.471313e+000; foe(n+1)=8.242572e+000; krok(n+1)=5.266571e-002; ng(n+1)=5.696634e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=3.469390e+000; foe(n+1)=8.119054e+000; krok(n+1)=2.462151e-005; ng(n+1)=2.654664e+002;
n=52; farx(n+1)=3.461632e+000; foe(n+1)=8.065021e+000; krok(n+1)=1.649119e-005; ng(n+1)=2.324106e+002;
n=53; farx(n+1)=3.432233e+000; foe(n+1)=7.817806e+000; krok(n+1)=1.451493e-004; ng(n+1)=2.065080e+002;
n=54; farx(n+1)=3.329900e+000; foe(n+1)=7.582965e+000; krok(n+1)=2.594055e-004; ng(n+1)=1.410217e+002;
n=55; farx(n+1)=3.270143e+000; foe(n+1)=7.313112e+000; krok(n+1)=8.897443e-004; ng(n+1)=7.764202e+001;
n=56; farx(n+1)=3.065019e+000; foe(n+1)=6.905372e+000; krok(n+1)=2.647552e-003; ng(n+1)=7.850517e+001;
n=57; farx(n+1)=3.022451e+000; foe(n+1)=6.609995e+000; krok(n+1)=3.625910e-003; ng(n+1)=6.918477e+001;
n=58; farx(n+1)=2.809392e+000; foe(n+1)=6.228120e+000; krok(n+1)=1.017451e-002; ng(n+1)=2.234607e+002;
n=59; farx(n+1)=2.654364e+000; foe(n+1)=5.818075e+000; krok(n+1)=1.021537e-002; ng(n+1)=2.543922e+002;
n=60; farx(n+1)=2.440664e+000; foe(n+1)=5.328842e+000; krok(n+1)=6.353015e-003; ng(n+1)=5.317276e+002;
n=61; farx(n+1)=2.376873e+000; foe(n+1)=5.015085e+000; krok(n+1)=5.000456e-003; ng(n+1)=4.522894e+002;
n=62; farx(n+1)=2.224559e+000; foe(n+1)=4.613008e+000; krok(n+1)=1.650507e-002; ng(n+1)=5.518600e+002;
n=63; farx(n+1)=2.094576e+000; foe(n+1)=4.430371e+000; krok(n+1)=7.632150e-003; ng(n+1)=8.528358e+002;
n=64; farx(n+1)=1.941039e+000; foe(n+1)=4.248939e+000; krok(n+1)=2.343419e-002; ng(n+1)=8.736617e+002;
n=65; farx(n+1)=1.757075e+000; foe(n+1)=4.011965e+000; krok(n+1)=1.226253e-002; ng(n+1)=6.402926e+002;
n=66; farx(n+1)=1.655423e+000; foe(n+1)=3.890510e+000; krok(n+1)=7.626318e-003; ng(n+1)=5.220945e+002;
n=67; farx(n+1)=1.595817e+000; foe(n+1)=3.817087e+000; krok(n+1)=9.041785e-003; ng(n+1)=3.861298e+002;
n=68; farx(n+1)=1.515561e+000; foe(n+1)=3.728067e+000; krok(n+1)=1.289397e-002; ng(n+1)=1.851656e+002;
n=69; farx(n+1)=1.420498e+000; foe(n+1)=3.587994e+000; krok(n+1)=2.357847e-002; ng(n+1)=2.091166e+002;
n=70; farx(n+1)=1.316246e+000; foe(n+1)=3.434230e+000; krok(n+1)=1.596372e-002; ng(n+1)=1.641083e+002;
n=71; farx(n+1)=1.244839e+000; foe(n+1)=3.294222e+000; krok(n+1)=5.709755e-002; ng(n+1)=2.687346e+002;
n=72; farx(n+1)=1.176469e+000; foe(n+1)=3.005248e+000; krok(n+1)=1.269153e-001; ng(n+1)=1.697603e+002;
n=73; farx(n+1)=1.158954e+000; foe(n+1)=2.832618e+000; krok(n+1)=4.524124e-002; ng(n+1)=2.852444e+002;
n=74; farx(n+1)=1.204312e+000; foe(n+1)=2.479361e+000; krok(n+1)=1.032026e-001; ng(n+1)=3.583464e+002;
n=75; farx(n+1)=1.115432e+000; foe(n+1)=2.007417e+000; krok(n+1)=1.209148e-001; ng(n+1)=3.687229e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.112312e+000; foe(n+1)=1.981107e+000; krok(n+1)=1.191391e-005; ng(n+1)=1.836703e+002;
n=77; farx(n+1)=1.112772e+000; foe(n+1)=1.955887e+000; krok(n+1)=8.786087e-006; ng(n+1)=2.386567e+002;
n=78; farx(n+1)=1.116217e+000; foe(n+1)=1.940882e+000; krok(n+1)=1.416371e-005; ng(n+1)=1.450013e+002;
n=79; farx(n+1)=1.107063e+000; foe(n+1)=1.889042e+000; krok(n+1)=7.741729e-004; ng(n+1)=4.021918e+001;
n=80; farx(n+1)=1.068008e+000; foe(n+1)=1.779836e+000; krok(n+1)=6.647718e-004; ng(n+1)=5.700803e+001;
n=81; farx(n+1)=1.048496e+000; foe(n+1)=1.728070e+000; krok(n+1)=2.372391e-004; ng(n+1)=8.968618e+001;
n=82; farx(n+1)=1.029221e+000; foe(n+1)=1.699550e+000; krok(n+1)=3.242576e-003; ng(n+1)=1.174154e+002;
n=83; farx(n+1)=9.837985e-001; foe(n+1)=1.646192e+000; krok(n+1)=2.659087e-003; ng(n+1)=1.686033e+002;
n=84; farx(n+1)=9.622343e-001; foe(n+1)=1.617821e+000; krok(n+1)=5.849663e-003; ng(n+1)=2.269076e+002;
n=85; farx(n+1)=9.255558e-001; foe(n+1)=1.539262e+000; krok(n+1)=1.360047e-002; ng(n+1)=2.298478e+002;
n=86; farx(n+1)=8.977275e-001; foe(n+1)=1.453404e+000; krok(n+1)=8.108299e-003; ng(n+1)=8.459421e+001;
n=87; farx(n+1)=8.568299e-001; foe(n+1)=1.369658e+000; krok(n+1)=2.025137e-002; ng(n+1)=2.880541e+002;
n=88; farx(n+1)=8.436038e-001; foe(n+1)=1.327113e+000; krok(n+1)=1.032555e-002; ng(n+1)=2.453865e+002;
n=89; farx(n+1)=8.217872e-001; foe(n+1)=1.287751e+000; krok(n+1)=1.830045e-002; ng(n+1)=4.860528e+001;
n=90; farx(n+1)=7.756478e-001; foe(n+1)=1.211289e+000; krok(n+1)=1.417782e-002; ng(n+1)=2.273779e+002;
n=91; farx(n+1)=7.564623e-001; foe(n+1)=1.181038e+000; krok(n+1)=1.000091e-002; ng(n+1)=1.764781e+002;
n=92; farx(n+1)=7.142442e-001; foe(n+1)=1.143941e+000; krok(n+1)=4.086147e-002; ng(n+1)=1.383291e+002;
n=93; farx(n+1)=6.975190e-001; foe(n+1)=1.123704e+000; krok(n+1)=7.729026e-002; ng(n+1)=4.471189e+001;
n=94; farx(n+1)=6.628008e-001; foe(n+1)=1.092798e+000; krok(n+1)=8.343059e-002; ng(n+1)=1.127117e+002;
n=95; farx(n+1)=6.189158e-001; foe(n+1)=1.061218e+000; krok(n+1)=3.881594e-002; ng(n+1)=1.305396e+002;
n=96; farx(n+1)=5.930801e-001; foe(n+1)=1.019429e+000; krok(n+1)=7.413078e-002; ng(n+1)=1.461828e+002;
n=97; farx(n+1)=6.349001e-001; foe(n+1)=9.749207e-001; krok(n+1)=6.345766e-002; ng(n+1)=1.840848e+002;
n=98; farx(n+1)=6.597105e-001; foe(n+1)=9.499161e-001; krok(n+1)=1.015974e-001; ng(n+1)=8.193636e+001;
n=99; farx(n+1)=6.558273e-001; foe(n+1)=8.997025e-001; krok(n+1)=1.630613e-001; ng(n+1)=2.309349e+002;
n=100; farx(n+1)=6.453045e-001; foe(n+1)=8.458685e-001; krok(n+1)=1.600146e-001; ng(n+1)=2.126258e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
