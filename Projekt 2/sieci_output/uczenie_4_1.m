%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.814894e+003; foe(n+1)=4.858754e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.808461e+003; foe(n+1)=3.891880e+003; krok(n+1)=3.301592e-004; ng(n+1)=2.596069e+003;
n=2; farx(n+1)=8.068778e+002; foe(n+1)=8.537259e+002; krok(n+1)=5.320104e-003; ng(n+1)=1.121482e+003;
n=3; farx(n+1)=8.307500e+002; foe(n+1)=7.660664e+002; krok(n+1)=1.858107e-004; ng(n+1)=2.429177e+003;
n=4; farx(n+1)=9.388778e+002; foe(n+1)=6.514995e+002; krok(n+1)=7.937301e-004; ng(n+1)=2.699868e+003;
n=5; farx(n+1)=6.558149e+002; foe(n+1)=5.498977e+002; krok(n+1)=1.517951e-003; ng(n+1)=1.007069e+003;
n=6; farx(n+1)=3.958348e+002; foe(n+1)=4.266562e+002; krok(n+1)=3.692938e-003; ng(n+1)=1.505129e+003;
n=7; farx(n+1)=3.978808e+002; foe(n+1)=4.199696e+002; krok(n+1)=1.387045e-003; ng(n+1)=4.498266e+002;
n=8; farx(n+1)=3.711641e+002; foe(n+1)=4.137134e+002; krok(n+1)=4.114874e-003; ng(n+1)=2.865183e+002;
n=9; farx(n+1)=3.161257e+002; foe(n+1)=3.901056e+002; krok(n+1)=7.385876e-003; ng(n+1)=6.544530e+002;
n=10; farx(n+1)=3.048679e+002; foe(n+1)=3.862122e+002; krok(n+1)=8.118300e-004; ng(n+1)=7.240240e+002;
n=11; farx(n+1)=2.603122e+002; foe(n+1)=3.699346e+002; krok(n+1)=6.818297e-003; ng(n+1)=7.252680e+002;
n=12; farx(n+1)=1.910887e+002; foe(n+1)=3.349629e+002; krok(n+1)=2.855545e-003; ng(n+1)=1.262097e+003;
n=13; farx(n+1)=1.441524e+002; foe(n+1)=3.135408e+002; krok(n+1)=6.791933e-004; ng(n+1)=1.369992e+003;
n=14; farx(n+1)=1.362791e+002; foe(n+1)=3.011489e+002; krok(n+1)=4.209206e-004; ng(n+1)=2.059902e+003;
n=15; farx(n+1)=1.333368e+002; foe(n+1)=2.849542e+002; krok(n+1)=1.606904e-003; ng(n+1)=1.882707e+003;
n=16; farx(n+1)=1.350566e+002; foe(n+1)=2.835042e+002; krok(n+1)=3.604935e-004; ng(n+1)=3.180151e+003;
n=17; farx(n+1)=1.448789e+002; foe(n+1)=2.405886e+002; krok(n+1)=4.509885e-002; ng(n+1)=2.944053e+003;
n=18; farx(n+1)=1.457248e+002; foe(n+1)=2.393165e+002; krok(n+1)=3.812868e-004; ng(n+1)=2.519114e+003;
n=19; farx(n+1)=1.200269e+002; foe(n+1)=2.254803e+002; krok(n+1)=1.170476e-002; ng(n+1)=2.485868e+003;
n=20; farx(n+1)=9.408804e+001; foe(n+1)=2.136259e+002; krok(n+1)=3.500452e-002; ng(n+1)=1.850400e+003;
n=21; farx(n+1)=9.002146e+001; foe(n+1)=2.062557e+002; krok(n+1)=7.981571e-003; ng(n+1)=1.932401e+003;
n=22; farx(n+1)=6.526156e+001; foe(n+1)=1.612325e+002; krok(n+1)=1.482879e-001; ng(n+1)=1.967475e+003;
n=23; farx(n+1)=7.980201e+001; foe(n+1)=1.524288e+002; krok(n+1)=1.490738e-001; ng(n+1)=1.563634e+003;
n=24; farx(n+1)=4.958129e+001; foe(n+1)=1.330848e+002; krok(n+1)=2.483631e-001; ng(n+1)=1.116999e+003;
n=25; farx(n+1)=4.027551e+001; foe(n+1)=9.260711e+001; krok(n+1)=7.101669e-001; ng(n+1)=7.612938e+002;
%odnowa zmiennej metryki
n=26; farx(n+1)=4.046192e+001; foe(n+1)=9.184477e+001; krok(n+1)=7.573766e-006; ng(n+1)=1.110222e+003;
n=27; farx(n+1)=3.823459e+001; foe(n+1)=9.011014e+001; krok(n+1)=4.448722e-004; ng(n+1)=2.289524e+002;
n=28; farx(n+1)=3.454933e+001; foe(n+1)=8.351164e+001; krok(n+1)=9.821596e-004; ng(n+1)=2.388926e+002;
n=29; farx(n+1)=3.317440e+001; foe(n+1)=8.243311e+001; krok(n+1)=6.756225e-005; ng(n+1)=1.003911e+003;
n=30; farx(n+1)=2.893016e+001; foe(n+1)=7.528333e+001; krok(n+1)=3.758973e-004; ng(n+1)=2.419708e+003;
n=31; farx(n+1)=2.888495e+001; foe(n+1)=7.521740e+001; krok(n+1)=1.522728e-005; ng(n+1)=3.297724e+004;
n=32; farx(n+1)=2.085142e+001; foe(n+1)=6.317426e+001; krok(n+1)=1.876224e-003; ng(n+1)=4.461137e+004;
n=33; farx(n+1)=2.002438e+001; foe(n+1)=5.828744e+001; krok(n+1)=1.469853e-003; ng(n+1)=1.502747e+003;
n=34; farx(n+1)=1.325154e+001; foe(n+1)=4.977383e+001; krok(n+1)=2.304922e-003; ng(n+1)=7.146790e+002;
n=35; farx(n+1)=9.632634e+000; foe(n+1)=3.894409e+001; krok(n+1)=1.903300e-003; ng(n+1)=6.881733e+002;
n=36; farx(n+1)=8.154893e+000; foe(n+1)=3.468812e+001; krok(n+1)=2.594055e-004; ng(n+1)=3.077582e+003;
n=37; farx(n+1)=8.155974e+000; foe(n+1)=3.422919e+001; krok(n+1)=2.666638e-004; ng(n+1)=4.281881e+003;
n=38; farx(n+1)=6.944229e+000; foe(n+1)=2.715992e+001; krok(n+1)=1.019447e-002; ng(n+1)=4.317502e+003;
n=39; farx(n+1)=7.004938e+000; foe(n+1)=2.640360e+001; krok(n+1)=2.064885e-004; ng(n+1)=4.258955e+003;
n=40; farx(n+1)=7.036497e+000; foe(n+1)=2.389809e+001; krok(n+1)=4.134280e-002; ng(n+1)=3.833200e+003;
n=41; farx(n+1)=6.983051e+000; foe(n+1)=2.360851e+001; krok(n+1)=1.238474e-002; ng(n+1)=3.206042e+003;
n=42; farx(n+1)=6.983051e+000; foe(n+1)=2.360851e+001; krok(n+1)=1.527361e-014; ng(n+1)=3.253007e+003;
n=43; farx(n+1)=6.306953e+000; foe(n+1)=1.883649e+001; krok(n+1)=5.281352e-001; ng(n+1)=3.253007e+003;
n=44; farx(n+1)=5.671620e+000; foe(n+1)=1.537449e+001; krok(n+1)=5.945329e-002; ng(n+1)=1.802085e+003;
n=45; farx(n+1)=5.548073e+000; foe(n+1)=1.451218e+001; krok(n+1)=6.134237e-002; ng(n+1)=5.915050e+002;
n=46; farx(n+1)=5.462983e+000; foe(n+1)=1.374663e+001; krok(n+1)=5.034993e-002; ng(n+1)=9.137834e+002;
n=47; farx(n+1)=5.090815e+000; foe(n+1)=1.270763e+001; krok(n+1)=3.355407e-001; ng(n+1)=8.139189e+002;
n=48; farx(n+1)=4.181155e+000; foe(n+1)=1.073099e+001; krok(n+1)=1.435731e+000; ng(n+1)=3.756250e+002;
n=49; farx(n+1)=3.756109e+000; foe(n+1)=1.006080e+001; krok(n+1)=2.981443e-001; ng(n+1)=5.318905e+002;
n=50; farx(n+1)=3.558944e+000; foe(n+1)=9.357174e+000; krok(n+1)=4.872449e-001; ng(n+1)=7.256858e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=3.537053e+000; foe(n+1)=9.310213e+000; krok(n+1)=2.952022e-005; ng(n+1)=1.358095e+002;
n=52; farx(n+1)=3.518592e+000; foe(n+1)=9.293464e+000; krok(n+1)=2.491142e-005; ng(n+1)=9.078373e+001;
n=53; farx(n+1)=3.497430e+000; foe(n+1)=9.259931e+000; krok(n+1)=2.462151e-005; ng(n+1)=1.231146e+002;
n=54; farx(n+1)=3.426419e+000; foe(n+1)=9.047117e+000; krok(n+1)=3.778882e-004; ng(n+1)=7.308862e+001;
n=55; farx(n+1)=3.311872e+000; foe(n+1)=8.905933e+000; krok(n+1)=3.032635e-003; ng(n+1)=2.635424e+001;
n=56; farx(n+1)=3.169950e+000; foe(n+1)=8.773254e+000; krok(n+1)=4.011567e-003; ng(n+1)=4.108885e+001;
n=57; farx(n+1)=3.096826e+000; foe(n+1)=8.595727e+000; krok(n+1)=1.995002e-003; ng(n+1)=1.045706e+002;
n=58; farx(n+1)=2.973704e+000; foe(n+1)=8.396795e+000; krok(n+1)=8.204311e-003; ng(n+1)=2.563454e+002;
n=59; farx(n+1)=2.935848e+000; foe(n+1)=8.280062e+000; krok(n+1)=9.683310e-003; ng(n+1)=5.344163e+002;
n=60; farx(n+1)=2.873257e+000; foe(n+1)=8.172047e+000; krok(n+1)=1.525264e-002; ng(n+1)=4.229748e+002;
n=61; farx(n+1)=2.819497e+000; foe(n+1)=8.124221e+000; krok(n+1)=1.305321e-002; ng(n+1)=4.216528e+002;
n=62; farx(n+1)=2.725132e+000; foe(n+1)=8.014601e+000; krok(n+1)=5.671128e-002; ng(n+1)=4.548933e+002;
n=63; farx(n+1)=2.571464e+000; foe(n+1)=7.707314e+000; krok(n+1)=5.157588e-002; ng(n+1)=3.725969e+002;
n=64; farx(n+1)=2.233385e+000; foe(n+1)=7.452196e+000; krok(n+1)=6.823373e-002; ng(n+1)=2.389351e+002;
n=65; farx(n+1)=2.060871e+000; foe(n+1)=7.198816e+000; krok(n+1)=7.581460e-002; ng(n+1)=4.138845e+002;
n=66; farx(n+1)=2.060717e+000; foe(n+1)=7.081450e+000; krok(n+1)=7.317146e-002; ng(n+1)=1.491806e+002;
n=67; farx(n+1)=1.905091e+000; foe(n+1)=6.319027e+000; krok(n+1)=8.612093e-001; ng(n+1)=1.215353e+002;
n=68; farx(n+1)=1.869671e+000; foe(n+1)=6.155555e+000; krok(n+1)=1.217616e-001; ng(n+1)=5.889635e+002;
n=69; farx(n+1)=1.867104e+000; foe(n+1)=5.937925e+000; krok(n+1)=1.871606e-001; ng(n+1)=4.124939e+002;
n=70; farx(n+1)=1.761354e+000; foe(n+1)=5.679333e+000; krok(n+1)=2.818076e-001; ng(n+1)=2.992397e+002;
n=71; farx(n+1)=1.723825e+000; foe(n+1)=5.537339e+000; krok(n+1)=1.677703e-001; ng(n+1)=2.873318e+002;
n=72; farx(n+1)=1.813034e+000; foe(n+1)=5.395195e+000; krok(n+1)=5.631821e-001; ng(n+1)=2.061914e+002;
n=73; farx(n+1)=1.774598e+000; foe(n+1)=5.266928e+000; krok(n+1)=3.614383e-001; ng(n+1)=1.600926e+002;
n=74; farx(n+1)=1.632362e+000; foe(n+1)=5.158032e+000; krok(n+1)=1.453256e-001; ng(n+1)=1.677145e+002;
n=75; farx(n+1)=1.480569e+000; foe(n+1)=4.985644e+000; krok(n+1)=1.845727e-001; ng(n+1)=3.261060e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.473088e+000; foe(n+1)=4.936105e+000; krok(n+1)=2.394869e-006; ng(n+1)=3.779438e+002;
n=77; farx(n+1)=1.465362e+000; foe(n+1)=4.905782e+000; krok(n+1)=1.928098e-005; ng(n+1)=1.387412e+002;
n=78; farx(n+1)=1.458683e+000; foe(n+1)=4.893125e+000; krok(n+1)=1.351677e-005; ng(n+1)=1.083269e+002;
n=79; farx(n+1)=1.445621e+000; foe(n+1)=4.865550e+000; krok(n+1)=5.405350e-004; ng(n+1)=2.329666e+001;
n=80; farx(n+1)=1.433019e+000; foe(n+1)=4.841273e+000; krok(n+1)=1.022424e-003; ng(n+1)=1.787513e+001;
n=81; farx(n+1)=1.438098e+000; foe(n+1)=4.830023e+000; krok(n+1)=1.021844e-003; ng(n+1)=1.187909e+001;
n=82; farx(n+1)=1.408477e+000; foe(n+1)=4.785919e+000; krok(n+1)=1.037360e-002; ng(n+1)=1.400425e+001;
n=83; farx(n+1)=1.353276e+000; foe(n+1)=4.731575e+000; krok(n+1)=4.606611e-003; ng(n+1)=8.273259e+001;
n=84; farx(n+1)=1.312627e+000; foe(n+1)=4.699765e+000; krok(n+1)=4.207904e-003; ng(n+1)=3.079181e+002;
n=85; farx(n+1)=1.264093e+000; foe(n+1)=4.655628e+000; krok(n+1)=3.229291e-002; ng(n+1)=5.028936e+002;
n=86; farx(n+1)=1.266176e+000; foe(n+1)=4.639142e+000; krok(n+1)=5.113139e-002; ng(n+1)=4.583023e+002;
n=87; farx(n+1)=1.240633e+000; foe(n+1)=4.615328e+000; krok(n+1)=6.503090e-002; ng(n+1)=3.459430e+002;
n=88; farx(n+1)=1.195775e+000; foe(n+1)=4.578036e+000; krok(n+1)=6.736795e-002; ng(n+1)=4.011333e+002;
n=89; farx(n+1)=1.142455e+000; foe(n+1)=4.484255e+000; krok(n+1)=9.327289e-002; ng(n+1)=2.844640e+002;
n=90; farx(n+1)=1.142197e+000; foe(n+1)=4.460883e+000; krok(n+1)=4.663645e-002; ng(n+1)=2.205023e+002;
n=91; farx(n+1)=1.260917e+000; foe(n+1)=4.420428e+000; krok(n+1)=1.237945e-001; ng(n+1)=9.807177e+001;
n=92; farx(n+1)=1.241465e+000; foe(n+1)=4.347692e+000; krok(n+1)=1.076512e-001; ng(n+1)=1.141783e+002;
n=93; farx(n+1)=1.235249e+000; foe(n+1)=4.301812e+000; krok(n+1)=1.400181e-001; ng(n+1)=3.703049e+002;
n=94; farx(n+1)=1.180274e+000; foe(n+1)=4.151892e+000; krok(n+1)=5.045340e-001; ng(n+1)=3.248581e+002;
n=95; farx(n+1)=1.176611e+000; foe(n+1)=4.103222e+000; krok(n+1)=9.048247e-002; ng(n+1)=4.880269e+002;
n=96; farx(n+1)=1.095462e+000; foe(n+1)=3.983228e+000; krok(n+1)=2.800362e-001; ng(n+1)=2.890526e+002;
n=97; farx(n+1)=9.645848e-001; foe(n+1)=3.879992e+000; krok(n+1)=1.348264e-001; ng(n+1)=3.379585e+002;
n=98; farx(n+1)=8.860149e-001; foe(n+1)=3.714287e+000; krok(n+1)=2.710777e-001; ng(n+1)=1.517952e+002;
n=99; farx(n+1)=8.552760e-001; foe(n+1)=3.590899e+000; krok(n+1)=2.615920e-001; ng(n+1)=6.536944e+002;
n=100; farx(n+1)=8.777454e-001; foe(n+1)=3.474671e+000; krok(n+1)=2.730638e-001; ng(n+1)=7.775593e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
