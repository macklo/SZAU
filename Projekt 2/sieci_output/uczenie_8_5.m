%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.107805e+003; foe(n+1)=3.996462e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.856696e+003; foe(n+1)=2.882846e+003; krok(n+1)=5.109218e-004; ng(n+1)=4.003925e+003;
n=2; farx(n+1)=1.039749e+003; foe(n+1)=9.327517e+002; krok(n+1)=2.530456e-003; ng(n+1)=2.633316e+003;
n=3; farx(n+1)=1.153562e+003; foe(n+1)=5.935502e+002; krok(n+1)=2.288097e-004; ng(n+1)=6.529680e+003;
n=4; farx(n+1)=1.456554e+003; foe(n+1)=5.370953e+002; krok(n+1)=8.523876e-004; ng(n+1)=3.996172e+003;
n=5; farx(n+1)=9.392302e+002; foe(n+1)=4.625547e+002; krok(n+1)=7.251821e-003; ng(n+1)=6.933188e+002;
n=6; farx(n+1)=3.351185e+002; foe(n+1)=3.440033e+002; krok(n+1)=6.118798e-003; ng(n+1)=1.557593e+003;
n=7; farx(n+1)=3.135200e+002; foe(n+1)=3.421816e+002; krok(n+1)=9.772716e-006; ng(n+1)=2.822083e+003;
n=8; farx(n+1)=2.073702e+002; foe(n+1)=3.157245e+002; krok(n+1)=4.652564e-004; ng(n+1)=3.446152e+003;
n=9; farx(n+1)=1.779734e+002; foe(n+1)=3.100971e+002; krok(n+1)=1.417047e-004; ng(n+1)=4.317859e+003;
n=10; farx(n+1)=1.723509e+002; foe(n+1)=3.066006e+002; krok(n+1)=2.201622e-003; ng(n+1)=4.183050e+003;
n=11; farx(n+1)=1.480054e+002; foe(n+1)=2.915870e+002; krok(n+1)=2.716773e-003; ng(n+1)=4.205091e+003;
n=12; farx(n+1)=1.412340e+002; foe(n+1)=2.849169e+002; krok(n+1)=2.075244e-003; ng(n+1)=1.777141e+003;
n=13; farx(n+1)=1.437105e+002; foe(n+1)=2.746384e+002; krok(n+1)=2.446921e-003; ng(n+1)=1.340564e+003;
n=14; farx(n+1)=1.401359e+002; foe(n+1)=2.659767e+002; krok(n+1)=3.505914e-003; ng(n+1)=1.155282e+003;
n=15; farx(n+1)=1.123686e+002; foe(n+1)=2.463197e+002; krok(n+1)=1.073668e-002; ng(n+1)=2.138469e+003;
n=16; farx(n+1)=1.084123e+002; foe(n+1)=2.368673e+002; krok(n+1)=3.643472e-004; ng(n+1)=3.406700e+003;
n=17; farx(n+1)=9.970181e+001; foe(n+1)=2.314286e+002; krok(n+1)=2.988526e-003; ng(n+1)=1.920660e+003;
n=18; farx(n+1)=9.484281e+001; foe(n+1)=2.265083e+002; krok(n+1)=1.013887e-003; ng(n+1)=1.635560e+003;
n=19; farx(n+1)=6.884477e+001; foe(n+1)=2.059946e+002; krok(n+1)=6.193383e-003; ng(n+1)=1.305025e+003;
n=20; farx(n+1)=5.577497e+001; foe(n+1)=1.966345e+002; krok(n+1)=2.630324e-004; ng(n+1)=3.567497e+003;
n=21; farx(n+1)=2.850282e+001; foe(n+1)=1.834158e+002; krok(n+1)=4.036613e-003; ng(n+1)=5.298724e+003;
n=22; farx(n+1)=2.716723e+001; foe(n+1)=1.294673e+002; krok(n+1)=1.200037e-004; ng(n+1)=1.149245e+004;
n=23; farx(n+1)=2.808033e+001; foe(n+1)=1.286654e+002; krok(n+1)=1.115447e-005; ng(n+1)=9.910484e+003;
n=24; farx(n+1)=2.688264e+001; foe(n+1)=1.275742e+002; krok(n+1)=2.261309e-003; ng(n+1)=8.988653e+003;
n=25; farx(n+1)=2.888979e+001; foe(n+1)=1.233969e+002; krok(n+1)=3.146871e-003; ng(n+1)=9.385010e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=2.169009e+001; foe(n+1)=9.201431e+001; krok(n+1)=1.882199e-005; ng(n+1)=7.293744e+003;
n=27; farx(n+1)=1.289530e+001; foe(n+1)=7.108081e+001; krok(n+1)=2.877521e-005; ng(n+1)=3.103760e+003;
n=28; farx(n+1)=1.148447e+001; foe(n+1)=4.723229e+001; krok(n+1)=4.449379e-004; ng(n+1)=1.008368e+003;
n=29; farx(n+1)=9.777062e+000; foe(n+1)=4.282936e+001; krok(n+1)=1.895397e-004; ng(n+1)=1.011764e+003;
n=30; farx(n+1)=9.163361e+000; foe(n+1)=3.344071e+001; krok(n+1)=8.861138e-004; ng(n+1)=1.602370e+003;
n=31; farx(n+1)=7.776074e+000; foe(n+1)=2.648311e+001; krok(n+1)=1.338974e-003; ng(n+1)=1.839498e+003;
n=32; farx(n+1)=6.263429e+000; foe(n+1)=1.622037e+001; krok(n+1)=3.910756e-003; ng(n+1)=1.122940e+003;
n=33; farx(n+1)=6.059143e+000; foe(n+1)=1.492371e+001; krok(n+1)=7.459822e-004; ng(n+1)=1.013199e+003;
n=34; farx(n+1)=5.907774e+000; foe(n+1)=1.361864e+001; krok(n+1)=2.973498e-003; ng(n+1)=4.597356e+002;
n=35; farx(n+1)=5.684648e+000; foe(n+1)=1.232670e+001; krok(n+1)=7.162643e-003; ng(n+1)=9.625995e+002;
n=36; farx(n+1)=5.180664e+000; foe(n+1)=1.102488e+001; krok(n+1)=3.485222e-003; ng(n+1)=5.117014e+002;
n=37; farx(n+1)=3.930146e+000; foe(n+1)=8.661528e+000; krok(n+1)=3.022871e-002; ng(n+1)=2.758058e+002;
n=38; farx(n+1)=3.721642e+000; foe(n+1)=8.347940e+000; krok(n+1)=1.043660e-003; ng(n+1)=5.675803e+002;
n=39; farx(n+1)=3.191218e+000; foe(n+1)=7.566253e+000; krok(n+1)=1.007731e-002; ng(n+1)=2.400391e+002;
n=40; farx(n+1)=2.736721e+000; foe(n+1)=6.967394e+000; krok(n+1)=8.053672e-003; ng(n+1)=2.743395e+002;
n=41; farx(n+1)=2.595445e+000; foe(n+1)=6.701287e+000; krok(n+1)=2.531421e-003; ng(n+1)=3.693567e+002;
n=42; farx(n+1)=2.526119e+000; foe(n+1)=6.339674e+000; krok(n+1)=5.637357e-003; ng(n+1)=4.918235e+002;
n=43; farx(n+1)=2.484343e+000; foe(n+1)=6.105143e+000; krok(n+1)=5.233259e-003; ng(n+1)=2.387514e+002;
n=44; farx(n+1)=2.455394e+000; foe(n+1)=5.901127e+000; krok(n+1)=6.753788e-003; ng(n+1)=2.843118e+002;
n=45; farx(n+1)=2.258522e+000; foe(n+1)=5.238497e+000; krok(n+1)=3.618094e-002; ng(n+1)=5.649228e+002;
n=46; farx(n+1)=2.203818e+000; foe(n+1)=5.042636e+000; krok(n+1)=6.764823e-003; ng(n+1)=4.275422e+002;
n=47; farx(n+1)=2.173322e+000; foe(n+1)=4.926611e+000; krok(n+1)=1.409339e-003; ng(n+1)=6.521997e+002;
n=48; farx(n+1)=2.136651e+000; foe(n+1)=4.661623e+000; krok(n+1)=1.004449e-002; ng(n+1)=4.360362e+002;
n=49; farx(n+1)=2.094009e+000; foe(n+1)=4.300586e+000; krok(n+1)=6.852302e-002; ng(n+1)=3.625961e+002;
n=50; farx(n+1)=2.062424e+000; foe(n+1)=4.095274e+000; krok(n+1)=5.919426e-002; ng(n+1)=5.325610e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=2.044415e+000; foe(n+1)=3.976586e+000; krok(n+1)=1.395478e-005; ng(n+1)=4.378949e+002;
n=52; farx(n+1)=2.039436e+000; foe(n+1)=3.955230e+000; krok(n+1)=7.204914e-006; ng(n+1)=2.455242e+002;
n=53; farx(n+1)=2.037592e+000; foe(n+1)=3.932669e+000; krok(n+1)=1.212530e-004; ng(n+1)=6.754802e+001;
n=54; farx(n+1)=2.039079e+000; foe(n+1)=3.886156e+000; krok(n+1)=6.791933e-004; ng(n+1)=4.089771e+001;
n=55; farx(n+1)=2.018348e+000; foe(n+1)=3.795824e+000; krok(n+1)=6.500859e-004; ng(n+1)=6.058715e+001;
n=56; farx(n+1)=1.999782e+000; foe(n+1)=3.747864e+000; krok(n+1)=4.994367e-004; ng(n+1)=4.105874e+001;
n=57; farx(n+1)=1.941561e+000; foe(n+1)=3.475565e+000; krok(n+1)=9.222377e-003; ng(n+1)=3.557706e+001;
n=58; farx(n+1)=1.872223e+000; foe(n+1)=3.323205e+000; krok(n+1)=9.793864e-003; ng(n+1)=1.130330e+002;
n=59; farx(n+1)=1.834839e+000; foe(n+1)=3.192038e+000; krok(n+1)=6.669234e-003; ng(n+1)=2.423098e+002;
n=60; farx(n+1)=1.739802e+000; foe(n+1)=3.020861e+000; krok(n+1)=1.314565e-002; ng(n+1)=2.066130e+002;
n=61; farx(n+1)=1.712882e+000; foe(n+1)=2.977944e+000; krok(n+1)=5.557598e-003; ng(n+1)=1.285421e+002;
n=62; farx(n+1)=1.648713e+000; foe(n+1)=2.886746e+000; krok(n+1)=2.055621e-002; ng(n+1)=3.183684e+002;
n=63; farx(n+1)=1.625388e+000; foe(n+1)=2.798026e+000; krok(n+1)=1.406421e-002; ng(n+1)=2.629235e+002;
n=64; farx(n+1)=1.629735e+000; foe(n+1)=2.687440e+000; krok(n+1)=2.047906e-002; ng(n+1)=3.456503e+002;
n=65; farx(n+1)=1.624480e+000; foe(n+1)=2.606075e+000; krok(n+1)=7.782134e-003; ng(n+1)=5.427410e+002;
n=66; farx(n+1)=1.591984e+000; foe(n+1)=2.517477e+000; krok(n+1)=4.858193e-002; ng(n+1)=1.815136e+002;
n=67; farx(n+1)=1.564925e+000; foe(n+1)=2.367951e+000; krok(n+1)=4.294671e-002; ng(n+1)=3.479353e+002;
n=68; farx(n+1)=1.558616e+000; foe(n+1)=2.227935e+000; krok(n+1)=6.158158e-002; ng(n+1)=3.087900e+002;
n=69; farx(n+1)=1.520562e+000; foe(n+1)=2.160419e+000; krok(n+1)=3.411686e-002; ng(n+1)=3.531442e+002;
n=70; farx(n+1)=1.445583e+000; foe(n+1)=2.092844e+000; krok(n+1)=4.524124e-002; ng(n+1)=2.646343e+002;
n=71; farx(n+1)=1.298168e+000; foe(n+1)=1.944426e+000; krok(n+1)=5.684082e-002; ng(n+1)=5.635222e+002;
n=72; farx(n+1)=1.233878e+000; foe(n+1)=1.886718e+000; krok(n+1)=2.617232e-002; ng(n+1)=2.941148e+002;
n=73; farx(n+1)=1.185920e+000; foe(n+1)=1.832088e+000; krok(n+1)=4.834060e-002; ng(n+1)=4.284112e+002;
n=74; farx(n+1)=1.134826e+000; foe(n+1)=1.747843e+000; krok(n+1)=6.306675e-002; ng(n+1)=1.368730e+002;
n=75; farx(n+1)=1.112146e+000; foe(n+1)=1.710991e+000; krok(n+1)=1.081188e-001; ng(n+1)=2.294861e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.111938e+000; foe(n+1)=1.703083e+000; krok(n+1)=4.858372e-006; ng(n+1)=1.851612e+002;
n=77; farx(n+1)=1.112285e+000; foe(n+1)=1.694434e+000; krok(n+1)=6.307179e-006; ng(n+1)=1.564786e+002;
n=78; farx(n+1)=1.109984e+000; foe(n+1)=1.663115e+000; krok(n+1)=3.584012e-005; ng(n+1)=1.378586e+002;
n=79; farx(n+1)=1.113742e+000; foe(n+1)=1.644322e+000; krok(n+1)=8.763006e-005; ng(n+1)=6.112065e+001;
n=80; farx(n+1)=1.114472e+000; foe(n+1)=1.641299e+000; krok(n+1)=2.084136e-004; ng(n+1)=2.557324e+001;
n=81; farx(n+1)=1.113292e+000; foe(n+1)=1.635929e+000; krok(n+1)=6.218716e-004; ng(n+1)=2.216243e+001;
n=82; farx(n+1)=1.112682e+000; foe(n+1)=1.623649e+000; krok(n+1)=6.764823e-003; ng(n+1)=1.343809e+001;
n=83; farx(n+1)=1.112394e+000; foe(n+1)=1.616048e+000; krok(n+1)=7.504897e-003; ng(n+1)=1.454962e+001;
n=84; farx(n+1)=1.114525e+000; foe(n+1)=1.608881e+000; krok(n+1)=5.390749e-003; ng(n+1)=2.018696e+001;
n=85; farx(n+1)=1.101938e+000; foe(n+1)=1.589361e+000; krok(n+1)=2.307159e-002; ng(n+1)=1.994697e+001;
n=86; farx(n+1)=1.068499e+000; foe(n+1)=1.544201e+000; krok(n+1)=1.460372e-001; ng(n+1)=7.611658e+001;
n=87; farx(n+1)=1.058028e+000; foe(n+1)=1.533363e+000; krok(n+1)=1.438322e-002; ng(n+1)=1.308112e+002;
n=88; farx(n+1)=1.040331e+000; foe(n+1)=1.500286e+000; krok(n+1)=3.722882e-002; ng(n+1)=1.332673e+002;
n=89; farx(n+1)=1.030902e+000; foe(n+1)=1.476537e+000; krok(n+1)=1.729712e-002; ng(n+1)=1.748997e+002;
n=90; farx(n+1)=1.037677e+000; foe(n+1)=1.437953e+000; krok(n+1)=3.535818e-002; ng(n+1)=1.273015e+002;
n=91; farx(n+1)=1.043101e+000; foe(n+1)=1.412393e+000; krok(n+1)=6.607631e-003; ng(n+1)=2.502025e+002;
n=92; farx(n+1)=1.056098e+000; foe(n+1)=1.390197e+000; krok(n+1)=2.477353e-002; ng(n+1)=8.612064e+001;
n=93; farx(n+1)=1.070906e+000; foe(n+1)=1.372261e+000; krok(n+1)=2.243330e-002; ng(n+1)=1.948609e+002;
n=94; farx(n+1)=1.080744e+000; foe(n+1)=1.351518e+000; krok(n+1)=1.977759e-002; ng(n+1)=1.501709e+002;
n=95; farx(n+1)=1.092541e+000; foe(n+1)=1.327115e+000; krok(n+1)=1.576669e-002; ng(n+1)=1.452465e+002;
n=96; farx(n+1)=1.101340e+000; foe(n+1)=1.287862e+000; krok(n+1)=5.411858e-002; ng(n+1)=1.016106e+002;
n=97; farx(n+1)=1.100997e+000; foe(n+1)=1.266916e+000; krok(n+1)=4.879939e-002; ng(n+1)=1.554358e+002;
n=98; farx(n+1)=1.109553e+000; foe(n+1)=1.231184e+000; krok(n+1)=6.407656e-002; ng(n+1)=2.344633e+002;
n=99; farx(n+1)=1.116918e+000; foe(n+1)=1.198511e+000; krok(n+1)=8.035594e-002; ng(n+1)=2.293927e+002;
n=100; farx(n+1)=1.098021e+000; foe(n+1)=1.155750e+000; krok(n+1)=1.139041e-001; ng(n+1)=1.675360e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
