%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.048753e+003; foe(n+1)=4.054766e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.007507e+003; foe(n+1)=3.010071e+003; krok(n+1)=4.988661e-004; ng(n+1)=2.099032e+003;
n=2; farx(n+1)=1.042132e+003; foe(n+1)=1.066453e+003; krok(n+1)=6.450159e-003; ng(n+1)=4.937793e+002;
n=3; farx(n+1)=1.071309e+003; foe(n+1)=9.895782e+002; krok(n+1)=1.976919e-005; ng(n+1)=5.218087e+003;
n=4; farx(n+1)=3.129968e+003; foe(n+1)=7.493082e+002; krok(n+1)=2.678603e-003; ng(n+1)=5.628008e+003;
n=5; farx(n+1)=2.408057e+003; foe(n+1)=5.473608e+002; krok(n+1)=3.771789e-003; ng(n+1)=1.725259e+003;
n=6; farx(n+1)=2.029472e+003; foe(n+1)=5.139620e+002; krok(n+1)=3.776610e-004; ng(n+1)=1.135496e+003;
n=7; farx(n+1)=1.540918e+003; foe(n+1)=4.921713e+002; krok(n+1)=2.659087e-003; ng(n+1)=7.570318e+002;
n=8; farx(n+1)=1.061120e+003; foe(n+1)=4.192855e+002; krok(n+1)=1.564230e-003; ng(n+1)=1.485316e+003;
n=9; farx(n+1)=9.986870e+002; foe(n+1)=4.055876e+002; krok(n+1)=5.637357e-003; ng(n+1)=7.199781e+002;
n=10; farx(n+1)=1.001779e+003; foe(n+1)=4.045730e+002; krok(n+1)=6.192369e-003; ng(n+1)=1.143295e+003;
n=11; farx(n+1)=1.001861e+003; foe(n+1)=4.019359e+002; krok(n+1)=4.289834e-003; ng(n+1)=1.030683e+003;
n=12; farx(n+1)=7.976937e+002; foe(n+1)=3.928589e+002; krok(n+1)=6.004768e-002; ng(n+1)=7.152279e+002;
n=13; farx(n+1)=4.285433e+002; foe(n+1)=3.547077e+002; krok(n+1)=2.192961e-001; ng(n+1)=2.570561e+002;
n=14; farx(n+1)=3.822274e+002; foe(n+1)=3.507743e+002; krok(n+1)=1.027490e-002; ng(n+1)=5.190621e+003;
n=15; farx(n+1)=3.421543e+002; foe(n+1)=3.411655e+002; krok(n+1)=4.614317e-002; ng(n+1)=7.795545e+003;
n=16; farx(n+1)=1.933759e+002; foe(n+1)=3.099284e+002; krok(n+1)=1.742650e-001; ng(n+1)=1.010190e+004;
n=17; farx(n+1)=1.401956e+002; foe(n+1)=2.777842e+002; krok(n+1)=7.255280e-002; ng(n+1)=2.086296e+004;
n=18; farx(n+1)=1.360481e+002; foe(n+1)=2.327686e+002; krok(n+1)=2.418674e-001; ng(n+1)=1.486749e+004;
n=19; farx(n+1)=1.180488e+002; foe(n+1)=2.013168e+002; krok(n+1)=1.630613e-001; ng(n+1)=7.533716e+003;
n=20; farx(n+1)=8.830854e+001; foe(n+1)=1.700662e+002; krok(n+1)=3.015923e-001; ng(n+1)=1.879048e+003;
n=21; farx(n+1)=2.809995e+001; foe(n+1)=1.402544e+002; krok(n+1)=6.322648e-001; ng(n+1)=1.386348e+003;
n=22; farx(n+1)=1.823208e+001; foe(n+1)=1.252877e+002; krok(n+1)=3.923690e-001; ng(n+1)=1.076137e+003;
n=23; farx(n+1)=1.240062e+001; foe(n+1)=1.096456e+002; krok(n+1)=9.270274e-001; ng(n+1)=3.000085e+003;
n=24; farx(n+1)=8.694328e+000; foe(n+1)=9.016549e+001; krok(n+1)=4.288049e-001; ng(n+1)=3.084872e+003;
n=25; farx(n+1)=9.380780e+000; foe(n+1)=8.091199e+001; krok(n+1)=4.837349e-001; ng(n+1)=1.020663e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=9.189150e+000; foe(n+1)=8.015884e+001; krok(n+1)=1.337356e-006; ng(n+1)=1.629303e+003;
n=27; farx(n+1)=8.918234e+000; foe(n+1)=7.952946e+001; krok(n+1)=2.077412e-005; ng(n+1)=3.615750e+002;
n=28; farx(n+1)=8.508297e+000; foe(n+1)=7.862437e+001; krok(n+1)=1.773994e-005; ng(n+1)=5.332377e+002;
n=29; farx(n+1)=6.922548e+000; foe(n+1)=7.588768e+001; krok(n+1)=8.202462e-003; ng(n+1)=1.553959e+002;
n=30; farx(n+1)=7.138238e+000; foe(n+1)=7.537870e+001; krok(n+1)=5.792131e-003; ng(n+1)=1.620683e+002;
n=31; farx(n+1)=6.733950e+000; foe(n+1)=7.515801e+001; krok(n+1)=1.367655e-002; ng(n+1)=1.910322e+002;
n=32; farx(n+1)=5.793721e+000; foe(n+1)=7.365255e+001; krok(n+1)=4.614317e-002; ng(n+1)=2.374087e+002;
n=33; farx(n+1)=6.682928e+000; foe(n+1)=7.333516e+001; krok(n+1)=3.920662e-002; ng(n+1)=1.666279e+002;
n=34; farx(n+1)=7.100410e+000; foe(n+1)=7.302529e+001; krok(n+1)=3.923690e-001; ng(n+1)=1.712383e+002;
n=35; farx(n+1)=6.177488e+000; foe(n+1)=7.277156e+001; krok(n+1)=1.981558e-001; ng(n+1)=2.629148e+002;
n=36; farx(n+1)=5.952332e+000; foe(n+1)=7.247432e+001; krok(n+1)=6.614848e-001; ng(n+1)=4.148060e+002;
n=37; farx(n+1)=5.816614e+000; foe(n+1)=7.231879e+001; krok(n+1)=1.008497e-001; ng(n+1)=4.475599e+002;
n=38; farx(n+1)=6.125902e+000; foe(n+1)=7.193338e+001; krok(n+1)=5.180995e-001; ng(n+1)=7.222438e+002;
n=39; farx(n+1)=6.317442e+000; foe(n+1)=7.180726e+001; krok(n+1)=2.343012e-001; ng(n+1)=2.207901e+002;
n=40; farx(n+1)=6.456274e+000; foe(n+1)=7.172012e+001; krok(n+1)=8.536620e-001; ng(n+1)=2.797407e+002;
n=41; farx(n+1)=6.552806e+000; foe(n+1)=7.159687e+001; krok(n+1)=7.970783e-001; ng(n+1)=2.670089e+002;
n=42; farx(n+1)=6.519887e+000; foe(n+1)=7.152408e+001; krok(n+1)=4.936410e-001; ng(n+1)=4.716592e+002;
n=43; farx(n+1)=6.258988e+000; foe(n+1)=7.143574e+001; krok(n+1)=4.635137e-001; ng(n+1)=4.342558e+002;
n=44; farx(n+1)=6.061218e+000; foe(n+1)=7.138937e+001; krok(n+1)=6.809733e-001; ng(n+1)=7.686602e+001;
n=45; farx(n+1)=6.082960e+000; foe(n+1)=7.134762e+001; krok(n+1)=1.447680e+000; ng(n+1)=2.031665e+002;
n=46; farx(n+1)=5.950840e+000; foe(n+1)=7.126276e+001; krok(n+1)=1.591497e+000; ng(n+1)=7.907696e+001;
n=47; farx(n+1)=5.730424e+000; foe(n+1)=7.117994e+001; krok(n+1)=5.870903e-001; ng(n+1)=1.570697e+002;
n=48; farx(n+1)=5.523648e+000; foe(n+1)=7.099112e+001; krok(n+1)=7.487569e-001; ng(n+1)=3.711137e+002;
n=49; farx(n+1)=5.341459e+000; foe(n+1)=7.086228e+001; krok(n+1)=4.275922e-001; ng(n+1)=3.338514e+002;
n=50; farx(n+1)=5.092436e+000; foe(n+1)=7.079282e+001; krok(n+1)=3.791860e-001; ng(n+1)=4.310242e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=5.102340e+000; foe(n+1)=7.077075e+001; krok(n+1)=6.820257e-007; ng(n+1)=4.140488e+002;
n=52; farx(n+1)=5.093614e+000; foe(n+1)=7.075891e+001; krok(n+1)=2.134892e-006; ng(n+1)=1.570243e+002;
n=53; farx(n+1)=5.077607e+000; foe(n+1)=7.075121e+001; krok(n+1)=2.092147e-005; ng(n+1)=4.009451e+001;
n=54; farx(n+1)=5.016653e+000; foe(n+1)=7.069632e+001; krok(n+1)=6.920879e-003; ng(n+1)=9.006353e+000;
n=55; farx(n+1)=5.126626e+000; foe(n+1)=7.066493e+001; krok(n+1)=4.708341e-003; ng(n+1)=1.725900e+001;
n=56; farx(n+1)=5.399987e+000; foe(n+1)=7.057977e+001; krok(n+1)=9.454744e-003; ng(n+1)=4.951476e+001;
n=57; farx(n+1)=5.272932e+000; foe(n+1)=7.055008e+001; krok(n+1)=1.308616e-002; ng(n+1)=3.264159e+002;
n=58; farx(n+1)=5.155451e+000; foe(n+1)=7.050726e+001; krok(n+1)=6.225707e-002; ng(n+1)=4.799469e+002;
n=59; farx(n+1)=5.223456e+000; foe(n+1)=7.047793e+001; krok(n+1)=5.623415e-002; ng(n+1)=4.944731e+002;
n=60; farx(n+1)=5.498035e+000; foe(n+1)=7.040705e+001; krok(n+1)=3.091610e-001; ng(n+1)=5.532925e+002;
n=61; farx(n+1)=5.202970e+000; foe(n+1)=7.019890e+001; krok(n+1)=1.686145e-001; ng(n+1)=3.684362e+002;
n=62; farx(n+1)=5.195907e+000; foe(n+1)=6.996760e+001; krok(n+1)=8.198208e-001; ng(n+1)=2.713657e+002;
n=63; farx(n+1)=5.396579e+000; foe(n+1)=6.963337e+001; krok(n+1)=1.420334e+000; ng(n+1)=6.784010e+002;
n=64; farx(n+1)=5.421024e+000; foe(n+1)=6.955724e+001; krok(n+1)=2.230612e-001; ng(n+1)=7.282754e+002;
n=65; farx(n+1)=5.273896e+000; foe(n+1)=6.946816e+001; krok(n+1)=5.467184e-001; ng(n+1)=3.332929e+002;
n=66; farx(n+1)=5.256956e+000; foe(n+1)=6.931685e+001; krok(n+1)=5.766634e-001; ng(n+1)=1.124979e+002;
n=67; farx(n+1)=5.320660e+000; foe(n+1)=6.899198e+001; krok(n+1)=1.668247e+000; ng(n+1)=6.536216e+002;
n=68; farx(n+1)=5.463315e+000; foe(n+1)=6.873050e+001; krok(n+1)=7.310195e-001; ng(n+1)=1.081663e+003;
n=69; farx(n+1)=5.391822e+000; foe(n+1)=6.849591e+001; krok(n+1)=5.249576e-001; ng(n+1)=2.678691e+002;
n=70; farx(n+1)=5.107498e+000; foe(n+1)=6.815725e+001; krok(n+1)=5.242150e-001; ng(n+1)=6.874107e+002;
n=71; farx(n+1)=4.849531e+000; foe(n+1)=6.798794e+001; krok(n+1)=3.614383e-001; ng(n+1)=4.038826e+002;
n=72; farx(n+1)=4.484045e+000; foe(n+1)=6.766619e+001; krok(n+1)=3.268918e-001; ng(n+1)=8.066776e+002;
n=73; farx(n+1)=4.205324e+000; foe(n+1)=6.746455e+001; krok(n+1)=2.073757e-001; ng(n+1)=5.463795e+002;
n=74; farx(n+1)=3.987687e+000; foe(n+1)=6.719689e+001; krok(n+1)=3.311156e-001; ng(n+1)=4.834712e+002;
n=75; farx(n+1)=3.991045e+000; foe(n+1)=6.686141e+001; krok(n+1)=2.863617e-001; ng(n+1)=1.478601e+003;
%odnowa zmiennej metryki
n=76; farx(n+1)=4.006848e+000; foe(n+1)=6.674169e+001; krok(n+1)=2.843325e-007; ng(n+1)=1.317881e+003;
n=77; farx(n+1)=3.992897e+000; foe(n+1)=6.670162e+001; krok(n+1)=6.463557e-007; ng(n+1)=5.699166e+002;
n=78; farx(n+1)=3.990785e+000; foe(n+1)=6.670082e+001; krok(n+1)=2.020969e-005; ng(n+1)=1.957414e+001;
n=79; farx(n+1)=3.997930e+000; foe(n+1)=6.667878e+001; krok(n+1)=3.185772e-004; ng(n+1)=2.448129e+001;
n=80; farx(n+1)=4.303633e+000; foe(n+1)=6.656463e+001; krok(n+1)=4.953534e-003; ng(n+1)=2.754746e+001;
n=81; farx(n+1)=4.269371e+000; foe(n+1)=6.645929e+001; krok(n+1)=4.751121e-003; ng(n+1)=2.168542e+002;
n=82; farx(n+1)=4.348954e+000; foe(n+1)=6.641470e+001; krok(n+1)=5.852382e-003; ng(n+1)=6.152949e+002;
n=83; farx(n+1)=4.225303e+000; foe(n+1)=6.639908e+001; krok(n+1)=1.703305e-002; ng(n+1)=9.418504e+002;
n=84; farx(n+1)=3.889223e+000; foe(n+1)=6.624700e+001; krok(n+1)=2.977641e-002; ng(n+1)=9.858655e+002;
n=85; farx(n+1)=3.771971e+000; foe(n+1)=6.622496e+001; krok(n+1)=2.047906e-002; ng(n+1)=9.518294e+002;
n=86; farx(n+1)=3.739831e+000; foe(n+1)=6.613387e+001; krok(n+1)=1.138873e-001; ng(n+1)=1.052641e+003;
n=87; farx(n+1)=3.657395e+000; foe(n+1)=6.590331e+001; krok(n+1)=8.008859e-001; ng(n+1)=7.303876e+002;
n=88; farx(n+1)=3.674327e+000; foe(n+1)=6.569309e+001; krok(n+1)=7.545111e-001; ng(n+1)=8.881978e+002;
n=89; farx(n+1)=3.611013e+000; foe(n+1)=6.548348e+001; krok(n+1)=3.589369e-001; ng(n+1)=2.175251e+002;
n=90; farx(n+1)=3.543309e+000; foe(n+1)=6.537186e+001; krok(n+1)=2.944474e-001; ng(n+1)=1.442244e+003;
n=91; farx(n+1)=3.400876e+000; foe(n+1)=6.521015e+001; krok(n+1)=2.515191e-001; ng(n+1)=1.118166e+003;
n=92; farx(n+1)=3.364835e+000; foe(n+1)=6.501691e+001; krok(n+1)=1.257512e+000; ng(n+1)=8.956900e+002;
n=93; farx(n+1)=3.332148e+000; foe(n+1)=6.485344e+001; krok(n+1)=2.974827e-001; ng(n+1)=1.009389e+003;
n=94; farx(n+1)=3.311663e+000; foe(n+1)=6.452569e+001; krok(n+1)=7.382907e-001; ng(n+1)=8.287815e+002;
n=95; farx(n+1)=3.303845e+000; foe(n+1)=6.429636e+001; krok(n+1)=1.528330e-001; ng(n+1)=1.580172e+003;
n=96; farx(n+1)=3.275656e+000; foe(n+1)=6.414593e+001; krok(n+1)=1.530597e-001; ng(n+1)=1.535484e+003;
n=97; farx(n+1)=3.237438e+000; foe(n+1)=6.383836e+001; krok(n+1)=3.762824e-001; ng(n+1)=1.257354e+003;
n=98; farx(n+1)=3.160976e+000; foe(n+1)=6.353977e+001; krok(n+1)=2.730900e-001; ng(n+1)=1.750125e+003;
n=99; farx(n+1)=3.175109e+000; foe(n+1)=6.330941e+001; krok(n+1)=3.214125e-001; ng(n+1)=6.842168e+002;
n=100; farx(n+1)=3.209343e+000; foe(n+1)=6.304068e+001; krok(n+1)=2.226006e-001; ng(n+1)=1.707099e+003;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
