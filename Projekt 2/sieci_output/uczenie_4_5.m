%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.631458e+003; foe(n+1)=4.809025e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.386650e+003; foe(n+1)=3.631898e+003; krok(n+1)=3.467612e-004; ng(n+1)=3.832883e+003;
n=2; farx(n+1)=9.229092e+002; foe(n+1)=1.186364e+003; krok(n+1)=1.171923e-003; ng(n+1)=3.395992e+003;
n=3; farx(n+1)=1.055154e+003; foe(n+1)=8.876355e+002; krok(n+1)=1.366331e-004; ng(n+1)=7.255884e+003;
n=4; farx(n+1)=1.468622e+003; foe(n+1)=8.275979e+002; krok(n+1)=1.043660e-003; ng(n+1)=4.999914e+003;
n=5; farx(n+1)=1.044060e+003; foe(n+1)=6.749459e+002; krok(n+1)=6.854056e-003; ng(n+1)=1.337910e+003;
n=6; farx(n+1)=6.167256e+002; foe(n+1)=5.673887e+002; krok(n+1)=5.045766e-004; ng(n+1)=2.213380e+003;
n=7; farx(n+1)=3.703951e+002; foe(n+1)=4.975146e+002; krok(n+1)=4.129769e-004; ng(n+1)=4.404876e+003;
n=8; farx(n+1)=3.120800e+002; foe(n+1)=4.760347e+002; krok(n+1)=6.145526e-004; ng(n+1)=4.808395e+003;
n=9; farx(n+1)=2.468932e+002; foe(n+1)=4.492440e+002; krok(n+1)=4.872731e-004; ng(n+1)=4.565959e+003;
n=10; farx(n+1)=2.292904e+002; foe(n+1)=4.443877e+002; krok(n+1)=1.974847e-004; ng(n+1)=3.182182e+003;
n=11; farx(n+1)=2.180614e+002; foe(n+1)=3.937198e+002; krok(n+1)=3.127439e-003; ng(n+1)=2.129038e+003;
n=12; farx(n+1)=1.449095e+002; foe(n+1)=3.690669e+002; krok(n+1)=6.728198e-003; ng(n+1)=2.293483e+003;
n=13; farx(n+1)=9.344735e+001; foe(n+1)=3.505460e+002; krok(n+1)=1.303632e-003; ng(n+1)=1.829356e+003;
n=14; farx(n+1)=8.777296e+001; foe(n+1)=3.480338e+002; krok(n+1)=3.794876e-004; ng(n+1)=2.270290e+003;
n=15; farx(n+1)=5.920316e+001; foe(n+1)=3.191262e+002; krok(n+1)=5.318174e-003; ng(n+1)=2.946071e+003;
n=16; farx(n+1)=5.505780e+001; foe(n+1)=2.992005e+002; krok(n+1)=1.406421e-002; ng(n+1)=5.434308e+003;
n=17; farx(n+1)=5.014856e+001; foe(n+1)=2.721802e+002; krok(n+1)=3.174920e-003; ng(n+1)=3.107599e+003;
n=18; farx(n+1)=4.751372e+001; foe(n+1)=2.584241e+002; krok(n+1)=8.956895e-003; ng(n+1)=3.466115e+003;
n=19; farx(n+1)=4.865803e+001; foe(n+1)=2.476026e+002; krok(n+1)=5.548179e-003; ng(n+1)=3.067110e+003;
n=20; farx(n+1)=4.577829e+001; foe(n+1)=2.136907e+002; krok(n+1)=4.614317e-002; ng(n+1)=3.991324e+003;
n=21; farx(n+1)=4.312374e+001; foe(n+1)=2.033585e+002; krok(n+1)=3.276344e-002; ng(n+1)=3.870794e+003;
n=22; farx(n+1)=4.217841e+001; foe(n+1)=1.801135e+002; krok(n+1)=6.306675e-002; ng(n+1)=2.572436e+003;
n=23; farx(n+1)=4.041841e+001; foe(n+1)=1.571541e+002; krok(n+1)=9.326072e-002; ng(n+1)=3.240441e+003;
n=24; farx(n+1)=3.269816e+001; foe(n+1)=1.197933e+002; krok(n+1)=3.614383e-001; ng(n+1)=1.262858e+003;
n=25; farx(n+1)=3.282167e+001; foe(n+1)=1.123661e+002; krok(n+1)=1.775417e-001; ng(n+1)=1.045586e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=3.326412e+001; foe(n+1)=1.112268e+002; krok(n+1)=6.129763e-006; ng(n+1)=1.463975e+003;
n=27; farx(n+1)=3.262591e+001; foe(n+1)=1.101563e+002; krok(n+1)=1.891762e-004; ng(n+1)=2.649760e+002;
n=28; farx(n+1)=3.241646e+001; foe(n+1)=1.092616e+002; krok(n+1)=5.429951e-005; ng(n+1)=4.480233e+002;
n=29; farx(n+1)=3.136541e+001; foe(n+1)=1.074693e+002; krok(n+1)=3.505243e-004; ng(n+1)=4.814335e+002;
n=30; farx(n+1)=3.233426e+001; foe(n+1)=1.046952e+002; krok(n+1)=6.647718e-004; ng(n+1)=5.121948e+002;
n=31; farx(n+1)=3.016525e+001; foe(n+1)=9.477874e+001; krok(n+1)=7.251821e-003; ng(n+1)=1.625445e+003;
n=32; farx(n+1)=2.922280e+001; foe(n+1)=9.078506e+001; krok(n+1)=1.679702e-003; ng(n+1)=1.854926e+003;
n=33; farx(n+1)=2.938209e+001; foe(n+1)=8.711428e+001; krok(n+1)=1.387045e-003; ng(n+1)=2.309656e+003;
n=34; farx(n+1)=2.837702e+001; foe(n+1)=8.032256e+001; krok(n+1)=1.296700e-003; ng(n+1)=2.666441e+003;
n=35; farx(n+1)=2.826733e+001; foe(n+1)=7.575883e+001; krok(n+1)=2.849598e-003; ng(n+1)=3.232877e+003;
n=36; farx(n+1)=2.756061e+001; foe(n+1)=7.311082e+001; krok(n+1)=2.926191e-003; ng(n+1)=2.940795e+003;
n=37; farx(n+1)=2.427253e+001; foe(n+1)=6.012971e+001; krok(n+1)=3.188116e-002; ng(n+1)=3.229217e+003;
n=38; farx(n+1)=2.403610e+001; foe(n+1)=5.763924e+001; krok(n+1)=2.580064e-002; ng(n+1)=1.683906e+003;
n=39; farx(n+1)=2.069752e+001; foe(n+1)=4.856067e+001; krok(n+1)=2.097129e-002; ng(n+1)=7.644270e+002;
n=40; farx(n+1)=1.936036e+001; foe(n+1)=4.548782e+001; krok(n+1)=1.703305e-002; ng(n+1)=1.591529e+003;
n=41; farx(n+1)=1.243396e+001; foe(n+1)=3.313966e+001; krok(n+1)=9.946858e-002; ng(n+1)=2.314432e+003;
n=42; farx(n+1)=7.233983e+000; foe(n+1)=2.520078e+001; krok(n+1)=9.365437e-002; ng(n+1)=2.115542e+003;
n=43; farx(n+1)=4.999871e+000; foe(n+1)=2.072565e+001; krok(n+1)=1.344556e-001; ng(n+1)=1.245001e+003;
n=44; farx(n+1)=4.385811e+000; foe(n+1)=1.751650e+001; krok(n+1)=6.273346e-001; ng(n+1)=9.333135e+002;
n=45; farx(n+1)=4.574007e+000; foe(n+1)=1.663256e+001; krok(n+1)=2.740921e-001; ng(n+1)=2.927917e+002;
n=46; farx(n+1)=4.781265e+000; foe(n+1)=1.604749e+001; krok(n+1)=4.033988e-001; ng(n+1)=3.369468e+002;
n=47; farx(n+1)=5.065188e+000; foe(n+1)=1.522265e+001; krok(n+1)=4.102061e-001; ng(n+1)=1.608476e+002;
n=48; farx(n+1)=4.154611e+000; foe(n+1)=1.439850e+001; krok(n+1)=9.859540e-001; ng(n+1)=2.226348e+002;
n=49; farx(n+1)=3.837804e+000; foe(n+1)=1.388940e+001; krok(n+1)=6.065168e-001; ng(n+1)=2.986177e+002;
n=50; farx(n+1)=3.266552e+000; foe(n+1)=1.318085e+001; krok(n+1)=9.685429e-001; ng(n+1)=4.941240e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=3.243917e+000; foe(n+1)=1.306976e+001; krok(n+1)=5.874564e-006; ng(n+1)=4.534815e+002;
n=52; farx(n+1)=3.226540e+000; foe(n+1)=1.300113e+001; krok(n+1)=3.126325e-005; ng(n+1)=1.263682e+002;
n=53; farx(n+1)=3.215691e+000; foe(n+1)=1.295363e+001; krok(n+1)=4.525103e-005; ng(n+1)=1.092274e+002;
n=54; farx(n+1)=3.351777e+000; foe(n+1)=1.265047e+001; krok(n+1)=4.129769e-004; ng(n+1)=8.421441e+001;
n=55; farx(n+1)=3.341057e+000; foe(n+1)=1.259719e+001; krok(n+1)=4.890003e-004; ng(n+1)=5.301712e+001;
n=56; farx(n+1)=3.250873e+000; foe(n+1)=1.246040e+001; krok(n+1)=2.804162e-003; ng(n+1)=4.147714e+001;
n=57; farx(n+1)=3.167250e+000; foe(n+1)=1.225406e+001; krok(n+1)=9.175071e-003; ng(n+1)=9.568352e+001;
n=58; farx(n+1)=3.095548e+000; foe(n+1)=1.212470e+001; krok(n+1)=5.999763e-003; ng(n+1)=4.179416e+002;
n=59; farx(n+1)=3.018164e+000; foe(n+1)=1.201172e+001; krok(n+1)=6.669234e-003; ng(n+1)=6.461023e+002;
n=60; farx(n+1)=2.803971e+000; foe(n+1)=1.176481e+001; krok(n+1)=2.788178e-002; ng(n+1)=4.858528e+002;
n=61; farx(n+1)=2.566405e+000; foe(n+1)=1.144770e+001; krok(n+1)=8.877087e-002; ng(n+1)=4.689298e+002;
n=62; farx(n+1)=2.344412e+000; foe(n+1)=1.116111e+001; krok(n+1)=1.957537e-002; ng(n+1)=7.052953e+002;
n=63; farx(n+1)=2.150991e+000; foe(n+1)=1.093371e+001; krok(n+1)=6.690246e-002; ng(n+1)=6.584935e+002;
n=64; farx(n+1)=1.956747e+000; foe(n+1)=1.072477e+001; krok(n+1)=3.550149e-002; ng(n+1)=5.922554e+002;
n=65; farx(n+1)=1.729471e+000; foe(n+1)=1.036399e+001; krok(n+1)=7.197563e-002; ng(n+1)=4.273825e+002;
n=66; farx(n+1)=1.382185e+000; foe(n+1)=9.829224e+000; krok(n+1)=2.022243e-001; ng(n+1)=4.935214e+002;
n=67; farx(n+1)=1.221588e+000; foe(n+1)=9.530248e+000; krok(n+1)=4.111242e-002; ng(n+1)=3.155173e+002;
n=68; farx(n+1)=1.266233e+000; foe(n+1)=8.972725e+000; krok(n+1)=8.005809e-002; ng(n+1)=9.346912e+002;
n=69; farx(n+1)=1.373711e+000; foe(n+1)=8.416427e+000; krok(n+1)=2.031949e-001; ng(n+1)=9.395864e+002;
n=70; farx(n+1)=1.353997e+000; foe(n+1)=8.146701e+000; krok(n+1)=2.223238e-001; ng(n+1)=3.459587e+002;
n=71; farx(n+1)=1.280887e+000; foe(n+1)=7.883145e+000; krok(n+1)=1.399123e-001; ng(n+1)=5.862334e+002;
n=72; farx(n+1)=1.363707e+000; foe(n+1)=7.597825e+000; krok(n+1)=1.928435e-001; ng(n+1)=3.351838e+002;
n=73; farx(n+1)=1.196310e+000; foe(n+1)=7.298446e+000; krok(n+1)=1.410702e-001; ng(n+1)=3.863617e+002;
n=74; farx(n+1)=1.057646e+000; foe(n+1)=6.911437e+000; krok(n+1)=2.800362e-001; ng(n+1)=3.822279e+002;
n=75; farx(n+1)=1.030081e+000; foe(n+1)=6.315733e+000; krok(n+1)=2.093785e-001; ng(n+1)=3.374804e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.030569e+000; foe(n+1)=6.274077e+000; krok(n+1)=1.086355e-006; ng(n+1)=6.037098e+002;
n=77; farx(n+1)=1.031021e+000; foe(n+1)=6.240593e+000; krok(n+1)=2.417991e-006; ng(n+1)=4.072073e+002;
n=78; farx(n+1)=1.030567e+000; foe(n+1)=6.230465e+000; krok(n+1)=1.856604e-006; ng(n+1)=2.613531e+002;
n=79; farx(n+1)=1.027202e+000; foe(n+1)=6.104547e+000; krok(n+1)=1.367364e-004; ng(n+1)=1.094649e+002;
n=80; farx(n+1)=1.032267e+000; foe(n+1)=6.047755e+000; krok(n+1)=9.675577e-005; ng(n+1)=9.294669e+001;
n=81; farx(n+1)=1.036457e+000; foe(n+1)=5.974311e+000; krok(n+1)=1.355019e-004; ng(n+1)=6.901623e+001;
n=82; farx(n+1)=1.037051e+000; foe(n+1)=5.774054e+000; krok(n+1)=4.264608e-003; ng(n+1)=8.368150e+001;
n=83; farx(n+1)=1.033025e+000; foe(n+1)=5.736848e+000; krok(n+1)=2.275225e-003; ng(n+1)=7.960100e+002;
n=84; farx(n+1)=1.029116e+000; foe(n+1)=5.654782e+000; krok(n+1)=2.438216e-003; ng(n+1)=8.328935e+002;
n=85; farx(n+1)=1.069891e+000; foe(n+1)=5.590570e+000; krok(n+1)=1.522640e-002; ng(n+1)=9.151438e+002;
n=86; farx(n+1)=1.132116e+000; foe(n+1)=5.363936e+000; krok(n+1)=2.679765e-002; ng(n+1)=1.012078e+003;
n=87; farx(n+1)=1.139612e+000; foe(n+1)=5.319559e+000; krok(n+1)=3.966104e-003; ng(n+1)=2.261221e+002;
n=88; farx(n+1)=1.158120e+000; foe(n+1)=5.192938e+000; krok(n+1)=5.027940e-003; ng(n+1)=6.792890e+002;
n=89; farx(n+1)=1.161403e+000; foe(n+1)=5.155796e+000; krok(n+1)=4.640059e-003; ng(n+1)=8.015214e+002;
n=90; farx(n+1)=1.188521e+000; foe(n+1)=5.033335e+000; krok(n+1)=9.464413e-002; ng(n+1)=2.004830e+002;
n=91; farx(n+1)=1.129468e+000; foe(n+1)=4.882450e+000; krok(n+1)=1.170340e-001; ng(n+1)=3.996513e+002;
n=92; farx(n+1)=1.098024e+000; foe(n+1)=4.784409e+000; krok(n+1)=1.390811e-001; ng(n+1)=3.993290e+002;
n=93; farx(n+1)=1.041377e+000; foe(n+1)=4.678581e+000; krok(n+1)=6.345766e-002; ng(n+1)=1.347383e+003;
n=94; farx(n+1)=1.041691e+000; foe(n+1)=4.624715e+000; krok(n+1)=2.502130e-002; ng(n+1)=1.248146e+003;
n=95; farx(n+1)=1.103435e+000; foe(n+1)=4.533053e+000; krok(n+1)=1.341318e-001; ng(n+1)=1.257420e+003;
n=96; farx(n+1)=9.656585e-001; foe(n+1)=4.395146e+000; krok(n+1)=1.380129e-001; ng(n+1)=6.251066e+002;
n=97; farx(n+1)=9.444216e-001; foe(n+1)=4.174806e+000; krok(n+1)=5.378226e-001; ng(n+1)=6.112437e+002;
n=98; farx(n+1)=9.350498e-001; foe(n+1)=4.051760e+000; krok(n+1)=3.456460e-001; ng(n+1)=2.132482e+002;
n=99; farx(n+1)=9.404091e-001; foe(n+1)=3.921033e+000; krok(n+1)=3.054682e-001; ng(n+1)=1.332467e+003;
n=100; farx(n+1)=9.292755e-001; foe(n+1)=3.873767e+000; krok(n+1)=1.818254e-001; ng(n+1)=8.956009e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
