%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.157627e+003; foe(n+1)=4.060044e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.994647e+003; foe(n+1)=2.935241e+003; krok(n+1)=5.236811e-004; ng(n+1)=3.988625e+003;
n=2; farx(n+1)=1.042492e+003; foe(n+1)=9.118507e+002; krok(n+1)=3.271540e-003; ng(n+1)=1.898583e+003;
n=3; farx(n+1)=1.131943e+003; foe(n+1)=7.056589e+002; krok(n+1)=1.306299e-004; ng(n+1)=6.170844e+003;
n=4; farx(n+1)=1.792592e+003; foe(n+1)=5.784896e+002; krok(n+1)=7.968363e-004; ng(n+1)=6.926498e+003;
n=5; farx(n+1)=1.463173e+003; foe(n+1)=5.148574e+002; krok(n+1)=6.391423e-003; ng(n+1)=1.028294e+003;
n=6; farx(n+1)=1.212083e+003; foe(n+1)=4.746528e+002; krok(n+1)=1.786211e-003; ng(n+1)=1.558818e+003;
n=7; farx(n+1)=1.072867e+003; foe(n+1)=4.629611e+002; krok(n+1)=3.968650e-004; ng(n+1)=1.728727e+003;
n=8; farx(n+1)=5.938549e+002; foe(n+1)=3.804118e+002; krok(n+1)=1.160015e-003; ng(n+1)=5.111471e+003;
n=9; farx(n+1)=3.597254e+002; foe(n+1)=3.341717e+002; krok(n+1)=1.626565e-003; ng(n+1)=2.742237e+003;
n=10; farx(n+1)=3.187672e+002; foe(n+1)=3.239340e+002; krok(n+1)=1.846469e-003; ng(n+1)=1.815067e+003;
n=11; farx(n+1)=3.121266e+002; foe(n+1)=3.222610e+002; krok(n+1)=1.024167e-003; ng(n+1)=7.944536e+002;
n=12; farx(n+1)=2.931579e+002; foe(n+1)=3.111679e+002; krok(n+1)=7.242402e-003; ng(n+1)=4.429434e+002;
n=13; farx(n+1)=2.836482e+002; foe(n+1)=3.006488e+002; krok(n+1)=1.006932e-003; ng(n+1)=1.254582e+003;
n=14; farx(n+1)=2.841278e+002; foe(n+1)=2.860707e+002; krok(n+1)=1.957529e-003; ng(n+1)=3.041772e+003;
n=15; farx(n+1)=2.990113e+002; foe(n+1)=2.802395e+002; krok(n+1)=1.123366e-003; ng(n+1)=1.964126e+003;
n=16; farx(n+1)=2.990814e+002; foe(n+1)=2.802230e+002; krok(n+1)=1.765316e-007; ng(n+1)=4.909034e+003;
n=17; farx(n+1)=2.981512e+002; foe(n+1)=2.800252e+002; krok(n+1)=1.315162e-004; ng(n+1)=2.829372e+003;
n=18; farx(n+1)=2.960402e+002; foe(n+1)=2.790010e+002; krok(n+1)=5.391128e-004; ng(n+1)=2.333339e+003;
n=19; farx(n+1)=3.138457e+002; foe(n+1)=2.722914e+002; krok(n+1)=6.130766e-003; ng(n+1)=5.674979e+003;
n=20; farx(n+1)=2.978521e+002; foe(n+1)=2.669809e+002; krok(n+1)=1.570350e-003; ng(n+1)=2.247470e+003;
n=21; farx(n+1)=2.031814e+002; foe(n+1)=2.255573e+002; krok(n+1)=1.853599e-002; ng(n+1)=2.714668e+003;
n=22; farx(n+1)=1.581896e+002; foe(n+1)=2.175565e+002; krok(n+1)=1.155113e-003; ng(n+1)=1.085629e+003;
n=23; farx(n+1)=8.275838e+001; foe(n+1)=1.971595e+002; krok(n+1)=6.099924e-003; ng(n+1)=1.118776e+003;
n=24; farx(n+1)=7.612454e+001; foe(n+1)=1.953568e+002; krok(n+1)=5.923015e-004; ng(n+1)=1.005128e+003;
n=25; farx(n+1)=8.102569e+001; foe(n+1)=1.929210e+002; krok(n+1)=5.548179e-003; ng(n+1)=1.618620e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=8.311854e+001; foe(n+1)=1.860805e+002; krok(n+1)=5.561724e-005; ng(n+1)=1.362326e+003;
n=27; farx(n+1)=5.682611e+001; foe(n+1)=1.691990e+002; krok(n+1)=4.838020e-004; ng(n+1)=7.670271e+002;
n=28; farx(n+1)=5.085791e+001; foe(n+1)=1.030430e+002; krok(n+1)=5.604282e-004; ng(n+1)=1.152106e+003;
n=29; farx(n+1)=5.187900e+001; foe(n+1)=1.007222e+002; krok(n+1)=2.860122e-005; ng(n+1)=3.365760e+003;
n=30; farx(n+1)=5.416252e+001; foe(n+1)=8.822642e+001; krok(n+1)=1.906580e-003; ng(n+1)=3.426198e+003;
n=31; farx(n+1)=6.183801e+001; foe(n+1)=8.371378e+001; krok(n+1)=1.401071e-004; ng(n+1)=3.876127e+003;
n=32; farx(n+1)=6.237279e+001; foe(n+1)=6.365864e+001; krok(n+1)=8.242481e-004; ng(n+1)=2.980155e+003;
n=33; farx(n+1)=6.250120e+001; foe(n+1)=5.892368e+001; krok(n+1)=1.310706e-003; ng(n+1)=2.011363e+003;
n=34; farx(n+1)=4.073690e+001; foe(n+1)=5.182240e+001; krok(n+1)=2.383354e-003; ng(n+1)=8.838252e+002;
n=35; farx(n+1)=3.465735e+001; foe(n+1)=4.909325e+001; krok(n+1)=2.855545e-003; ng(n+1)=1.107336e+003;
n=36; farx(n+1)=3.513491e+001; foe(n+1)=3.875453e+001; krok(n+1)=5.296999e-003; ng(n+1)=5.115601e+002;
n=37; farx(n+1)=3.439782e+001; foe(n+1)=3.839282e+001; krok(n+1)=6.518159e-004; ng(n+1)=1.683338e+003;
n=38; farx(n+1)=3.153559e+001; foe(n+1)=3.569349e+001; krok(n+1)=3.376366e-004; ng(n+1)=8.020111e+003;
n=39; farx(n+1)=3.190498e+001; foe(n+1)=3.312555e+001; krok(n+1)=1.006709e-003; ng(n+1)=5.172093e+003;
n=40; farx(n+1)=2.827240e+001; foe(n+1)=3.174363e+001; krok(n+1)=5.548179e-003; ng(n+1)=8.002427e+002;
n=41; farx(n+1)=2.444304e+001; foe(n+1)=2.945846e+001; krok(n+1)=2.855545e-003; ng(n+1)=3.861636e+003;
n=42; farx(n+1)=2.273366e+001; foe(n+1)=2.789759e+001; krok(n+1)=9.416683e-003; ng(n+1)=8.091657e+002;
n=43; farx(n+1)=2.078439e+001; foe(n+1)=2.674259e+001; krok(n+1)=3.272561e-003; ng(n+1)=9.792182e+002;
n=44; farx(n+1)=1.587386e+001; foe(n+1)=2.323358e+001; krok(n+1)=2.855545e-003; ng(n+1)=1.013727e+003;
n=45; farx(n+1)=1.517472e+001; foe(n+1)=2.260740e+001; krok(n+1)=4.884963e-003; ng(n+1)=3.737117e+002;
n=46; farx(n+1)=1.488493e+001; foe(n+1)=2.164157e+001; krok(n+1)=7.557177e-003; ng(n+1)=7.571020e+002;
n=47; farx(n+1)=1.427640e+001; foe(n+1)=2.035494e+001; krok(n+1)=2.529683e-002; ng(n+1)=5.934607e+002;
n=48; farx(n+1)=1.307109e+001; foe(n+1)=1.944460e+001; krok(n+1)=1.477175e-002; ng(n+1)=4.258365e+002;
n=49; farx(n+1)=1.178849e+001; foe(n+1)=1.880347e+001; krok(n+1)=3.432277e-002; ng(n+1)=3.711808e+002;
n=50; farx(n+1)=1.109988e+001; foe(n+1)=1.791510e+001; krok(n+1)=2.001452e-002; ng(n+1)=5.010336e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=1.094949e+001; foe(n+1)=1.771362e+001; krok(n+1)=1.722049e-005; ng(n+1)=4.422728e+002;
n=52; farx(n+1)=1.079906e+001; foe(n+1)=1.738627e+001; krok(n+1)=9.848602e-005; ng(n+1)=2.855218e+002;
n=53; farx(n+1)=1.074634e+001; foe(n+1)=1.695343e+001; krok(n+1)=8.503167e-004; ng(n+1)=1.074568e+002;
n=54; farx(n+1)=9.909406e+000; foe(n+1)=1.498193e+001; krok(n+1)=8.523876e-004; ng(n+1)=2.263186e+002;
n=55; farx(n+1)=8.962632e+000; foe(n+1)=1.451218e+001; krok(n+1)=1.982946e-003; ng(n+1)=4.749529e+002;
n=56; farx(n+1)=8.767521e+000; foe(n+1)=1.396337e+001; krok(n+1)=3.621201e-003; ng(n+1)=4.440752e+002;
n=57; farx(n+1)=8.629646e+000; foe(n+1)=1.348463e+001; krok(n+1)=2.739044e-003; ng(n+1)=1.321566e+002;
n=58; farx(n+1)=8.065016e+000; foe(n+1)=1.205608e+001; krok(n+1)=1.865271e-002; ng(n+1)=3.224690e+002;
n=59; farx(n+1)=7.230952e+000; foe(n+1)=1.103460e+001; krok(n+1)=3.625910e-003; ng(n+1)=3.435377e+002;
n=60; farx(n+1)=6.934462e+000; foe(n+1)=1.071038e+001; krok(n+1)=1.100811e-003; ng(n+1)=4.980708e+002;
n=61; farx(n+1)=6.140962e+000; foe(n+1)=9.973539e+000; krok(n+1)=3.995494e-003; ng(n+1)=2.075060e+002;
n=62; farx(n+1)=5.706386e+000; foe(n+1)=9.374801e+000; krok(n+1)=1.382426e-002; ng(n+1)=5.848504e+002;
n=63; farx(n+1)=5.537022e+000; foe(n+1)=9.201619e+000; krok(n+1)=4.440041e-003; ng(n+1)=7.426214e+002;
n=64; farx(n+1)=5.034518e+000; foe(n+1)=8.856669e+000; krok(n+1)=8.704305e-003; ng(n+1)=7.260771e+002;
n=65; farx(n+1)=4.669720e+000; foe(n+1)=8.522188e+000; krok(n+1)=1.680696e-002; ng(n+1)=8.176551e+002;
n=66; farx(n+1)=4.258380e+000; foe(n+1)=8.068895e+000; krok(n+1)=9.784074e-003; ng(n+1)=4.218947e+002;
n=67; farx(n+1)=3.993502e+000; foe(n+1)=7.780292e+000; krok(n+1)=2.709456e-002; ng(n+1)=5.493608e+002;
n=68; farx(n+1)=3.864290e+000; foe(n+1)=7.671717e+000; krok(n+1)=1.417782e-002; ng(n+1)=4.974585e+002;
n=69; farx(n+1)=3.541605e+000; foe(n+1)=7.049254e+000; krok(n+1)=7.071636e-002; ng(n+1)=6.756311e+002;
n=70; farx(n+1)=3.404056e+000; foe(n+1)=6.735246e+000; krok(n+1)=1.729712e-002; ng(n+1)=2.841581e+002;
n=71; farx(n+1)=3.277199e+000; foe(n+1)=6.524845e+000; krok(n+1)=3.818352e-002; ng(n+1)=4.967763e+002;
n=72; farx(n+1)=3.063440e+000; foe(n+1)=6.162101e+000; krok(n+1)=4.382471e-002; ng(n+1)=2.939184e+002;
n=73; farx(n+1)=3.066289e+000; foe(n+1)=6.044220e+000; krok(n+1)=2.008299e-002; ng(n+1)=1.975792e+002;
n=74; farx(n+1)=2.930596e+000; foe(n+1)=5.729911e+000; krok(n+1)=6.288276e-002; ng(n+1)=3.057891e+002;
n=75; farx(n+1)=2.862114e+000; foe(n+1)=5.428893e+000; krok(n+1)=7.000904e-002; ng(n+1)=1.491257e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=2.838914e+000; foe(n+1)=5.346205e+000; krok(n+1)=9.054721e-006; ng(n+1)=4.248376e+002;
n=77; farx(n+1)=2.822304e+000; foe(n+1)=5.318816e+000; krok(n+1)=6.836820e-005; ng(n+1)=9.327068e+001;
n=78; farx(n+1)=2.809964e+000; foe(n+1)=5.302417e+000; krok(n+1)=8.370634e-005; ng(n+1)=6.480279e+001;
n=79; farx(n+1)=2.743105e+000; foe(n+1)=5.040044e+000; krok(n+1)=1.758026e-003; ng(n+1)=5.754029e+001;
n=80; farx(n+1)=2.646437e+000; foe(n+1)=4.883509e+000; krok(n+1)=1.441974e-003; ng(n+1)=1.174464e+002;
n=81; farx(n+1)=2.553814e+000; foe(n+1)=4.760445e+000; krok(n+1)=9.659422e-003; ng(n+1)=2.901170e+002;
n=82; farx(n+1)=2.487565e+000; foe(n+1)=4.637157e+000; krok(n+1)=6.970444e-003; ng(n+1)=2.388040e+002;
n=83; farx(n+1)=2.435145e+000; foe(n+1)=4.591815e+000; krok(n+1)=2.511123e-003; ng(n+1)=2.437528e+002;
n=84; farx(n+1)=2.228282e+000; foe(n+1)=4.373410e+000; krok(n+1)=2.967304e-002; ng(n+1)=2.730430e+002;
n=85; farx(n+1)=2.088016e+000; foe(n+1)=4.209837e+000; krok(n+1)=8.085176e-003; ng(n+1)=4.120444e+002;
n=86; farx(n+1)=2.023122e+000; foe(n+1)=4.111531e+000; krok(n+1)=4.191618e-003; ng(n+1)=5.125611e+002;
n=87; farx(n+1)=1.983223e+000; foe(n+1)=4.050428e+000; krok(n+1)=2.038894e-002; ng(n+1)=2.670697e+002;
n=88; farx(n+1)=1.913154e+000; foe(n+1)=3.892524e+000; krok(n+1)=3.640360e-002; ng(n+1)=4.729746e+002;
n=89; farx(n+1)=1.887451e+000; foe(n+1)=3.819888e+000; krok(n+1)=7.366454e-003; ng(n+1)=1.058764e+002;
n=90; farx(n+1)=1.847819e+000; foe(n+1)=3.700278e+000; krok(n+1)=2.243330e-002; ng(n+1)=3.037918e+002;
n=91; farx(n+1)=1.816207e+000; foe(n+1)=3.582984e+000; krok(n+1)=3.618094e-002; ng(n+1)=2.132313e+002;
n=92; farx(n+1)=1.754260e+000; foe(n+1)=3.465136e+000; krok(n+1)=5.335388e-002; ng(n+1)=1.979738e+002;
n=93; farx(n+1)=1.746456e+000; foe(n+1)=3.446205e+000; krok(n+1)=2.943935e-003; ng(n+1)=2.480572e+002;
n=94; farx(n+1)=1.676204e+000; foe(n+1)=3.342146e+000; krok(n+1)=4.895039e-002; ng(n+1)=1.749516e+002;
n=95; farx(n+1)=1.700739e+000; foe(n+1)=3.285507e+000; krok(n+1)=2.642647e-002; ng(n+1)=1.889246e+002;
n=96; farx(n+1)=1.855059e+000; foe(n+1)=3.091548e+000; krok(n+1)=2.785422e-002; ng(n+1)=2.935380e+002;
n=97; farx(n+1)=1.809531e+000; foe(n+1)=3.048018e+000; krok(n+1)=4.544795e-002; ng(n+1)=1.905807e+002;
n=98; farx(n+1)=1.724569e+000; foe(n+1)=2.945499e+000; krok(n+1)=1.827549e-001; ng(n+1)=9.534508e+001;
n=99; farx(n+1)=1.618057e+000; foe(n+1)=2.881056e+000; krok(n+1)=1.008497e-001; ng(n+1)=8.297389e+001;
n=100; farx(n+1)=1.518192e+000; foe(n+1)=2.824125e+000; krok(n+1)=1.702433e-001; ng(n+1)=1.268759e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)