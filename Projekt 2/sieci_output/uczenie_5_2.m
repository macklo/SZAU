%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.009182e+003; foe(n+1)=4.014593e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.928010e+003; foe(n+1)=2.948513e+003; krok(n+1)=5.033545e-004; ng(n+1)=2.909504e+003;
n=2; farx(n+1)=1.077350e+003; foe(n+1)=1.141451e+003; krok(n+1)=3.223493e-003; ng(n+1)=1.629853e+003;
n=3; farx(n+1)=1.174452e+003; foe(n+1)=9.770200e+002; krok(n+1)=6.503398e-005; ng(n+1)=7.886044e+003;
n=4; farx(n+1)=2.468950e+003; foe(n+1)=6.528351e+002; krok(n+1)=1.442683e-003; ng(n+1)=9.545292e+003;
n=5; farx(n+1)=1.611356e+003; foe(n+1)=5.491162e+002; krok(n+1)=5.829556e-003; ng(n+1)=1.312746e+003;
n=6; farx(n+1)=1.188920e+003; foe(n+1)=4.938850e+002; krok(n+1)=1.244258e-004; ng(n+1)=2.508106e+003;
n=7; farx(n+1)=1.050149e+003; foe(n+1)=4.810470e+002; krok(n+1)=6.564208e-005; ng(n+1)=2.819320e+003;
n=8; farx(n+1)=6.460402e+002; foe(n+1)=4.166411e+002; krok(n+1)=1.284763e-003; ng(n+1)=7.815823e+003;
n=9; farx(n+1)=4.547228e+002; foe(n+1)=3.676745e+002; krok(n+1)=2.129131e-003; ng(n+1)=3.388456e+003;
n=10; farx(n+1)=2.435363e+002; foe(n+1)=3.155787e+002; krok(n+1)=2.118327e-003; ng(n+1)=2.128071e+003;
n=11; farx(n+1)=1.393684e+002; foe(n+1)=2.347788e+002; krok(n+1)=2.521031e-003; ng(n+1)=1.528136e+003;
n=12; farx(n+1)=1.536136e+002; foe(n+1)=2.247237e+002; krok(n+1)=1.066152e-003; ng(n+1)=7.714993e+002;
n=13; farx(n+1)=1.729477e+002; foe(n+1)=2.181458e+002; krok(n+1)=1.358214e-004; ng(n+1)=2.800346e+003;
n=14; farx(n+1)=2.073265e+002; foe(n+1)=1.822184e+002; krok(n+1)=1.954804e-003; ng(n+1)=6.294689e+002;
n=15; farx(n+1)=2.168091e+002; foe(n+1)=1.524909e+002; krok(n+1)=5.368339e-003; ng(n+1)=2.344550e+003;
n=16; farx(n+1)=2.369331e+002; foe(n+1)=1.426138e+002; krok(n+1)=5.981703e-003; ng(n+1)=2.118350e+003;
n=17; farx(n+1)=2.425726e+002; foe(n+1)=1.391913e+002; krok(n+1)=7.740461e-004; ng(n+1)=1.877528e+003;
n=18; farx(n+1)=1.723495e+002; foe(n+1)=1.234274e+002; krok(n+1)=9.299198e-003; ng(n+1)=1.706140e+003;
n=19; farx(n+1)=6.409204e+001; foe(n+1)=9.777789e+001; krok(n+1)=1.056509e-002; ng(n+1)=1.619525e+003;
n=20; farx(n+1)=4.866845e+001; foe(n+1)=8.639094e+001; krok(n+1)=2.013864e-003; ng(n+1)=1.321906e+003;
n=21; farx(n+1)=5.070611e+001; foe(n+1)=8.257903e+001; krok(n+1)=7.692316e-003; ng(n+1)=1.327321e+003;
n=22; farx(n+1)=4.676449e+001; foe(n+1)=7.701182e+001; krok(n+1)=1.941426e-002; ng(n+1)=1.443225e+003;
n=23; farx(n+1)=4.050378e+001; foe(n+1)=6.618534e+001; krok(n+1)=6.272594e-002; ng(n+1)=1.201846e+003;
n=24; farx(n+1)=3.549803e+001; foe(n+1)=5.754732e+001; krok(n+1)=3.535818e-002; ng(n+1)=1.158093e+003;
n=25; farx(n+1)=2.740054e+001; foe(n+1)=4.549376e+001; krok(n+1)=6.384007e-002; ng(n+1)=5.612637e+002;
%odnowa zmiennej metryki
n=26; farx(n+1)=2.714596e+001; foe(n+1)=4.493318e+001; krok(n+1)=2.299949e-005; ng(n+1)=4.792591e+002;
n=27; farx(n+1)=2.700232e+001; foe(n+1)=4.172320e+001; krok(n+1)=2.178332e-004; ng(n+1)=4.435629e+002;
n=28; farx(n+1)=2.734021e+001; foe(n+1)=3.883444e+001; krok(n+1)=7.968363e-004; ng(n+1)=2.235125e+002;
n=29; farx(n+1)=2.564410e+001; foe(n+1)=3.643923e+001; krok(n+1)=5.608389e-003; ng(n+1)=5.858692e+002;
n=30; farx(n+1)=2.302841e+001; foe(n+1)=3.512828e+001; krok(n+1)=3.063157e-004; ng(n+1)=5.759612e+002;
n=31; farx(n+1)=1.887338e+001; foe(n+1)=3.196766e+001; krok(n+1)=7.032104e-003; ng(n+1)=2.349903e+002;
n=32; farx(n+1)=1.744235e+001; foe(n+1)=3.092902e+001; krok(n+1)=2.144917e-003; ng(n+1)=6.501386e+002;
n=33; farx(n+1)=1.670834e+001; foe(n+1)=3.010197e+001; krok(n+1)=4.944619e-003; ng(n+1)=1.530510e+002;
n=34; farx(n+1)=1.557597e+001; foe(n+1)=2.851681e+001; krok(n+1)=2.021294e-003; ng(n+1)=4.091304e+002;
n=35; farx(n+1)=1.489429e+001; foe(n+1)=2.760556e+001; krok(n+1)=6.319511e-003; ng(n+1)=4.928749e+002;
n=36; farx(n+1)=1.348440e+001; foe(n+1)=2.608640e+001; krok(n+1)=1.716139e-002; ng(n+1)=4.632815e+002;
n=37; farx(n+1)=1.319301e+001; foe(n+1)=2.558354e+001; krok(n+1)=1.095618e-002; ng(n+1)=8.095227e+002;
n=38; farx(n+1)=1.324552e+001; foe(n+1)=2.504530e+001; krok(n+1)=1.640492e-002; ng(n+1)=6.111086e+002;
n=39; farx(n+1)=1.303897e+001; foe(n+1)=2.395871e+001; krok(n+1)=2.432476e-002; ng(n+1)=1.046652e+003;
n=40; farx(n+1)=1.413421e+001; foe(n+1)=2.200517e+001; krok(n+1)=3.497807e-002; ng(n+1)=7.243932e+002;
n=41; farx(n+1)=1.572448e+001; foe(n+1)=2.053909e+001; krok(n+1)=1.820180e-002; ng(n+1)=4.354615e+002;
n=42; farx(n+1)=1.534813e+001; foe(n+1)=1.612551e+001; krok(n+1)=3.023171e-001; ng(n+1)=1.069133e+003;
n=43; farx(n+1)=1.526417e+001; foe(n+1)=1.537645e+001; krok(n+1)=9.899545e-003; ng(n+1)=1.188590e+003;
n=44; farx(n+1)=1.337230e+001; foe(n+1)=1.309099e+001; krok(n+1)=1.088405e-001; ng(n+1)=3.189521e+002;
n=45; farx(n+1)=1.324680e+001; foe(n+1)=1.195670e+001; krok(n+1)=5.609462e-002; ng(n+1)=3.181686e+002;
n=46; farx(n+1)=1.116988e+001; foe(n+1)=1.055593e+001; krok(n+1)=5.312625e-001; ng(n+1)=1.702141e+002;
n=47; farx(n+1)=1.083300e+001; foe(n+1)=1.027462e+001; krok(n+1)=4.834060e-002; ng(n+1)=1.770937e+002;
n=48; farx(n+1)=9.496315e+000; foe(n+1)=9.531904e+000; krok(n+1)=2.320583e-001; ng(n+1)=2.766857e+002;
n=49; farx(n+1)=7.955224e+000; foe(n+1)=8.992283e+000; krok(n+1)=4.423764e-001; ng(n+1)=1.140513e+002;
n=50; farx(n+1)=6.406338e+000; foe(n+1)=8.455457e+000; krok(n+1)=6.211061e-001; ng(n+1)=2.953890e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=6.373425e+000; foe(n+1)=8.277729e+000; krok(n+1)=1.792006e-005; ng(n+1)=3.745240e+002;
n=52; farx(n+1)=6.373342e+000; foe(n+1)=8.222155e+000; krok(n+1)=7.471755e-006; ng(n+1)=3.115581e+002;
n=53; farx(n+1)=6.402427e+000; foe(n+1)=8.186743e+000; krok(n+1)=2.263251e-004; ng(n+1)=4.569334e+001;
n=54; farx(n+1)=6.363962e+000; foe(n+1)=8.167908e+000; krok(n+1)=4.788092e-004; ng(n+1)=2.900825e+001;
n=55; farx(n+1)=6.449181e+000; foe(n+1)=8.043774e+000; krok(n+1)=1.214548e-002; ng(n+1)=2.824082e+001;
n=56; farx(n+1)=6.311026e+000; foe(n+1)=7.927376e+000; krok(n+1)=2.511123e-003; ng(n+1)=5.479387e+001;
n=57; farx(n+1)=6.374554e+000; foe(n+1)=7.850341e+000; krok(n+1)=8.349282e-003; ng(n+1)=1.251466e+002;
n=58; farx(n+1)=5.835953e+000; foe(n+1)=7.673062e+000; krok(n+1)=1.958773e-002; ng(n+1)=1.179661e+002;
n=59; farx(n+1)=4.676532e+000; foe(n+1)=7.123251e+000; krok(n+1)=6.741320e-002; ng(n+1)=1.338199e+002;
n=60; farx(n+1)=4.132851e+000; foe(n+1)=6.735911e+000; krok(n+1)=1.401574e-002; ng(n+1)=6.358146e+002;
n=61; farx(n+1)=3.701131e+000; foe(n+1)=6.603814e+000; krok(n+1)=1.000091e-002; ng(n+1)=4.865401e+002;
n=62; farx(n+1)=3.046735e+000; foe(n+1)=6.294255e+000; krok(n+1)=1.555551e-002; ng(n+1)=7.103203e+002;
n=63; farx(n+1)=2.795655e+000; foe(n+1)=6.115667e+000; krok(n+1)=3.043820e-002; ng(n+1)=6.973604e+002;
n=64; farx(n+1)=2.647755e+000; foe(n+1)=5.868368e+000; krok(n+1)=4.382471e-002; ng(n+1)=7.133301e+002;
n=65; farx(n+1)=2.483527e+000; foe(n+1)=5.654456e+000; krok(n+1)=8.875371e-003; ng(n+1)=3.577056e+002;
n=66; farx(n+1)=2.419440e+000; foe(n+1)=5.360583e+000; krok(n+1)=5.442027e-002; ng(n+1)=3.753166e+002;
n=67; farx(n+1)=2.351413e+000; foe(n+1)=5.033619e+000; krok(n+1)=2.191235e-002; ng(n+1)=3.812524e+002;
n=68; farx(n+1)=2.023887e+000; foe(n+1)=4.514604e+000; krok(n+1)=1.365319e-001; ng(n+1)=2.359703e+002;
n=69; farx(n+1)=1.897610e+000; foe(n+1)=4.257570e+000; krok(n+1)=3.571281e-002; ng(n+1)=5.034755e+002;
n=70; farx(n+1)=1.701762e+000; foe(n+1)=3.556639e+000; krok(n+1)=3.138666e-001; ng(n+1)=5.552539e+002;
n=71; farx(n+1)=1.744779e+000; foe(n+1)=3.320679e+000; krok(n+1)=1.482616e-001; ng(n+1)=1.974040e+002;
n=72; farx(n+1)=1.764384e+000; foe(n+1)=3.118201e+000; krok(n+1)=1.134226e-001; ng(n+1)=3.050327e+002;
n=73; farx(n+1)=1.632461e+000; foe(n+1)=2.888080e+000; krok(n+1)=3.319551e-001; ng(n+1)=2.226281e+002;
n=74; farx(n+1)=1.577748e+000; foe(n+1)=2.799920e+000; krok(n+1)=2.538306e-001; ng(n+1)=1.598522e+002;
n=75; farx(n+1)=1.519823e+000; foe(n+1)=2.690139e+000; krok(n+1)=1.794685e-001; ng(n+1)=1.424827e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.519057e+000; foe(n+1)=2.676202e+000; krok(n+1)=2.311050e-006; ng(n+1)=2.623860e+002;
n=77; farx(n+1)=1.515074e+000; foe(n+1)=2.657966e+000; krok(n+1)=5.658127e-005; ng(n+1)=6.630782e+001;
n=78; farx(n+1)=1.507851e+000; foe(n+1)=2.643391e+000; krok(n+1)=3.572825e-005; ng(n+1)=7.499210e+001;
n=79; farx(n+1)=1.495575e+000; foe(n+1)=2.623924e+000; krok(n+1)=2.157932e-004; ng(n+1)=3.794732e+001;
n=80; farx(n+1)=1.492718e+000; foe(n+1)=2.607881e+000; krok(n+1)=4.053221e-004; ng(n+1)=2.286929e+001;
n=81; farx(n+1)=1.474679e+000; foe(n+1)=2.568262e+000; krok(n+1)=9.832841e-003; ng(n+1)=2.192823e+001;
n=82; farx(n+1)=1.459910e+000; foe(n+1)=2.529796e+000; krok(n+1)=5.608324e-003; ng(n+1)=1.773509e+001;
n=83; farx(n+1)=1.361001e+000; foe(n+1)=2.447265e+000; krok(n+1)=1.276026e-002; ng(n+1)=5.520031e+001;
n=84; farx(n+1)=1.335613e+000; foe(n+1)=2.412740e+000; krok(n+1)=2.387143e-002; ng(n+1)=1.107803e+002;
n=85; farx(n+1)=1.327294e+000; foe(n+1)=2.392464e+000; krok(n+1)=1.556168e-002; ng(n+1)=1.527568e+002;
n=86; farx(n+1)=1.319698e+000; foe(n+1)=2.329401e+000; krok(n+1)=2.419175e-002; ng(n+1)=2.591309e+002;
n=87; farx(n+1)=1.332853e+000; foe(n+1)=2.239627e+000; krok(n+1)=3.394649e-002; ng(n+1)=3.138357e+002;
n=88; farx(n+1)=1.346389e+000; foe(n+1)=2.168961e+000; krok(n+1)=8.430727e-002; ng(n+1)=1.342289e+002;
n=89; farx(n+1)=1.302126e+000; foe(n+1)=2.092711e+000; krok(n+1)=6.306675e-002; ng(n+1)=7.902945e+001;
n=90; farx(n+1)=1.299545e+000; foe(n+1)=2.026718e+000; krok(n+1)=2.223039e-002; ng(n+1)=3.741884e+002;
n=91; farx(n+1)=1.283764e+000; foe(n+1)=1.978611e+000; krok(n+1)=3.616714e-002; ng(n+1)=1.585774e+002;
n=92; farx(n+1)=1.259324e+000; foe(n+1)=1.874145e+000; krok(n+1)=4.773765e-002; ng(n+1)=2.999679e+002;
n=93; farx(n+1)=1.246204e+000; foe(n+1)=1.819170e+000; krok(n+1)=8.139610e-002; ng(n+1)=3.219950e+002;
n=94; farx(n+1)=1.245816e+000; foe(n+1)=1.758675e+000; krok(n+1)=8.892157e-002; ng(n+1)=6.061154e+002;
n=95; farx(n+1)=1.250250e+000; foe(n+1)=1.692506e+000; krok(n+1)=1.845727e-001; ng(n+1)=2.051945e+002;
n=96; farx(n+1)=1.216957e+000; foe(n+1)=1.623004e+000; krok(n+1)=2.223238e-001; ng(n+1)=2.085341e+002;
n=97; farx(n+1)=1.236801e+000; foe(n+1)=1.559541e+000; krok(n+1)=2.522670e-001; ng(n+1)=1.465137e+002;
n=98; farx(n+1)=1.224502e+000; foe(n+1)=1.511733e+000; krok(n+1)=1.328156e-001; ng(n+1)=1.064338e+002;
n=99; farx(n+1)=1.185315e+000; foe(n+1)=1.441500e+000; krok(n+1)=3.550835e-001; ng(n+1)=1.638796e+002;
n=100; farx(n+1)=1.138079e+000; foe(n+1)=1.362762e+000; krok(n+1)=1.284644e-001; ng(n+1)=2.152463e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)