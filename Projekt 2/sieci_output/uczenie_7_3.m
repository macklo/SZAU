%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.900550e+003; foe(n+1)=4.844621e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.777687e+003; foe(n+1)=3.739398e+003; krok(n+1)=3.771986e-004; ng(n+1)=4.622693e+003;
n=2; farx(n+1)=9.015118e+002; foe(n+1)=8.593221e+002; krok(n+1)=1.772228e-003; ng(n+1)=3.243109e+003;
n=3; farx(n+1)=9.851008e+002; foe(n+1)=7.900242e+002; krok(n+1)=1.784716e-004; ng(n+1)=3.168171e+003;
n=4; farx(n+1)=6.361971e+002; foe(n+1)=7.363453e+002; krok(n+1)=1.636169e-003; ng(n+1)=8.860832e+002;
n=5; farx(n+1)=2.153841e+002; foe(n+1)=3.720750e+002; krok(n+1)=7.433745e-004; ng(n+1)=4.864494e+003;
n=6; farx(n+1)=2.077321e+002; foe(n+1)=3.476561e+002; krok(n+1)=4.644788e-004; ng(n+1)=1.869067e+003;
n=7; farx(n+1)=1.659300e+002; foe(n+1)=3.303307e+002; krok(n+1)=8.100612e-004; ng(n+1)=1.473137e+003;
n=8; farx(n+1)=1.309179e+002; foe(n+1)=3.010029e+002; krok(n+1)=1.256985e-003; ng(n+1)=8.858443e+002;
n=9; farx(n+1)=8.504755e+001; foe(n+1)=2.668950e+002; krok(n+1)=3.059399e-003; ng(n+1)=1.234847e+003;
n=10; farx(n+1)=8.084475e+001; foe(n+1)=2.626913e+002; krok(n+1)=1.075741e-004; ng(n+1)=1.608513e+003;
n=11; farx(n+1)=6.410702e+001; foe(n+1)=2.425711e+002; krok(n+1)=1.925691e-003; ng(n+1)=1.854815e+003;
n=12; farx(n+1)=6.085038e+001; foe(n+1)=2.393408e+002; krok(n+1)=1.169035e-004; ng(n+1)=3.509466e+003;
n=13; farx(n+1)=5.822865e+001; foe(n+1)=2.316267e+002; krok(n+1)=9.064776e-004; ng(n+1)=4.275410e+003;
n=14; farx(n+1)=5.263897e+001; foe(n+1)=2.258189e+002; krok(n+1)=1.420248e-003; ng(n+1)=5.552305e+003;
n=15; farx(n+1)=4.581702e+001; foe(n+1)=2.134603e+002; krok(n+1)=1.407045e-004; ng(n+1)=7.145513e+003;
n=16; farx(n+1)=3.891815e+001; foe(n+1)=1.806908e+002; krok(n+1)=3.026819e-003; ng(n+1)=9.589022e+003;
n=17; farx(n+1)=3.935506e+001; foe(n+1)=1.725593e+002; krok(n+1)=1.484824e-004; ng(n+1)=5.711193e+003;
n=18; farx(n+1)=3.684283e+001; foe(n+1)=1.399367e+002; krok(n+1)=1.317301e-003; ng(n+1)=2.810047e+003;
n=19; farx(n+1)=4.072292e+001; foe(n+1)=1.281515e+002; krok(n+1)=2.054008e-004; ng(n+1)=5.186495e+003;
n=20; farx(n+1)=3.461197e+001; foe(n+1)=1.153057e+002; krok(n+1)=9.053002e-004; ng(n+1)=1.901741e+003;
n=21; farx(n+1)=3.289414e+001; foe(n+1)=1.108421e+002; krok(n+1)=1.532691e-003; ng(n+1)=3.122295e+003;
n=22; farx(n+1)=3.467670e+001; foe(n+1)=1.043055e+002; krok(n+1)=5.188110e-004; ng(n+1)=2.623619e+003;
n=23; farx(n+1)=2.927335e+001; foe(n+1)=8.508105e+001; krok(n+1)=2.482161e-003; ng(n+1)=3.830536e+003;
n=24; farx(n+1)=2.992753e+001; foe(n+1)=7.882884e+001; krok(n+1)=4.748903e-004; ng(n+1)=1.613255e+003;
n=25; farx(n+1)=2.967636e+001; foe(n+1)=7.710364e+001; krok(n+1)=7.937301e-004; ng(n+1)=2.584678e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=2.942391e+001; foe(n+1)=7.405258e+001; krok(n+1)=2.769106e-005; ng(n+1)=1.484705e+003;
n=27; farx(n+1)=2.924116e+001; foe(n+1)=7.246851e+001; krok(n+1)=3.954660e-006; ng(n+1)=2.706737e+003;
n=28; farx(n+1)=2.821097e+001; foe(n+1)=7.043619e+001; krok(n+1)=3.014162e-005; ng(n+1)=1.250203e+003;
n=29; farx(n+1)=2.403848e+001; foe(n+1)=5.655155e+001; krok(n+1)=3.100233e-004; ng(n+1)=1.109037e+003;
n=30; farx(n+1)=2.098453e+001; foe(n+1)=3.879532e+001; krok(n+1)=1.908037e-003; ng(n+1)=5.749727e+003;
n=31; farx(n+1)=1.910005e+001; foe(n+1)=3.442225e+001; krok(n+1)=1.516248e-004; ng(n+1)=1.015814e+003;
n=32; farx(n+1)=1.832976e+001; foe(n+1)=3.058580e+001; krok(n+1)=7.968363e-004; ng(n+1)=2.379154e+003;
n=33; farx(n+1)=1.708312e+001; foe(n+1)=2.851828e+001; krok(n+1)=9.317073e-004; ng(n+1)=1.816092e+003;
n=34; farx(n+1)=1.486035e+001; foe(n+1)=2.491997e+001; krok(n+1)=1.667309e-003; ng(n+1)=1.537332e+003;
n=35; farx(n+1)=1.096930e+001; foe(n+1)=2.021196e+001; krok(n+1)=4.502543e-003; ng(n+1)=1.651276e+003;
n=36; farx(n+1)=6.978702e+000; foe(n+1)=1.477106e+001; krok(n+1)=5.915528e-003; ng(n+1)=1.563998e+003;
n=37; farx(n+1)=6.438599e+000; foe(n+1)=1.388139e+001; krok(n+1)=5.282547e-003; ng(n+1)=6.319201e+002;
n=38; farx(n+1)=5.340596e+000; foe(n+1)=1.167846e+001; krok(n+1)=2.827577e-003; ng(n+1)=1.911010e+003;
n=39; farx(n+1)=4.818199e+000; foe(n+1)=1.100947e+001; krok(n+1)=3.240245e-003; ng(n+1)=6.292085e+002;
n=40; farx(n+1)=3.994001e+000; foe(n+1)=9.649466e+000; krok(n+1)=7.251821e-003; ng(n+1)=1.615109e+003;
n=41; farx(n+1)=3.530898e+000; foe(n+1)=8.854263e+000; krok(n+1)=2.561707e-003; ng(n+1)=1.326844e+003;
n=42; farx(n+1)=2.834544e+000; foe(n+1)=8.160394e+000; krok(n+1)=4.026836e-003; ng(n+1)=4.199929e+002;
n=43; farx(n+1)=2.590929e+000; foe(n+1)=7.857869e+000; krok(n+1)=1.285709e-002; ng(n+1)=2.414498e+002;
n=44; farx(n+1)=2.306910e+000; foe(n+1)=7.353678e+000; krok(n+1)=2.118799e-002; ng(n+1)=4.779390e+002;
n=45; farx(n+1)=2.245519e+000; foe(n+1)=7.177017e+000; krok(n+1)=1.333847e-002; ng(n+1)=2.896837e+002;
n=46; farx(n+1)=2.092151e+000; foe(n+1)=6.677691e+000; krok(n+1)=4.663645e-002; ng(n+1)=4.755379e+002;
n=47; farx(n+1)=2.054416e+000; foe(n+1)=6.495011e+000; krok(n+1)=1.923658e-002; ng(n+1)=2.824815e+002;
n=48; farx(n+1)=2.058093e+000; foe(n+1)=6.382424e+000; krok(n+1)=1.219985e-002; ng(n+1)=3.753377e+002;
n=49; farx(n+1)=2.032059e+000; foe(n+1)=6.236357e+000; krok(n+1)=2.710756e-002; ng(n+1)=2.413591e+002;
n=50; farx(n+1)=1.901399e+000; foe(n+1)=5.897168e+000; krok(n+1)=5.026897e-002; ng(n+1)=6.736496e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=1.886095e+000; foe(n+1)=5.836905e+000; krok(n+1)=4.942298e-006; ng(n+1)=5.384024e+002;
n=52; farx(n+1)=1.886823e+000; foe(n+1)=5.662038e+000; krok(n+1)=6.703667e-005; ng(n+1)=1.934688e+002;
n=53; farx(n+1)=1.830052e+000; foe(n+1)=5.346383e+000; krok(n+1)=8.836179e-005; ng(n+1)=2.683074e+002;
n=54; farx(n+1)=1.829177e+000; foe(n+1)=5.314279e+000; krok(n+1)=2.357491e-005; ng(n+1)=1.883399e+002;
n=55; farx(n+1)=1.807610e+000; foe(n+1)=5.101480e+000; krok(n+1)=2.596069e-004; ng(n+1)=1.538697e+002;
n=56; farx(n+1)=1.783455e+000; foe(n+1)=4.954890e+000; krok(n+1)=5.817002e-004; ng(n+1)=2.252913e+002;
n=57; farx(n+1)=1.707509e+000; foe(n+1)=4.675083e+000; krok(n+1)=2.241713e-003; ng(n+1)=3.940845e+002;
n=58; farx(n+1)=1.685188e+000; foe(n+1)=4.512889e+000; krok(n+1)=5.167850e-003; ng(n+1)=7.850319e+002;
n=59; farx(n+1)=1.565139e+000; foe(n+1)=4.125706e+000; krok(n+1)=3.851383e-003; ng(n+1)=9.889506e+002;
n=60; farx(n+1)=1.552704e+000; foe(n+1)=3.860465e+000; krok(n+1)=1.321526e-002; ng(n+1)=2.954441e+002;
n=61; farx(n+1)=1.512197e+000; foe(n+1)=3.677333e+000; krok(n+1)=3.364099e-003; ng(n+1)=2.914515e+002;
n=62; farx(n+1)=1.467839e+000; foe(n+1)=3.594503e+000; krok(n+1)=3.758726e-003; ng(n+1)=3.689478e+002;
n=63; farx(n+1)=1.417003e+000; foe(n+1)=3.504441e+000; krok(n+1)=2.293321e-003; ng(n+1)=4.300807e+002;
n=64; farx(n+1)=1.233445e+000; foe(n+1)=3.215451e+000; krok(n+1)=3.360874e-002; ng(n+1)=2.978896e+002;
n=65; farx(n+1)=1.046203e+000; foe(n+1)=2.754976e+000; krok(n+1)=1.023953e-002; ng(n+1)=5.762141e+002;
n=66; farx(n+1)=1.019262e+000; foe(n+1)=2.679028e+000; krok(n+1)=7.088910e-003; ng(n+1)=3.098165e+002;
n=67; farx(n+1)=9.799532e-001; foe(n+1)=2.579132e+000; krok(n+1)=6.324208e-003; ng(n+1)=4.217780e+002;
n=68; farx(n+1)=9.846879e-001; foe(n+1)=2.511459e+000; krok(n+1)=9.428457e-003; ng(n+1)=2.280790e+002;
n=69; farx(n+1)=9.683155e-001; foe(n+1)=2.428478e+000; krok(n+1)=9.830492e-003; ng(n+1)=6.254592e+002;
n=70; farx(n+1)=9.130008e-001; foe(n+1)=2.247242e+000; krok(n+1)=3.703353e-002; ng(n+1)=5.415550e+002;
n=71; farx(n+1)=8.562446e-001; foe(n+1)=2.156950e+000; krok(n+1)=9.949946e-003; ng(n+1)=3.926221e+002;
n=72; farx(n+1)=8.237147e-001; foe(n+1)=2.037196e+000; krok(n+1)=1.352965e-002; ng(n+1)=2.757358e+002;
n=73; farx(n+1)=7.908001e-001; foe(n+1)=1.840209e+000; krok(n+1)=2.502768e-002; ng(n+1)=6.828825e+002;
n=74; farx(n+1)=7.615955e-001; foe(n+1)=1.714227e+000; krok(n+1)=5.753341e-002; ng(n+1)=5.393018e+002;
n=75; farx(n+1)=7.576309e-001; foe(n+1)=1.637162e+000; krok(n+1)=3.758726e-003; ng(n+1)=6.896181e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=7.551114e-001; foe(n+1)=1.596237e+000; krok(n+1)=1.438270e-006; ng(n+1)=6.762296e+002;
n=77; farx(n+1)=7.550892e-001; foe(n+1)=1.585556e+000; krok(n+1)=3.880465e-006; ng(n+1)=2.389353e+002;
n=78; farx(n+1)=7.536010e-001; foe(n+1)=1.506520e+000; krok(n+1)=6.783250e-005; ng(n+1)=1.615325e+002;
n=79; farx(n+1)=7.549162e-001; foe(n+1)=1.480992e+000; krok(n+1)=1.787641e-005; ng(n+1)=1.623882e+002;
n=80; farx(n+1)=7.553228e-001; foe(n+1)=1.468407e+000; krok(n+1)=1.016603e-004; ng(n+1)=5.025283e+001;
n=81; farx(n+1)=7.462225e-001; foe(n+1)=1.438750e+000; krok(n+1)=1.213392e-003; ng(n+1)=3.474620e+001;
n=82; farx(n+1)=7.448076e-001; foe(n+1)=1.423100e+000; krok(n+1)=1.772025e-003; ng(n+1)=6.022656e+001;
n=83; farx(n+1)=6.809877e-001; foe(n+1)=1.236389e+000; krok(n+1)=7.453659e-003; ng(n+1)=5.000874e+001;
n=84; farx(n+1)=6.742649e-001; foe(n+1)=1.181709e+000; krok(n+1)=1.505777e-003; ng(n+1)=2.203370e+002;
n=85; farx(n+1)=6.649236e-001; foe(n+1)=1.154230e+000; krok(n+1)=1.588254e-003; ng(n+1)=3.224444e+002;
n=86; farx(n+1)=6.591535e-001; foe(n+1)=1.119436e+000; krok(n+1)=5.167850e-003; ng(n+1)=2.349445e+002;
n=87; farx(n+1)=6.075013e-001; foe(n+1)=1.084223e+000; krok(n+1)=1.532510e-002; ng(n+1)=3.032267e+002;
n=88; farx(n+1)=5.778300e-001; foe(n+1)=1.058532e+000; krok(n+1)=5.557598e-003; ng(n+1)=1.505422e+002;
n=89; farx(n+1)=5.201605e-001; foe(n+1)=1.002715e+000; krok(n+1)=2.011176e-002; ng(n+1)=5.432067e+002;
n=90; farx(n+1)=4.821369e-001; foe(n+1)=9.553318e-001; krok(n+1)=8.579668e-003; ng(n+1)=7.021654e+002;
n=91; farx(n+1)=4.686930e-001; foe(n+1)=9.258768e-001; krok(n+1)=6.349841e-003; ng(n+1)=2.968807e+002;
n=92; farx(n+1)=4.710893e-001; foe(n+1)=8.903132e-001; krok(n+1)=7.585588e-003; ng(n+1)=2.907202e+002;
n=93; farx(n+1)=4.755164e-001; foe(n+1)=8.649139e-001; krok(n+1)=1.799391e-002; ng(n+1)=3.079270e+002;
n=94; farx(n+1)=4.710539e-001; foe(n+1)=8.512045e-001; krok(n+1)=2.219272e-002; ng(n+1)=3.361798e+002;
n=95; farx(n+1)=4.568173e-001; foe(n+1)=8.255255e-001; krok(n+1)=2.602503e-002; ng(n+1)=3.215984e+002;
n=96; farx(n+1)=4.398277e-001; foe(n+1)=8.124145e-001; krok(n+1)=8.625805e-003; ng(n+1)=4.873621e+002;
n=97; farx(n+1)=4.428600e-001; foe(n+1)=7.758076e-001; krok(n+1)=6.257219e-002; ng(n+1)=2.754360e+002;
n=98; farx(n+1)=4.492691e-001; foe(n+1)=7.392025e-001; krok(n+1)=4.509885e-002; ng(n+1)=2.460658e+002;
n=99; farx(n+1)=4.382707e-001; foe(n+1)=6.877233e-001; krok(n+1)=1.467726e-001; ng(n+1)=2.442710e+002;
n=100; farx(n+1)=4.343018e-001; foe(n+1)=6.749646e-001; krok(n+1)=7.709486e-002; ng(n+1)=2.970473e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
