%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.537746e+003; foe(n+1)=4.552377e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.807449e+003; foe(n+1)=2.885759e+003; krok(n+1)=4.863026e-004; ng(n+1)=4.360486e+003;
n=2; farx(n+1)=1.066313e+003; foe(n+1)=1.029176e+003; krok(n+1)=2.586841e-003; ng(n+1)=2.350658e+003;
n=3; farx(n+1)=1.284968e+003; foe(n+1)=6.462053e+002; krok(n+1)=1.954789e-004; ng(n+1)=7.830025e+003;
n=4; farx(n+1)=1.611247e+003; foe(n+1)=5.360461e+002; krok(n+1)=5.970040e-004; ng(n+1)=5.336901e+003;
n=5; farx(n+1)=6.453912e+002; foe(n+1)=3.720381e+002; krok(n+1)=8.427741e-003; ng(n+1)=5.431444e+002;
n=6; farx(n+1)=5.520579e+002; foe(n+1)=3.585088e+002; krok(n+1)=1.828869e-004; ng(n+1)=1.297648e+003;
n=7; farx(n+1)=3.683305e+002; foe(n+1)=3.153609e+002; krok(n+1)=6.352069e-005; ng(n+1)=4.462622e+003;
n=8; farx(n+1)=3.134731e+002; foe(n+1)=3.036230e+002; krok(n+1)=1.491801e-003; ng(n+1)=1.516858e+003;
n=9; farx(n+1)=1.339000e+002; foe(n+1)=2.599277e+002; krok(n+1)=2.741623e-002; ng(n+1)=4.391492e+002;
n=10; farx(n+1)=1.310305e+002; foe(n+1)=2.591345e+002; krok(n+1)=3.318147e-005; ng(n+1)=1.607553e+003;
n=11; farx(n+1)=1.200707e+002; foe(n+1)=2.568078e+002; krok(n+1)=1.282071e-003; ng(n+1)=9.821551e+002;
n=12; farx(n+1)=1.152272e+002; foe(n+1)=2.525605e+002; krok(n+1)=1.109421e-003; ng(n+1)=8.830043e+002;
n=13; farx(n+1)=1.062803e+002; foe(n+1)=2.489446e+002; krok(n+1)=9.208067e-004; ng(n+1)=1.287882e+003;
n=14; farx(n+1)=9.156899e+001; foe(n+1)=2.427234e+002; krok(n+1)=8.687921e-004; ng(n+1)=3.873033e+003;
n=15; farx(n+1)=7.513973e+001; foe(n+1)=2.336570e+002; krok(n+1)=1.181843e-003; ng(n+1)=1.270986e+003;
n=16; farx(n+1)=6.238418e+001; foe(n+1)=2.282206e+002; krok(n+1)=2.102562e-004; ng(n+1)=6.417884e+003;
n=17; farx(n+1)=5.978167e+001; foe(n+1)=2.263301e+002; krok(n+1)=6.254538e-004; ng(n+1)=6.363487e+003;
n=18; farx(n+1)=5.744544e+001; foe(n+1)=2.248070e+002; krok(n+1)=1.969720e-004; ng(n+1)=4.855054e+003;
n=19; farx(n+1)=5.098752e+001; foe(n+1)=2.185322e+002; krok(n+1)=2.914778e-003; ng(n+1)=1.174810e+004;
n=20; farx(n+1)=4.998712e+001; foe(n+1)=2.173203e+002; krok(n+1)=3.081200e-004; ng(n+1)=3.535787e+003;
n=21; farx(n+1)=4.849869e+001; foe(n+1)=2.124654e+002; krok(n+1)=1.812955e-003; ng(n+1)=9.151323e+003;
n=22; farx(n+1)=5.166205e+001; foe(n+1)=2.070104e+002; krok(n+1)=3.764443e-004; ng(n+1)=1.171017e+004;
n=23; farx(n+1)=5.306795e+001; foe(n+1)=2.061267e+002; krok(n+1)=3.909299e-004; ng(n+1)=8.713120e+003;
n=24; farx(n+1)=7.071165e+001; foe(n+1)=1.953240e+002; krok(n+1)=4.611188e-003; ng(n+1)=2.905223e+003;
n=25; farx(n+1)=7.269184e+001; foe(n+1)=1.946788e+002; krok(n+1)=7.188658e-005; ng(n+1)=1.109545e+004;
%odnowa zmiennej metryki
n=26; farx(n+1)=7.257151e+001; foe(n+1)=1.945132e+002; krok(n+1)=1.831484e-006; ng(n+1)=1.297523e+003;
n=27; farx(n+1)=6.873456e+001; foe(n+1)=1.730730e+002; krok(n+1)=1.104522e-005; ng(n+1)=5.311459e+003;
n=28; farx(n+1)=3.764089e+001; foe(n+1)=1.206716e+002; krok(n+1)=1.915864e-004; ng(n+1)=6.209548e+003;
n=29; farx(n+1)=2.956777e+001; foe(n+1)=1.026097e+002; krok(n+1)=4.651747e-004; ng(n+1)=1.784311e+003;
n=30; farx(n+1)=2.415652e+001; foe(n+1)=8.620894e+001; krok(n+1)=8.336543e-004; ng(n+1)=1.490655e+003;
n=31; farx(n+1)=1.759374e+001; foe(n+1)=6.045198e+001; krok(n+1)=9.576184e-004; ng(n+1)=2.309116e+003;
n=32; farx(n+1)=1.817110e+001; foe(n+1)=5.184950e+001; krok(n+1)=8.733921e-004; ng(n+1)=2.280038e+003;
n=33; farx(n+1)=2.087515e+001; foe(n+1)=4.437662e+001; krok(n+1)=1.034888e-003; ng(n+1)=1.602601e+003;
n=34; farx(n+1)=2.156321e+001; foe(n+1)=4.296378e+001; krok(n+1)=9.377529e-006; ng(n+1)=1.881828e+004;
n=35; farx(n+1)=2.290235e+001; foe(n+1)=3.909430e+001; krok(n+1)=1.670732e-004; ng(n+1)=1.988940e+004;
n=36; farx(n+1)=2.146080e+001; foe(n+1)=3.545755e+001; krok(n+1)=3.206931e-003; ng(n+1)=1.068101e+003;
n=37; farx(n+1)=1.805867e+001; foe(n+1)=3.031639e+001; krok(n+1)=3.930173e-003; ng(n+1)=1.326412e+003;
n=38; farx(n+1)=1.182660e+001; foe(n+1)=2.374710e+001; krok(n+1)=1.219985e-002; ng(n+1)=7.225783e+002;
n=39; farx(n+1)=9.509911e+000; foe(n+1)=2.030859e+001; krok(n+1)=4.290347e-003; ng(n+1)=6.871070e+002;
n=40; farx(n+1)=9.320407e+000; foe(n+1)=1.969068e+001; krok(n+1)=2.996980e-004; ng(n+1)=1.641169e+003;
n=41; farx(n+1)=8.586148e+000; foe(n+1)=1.732497e+001; krok(n+1)=1.214548e-002; ng(n+1)=5.774730e+002;
n=42; farx(n+1)=8.473075e+000; foe(n+1)=1.687234e+001; krok(n+1)=8.300976e-003; ng(n+1)=4.972082e+002;
n=43; farx(n+1)=8.413657e+000; foe(n+1)=1.498944e+001; krok(n+1)=1.952078e-002; ng(n+1)=8.808526e+002;
n=44; farx(n+1)=8.531197e+000; foe(n+1)=1.445681e+001; krok(n+1)=1.204622e-002; ng(n+1)=3.388184e+002;
n=45; farx(n+1)=8.176719e+000; foe(n+1)=1.395759e+001; krok(n+1)=2.680031e-002; ng(n+1)=3.380121e+002;
n=46; farx(n+1)=8.176402e+000; foe(n+1)=1.329302e+001; krok(n+1)=2.142359e-002; ng(n+1)=6.441233e+002;
n=47; farx(n+1)=7.849264e+000; foe(n+1)=1.300127e+001; krok(n+1)=1.605804e-002; ng(n+1)=3.805012e+002;
n=48; farx(n+1)=7.453479e+000; foe(n+1)=1.259176e+001; krok(n+1)=6.423218e-002; ng(n+1)=2.736202e+002;
n=49; farx(n+1)=6.826282e+000; foe(n+1)=1.193449e+001; krok(n+1)=7.224063e-002; ng(n+1)=3.147489e+002;
n=50; farx(n+1)=6.593109e+000; foe(n+1)=1.154440e+001; krok(n+1)=5.873647e-002; ng(n+1)=2.578747e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=6.593023e+000; foe(n+1)=1.136889e+001; krok(n+1)=4.645269e-005; ng(n+1)=2.805923e+002;
n=52; farx(n+1)=6.691674e+000; foe(n+1)=1.120422e+001; krok(n+1)=2.323045e-005; ng(n+1)=4.222837e+002;
n=53; farx(n+1)=6.047365e+000; foe(n+1)=1.067760e+001; krok(n+1)=6.518159e-004; ng(n+1)=1.268645e+002;
n=54; farx(n+1)=6.014910e+000; foe(n+1)=1.064834e+001; krok(n+1)=1.758667e-004; ng(n+1)=5.635458e+001;
n=55; farx(n+1)=5.786770e+000; foe(n+1)=1.041298e+001; krok(n+1)=2.909519e-003; ng(n+1)=4.574852e+001;
n=56; farx(n+1)=5.584035e+000; foe(n+1)=1.003706e+001; krok(n+1)=5.485375e-003; ng(n+1)=4.581621e+001;
n=57; farx(n+1)=5.286725e+000; foe(n+1)=9.724946e+000; krok(n+1)=5.711090e-003; ng(n+1)=1.784437e+002;
n=58; farx(n+1)=4.672345e+000; foe(n+1)=9.063932e+000; krok(n+1)=1.576669e-002; ng(n+1)=2.985368e+002;
n=59; farx(n+1)=4.125252e+000; foe(n+1)=8.702535e+000; krok(n+1)=7.162643e-003; ng(n+1)=4.377924e+002;
n=60; farx(n+1)=3.855568e+000; foe(n+1)=8.450674e+000; krok(n+1)=1.274938e-002; ng(n+1)=4.553456e+002;
n=61; farx(n+1)=3.554571e+000; foe(n+1)=8.152901e+000; krok(n+1)=2.074719e-002; ng(n+1)=1.840257e+002;
n=62; farx(n+1)=3.262594e+000; foe(n+1)=7.792554e+000; krok(n+1)=5.977052e-003; ng(n+1)=5.268413e+002;
n=63; farx(n+1)=3.138375e+000; foe(n+1)=7.517387e+000; krok(n+1)=9.248075e-003; ng(n+1)=3.952525e+002;
n=64; farx(n+1)=3.063547e+000; foe(n+1)=7.306104e+000; krok(n+1)=8.483440e-003; ng(n+1)=2.387475e+002;
n=65; farx(n+1)=3.057162e+000; foe(n+1)=7.214889e+000; krok(n+1)=7.119007e-003; ng(n+1)=1.868754e+002;
n=66; farx(n+1)=2.990760e+000; foe(n+1)=7.111318e+000; krok(n+1)=1.626130e-002; ng(n+1)=1.384428e+002;
n=67; farx(n+1)=2.860423e+000; foe(n+1)=6.860525e+000; krok(n+1)=3.043820e-002; ng(n+1)=1.955968e+002;
n=68; farx(n+1)=2.811506e+000; foe(n+1)=6.669455e+000; krok(n+1)=9.866494e-003; ng(n+1)=2.701224e+002;
n=69; farx(n+1)=2.799042e+000; foe(n+1)=6.461939e+000; krok(n+1)=1.401574e-002; ng(n+1)=3.716976e+002;
n=70; farx(n+1)=2.800835e+000; foe(n+1)=6.327182e+000; krok(n+1)=2.570124e-002; ng(n+1)=1.365555e+002;
n=71; farx(n+1)=2.706928e+000; foe(n+1)=5.970035e+000; krok(n+1)=1.110632e-001; ng(n+1)=2.779858e+002;
n=72; farx(n+1)=2.601126e+000; foe(n+1)=5.743373e+000; krok(n+1)=3.758985e-002; ng(n+1)=3.818544e+002;
n=73; farx(n+1)=2.524809e+000; foe(n+1)=5.577623e+000; krok(n+1)=8.333356e-002; ng(n+1)=1.995906e+002;
n=74; farx(n+1)=2.507803e+000; foe(n+1)=5.544184e+000; krok(n+1)=2.928765e-002; ng(n+1)=1.641934e+002;
n=75; farx(n+1)=2.452520e+000; foe(n+1)=5.451347e+000; krok(n+1)=3.416990e-002; ng(n+1)=1.471041e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=2.445059e+000; foe(n+1)=5.400272e+000; krok(n+1)=6.925078e-006; ng(n+1)=3.500363e+002;
n=77; farx(n+1)=2.438892e+000; foe(n+1)=5.335657e+000; krok(n+1)=8.876552e-005; ng(n+1)=1.340424e+002;
n=78; farx(n+1)=2.430532e+000; foe(n+1)=5.319022e+000; krok(n+1)=4.176829e-005; ng(n+1)=9.105984e+001;
n=79; farx(n+1)=2.428752e+000; foe(n+1)=5.286118e+000; krok(n+1)=1.784716e-004; ng(n+1)=6.313091e+001;
n=80; farx(n+1)=2.430473e+000; foe(n+1)=5.247443e+000; krok(n+1)=1.588254e-003; ng(n+1)=2.618479e+001;
n=81; farx(n+1)=2.349727e+000; foe(n+1)=5.082019e+000; krok(n+1)=6.583213e-003; ng(n+1)=2.616635e+001;
n=82; farx(n+1)=2.343626e+000; foe(n+1)=4.909077e+000; krok(n+1)=6.193383e-003; ng(n+1)=9.200818e+001;
n=83; farx(n+1)=2.314825e+000; foe(n+1)=4.837577e+000; krok(n+1)=2.959713e-002; ng(n+1)=3.470436e+002;
n=84; farx(n+1)=2.279521e+000; foe(n+1)=4.745539e+000; krok(n+1)=8.705227e-003; ng(n+1)=2.718888e+002;
n=85; farx(n+1)=2.219380e+000; foe(n+1)=4.709431e+000; krok(n+1)=1.071180e-002; ng(n+1)=2.369411e+002;
n=86; farx(n+1)=2.171588e+000; foe(n+1)=4.680067e+000; krok(n+1)=2.069473e-002; ng(n+1)=2.362571e+002;
n=87; farx(n+1)=2.155812e+000; foe(n+1)=4.659944e+000; krok(n+1)=1.734425e-002; ng(n+1)=2.102000e+002;
n=88; farx(n+1)=2.125880e+000; foe(n+1)=4.612428e+000; krok(n+1)=3.144138e-002; ng(n+1)=1.988289e+002;
n=89; farx(n+1)=2.110494e+000; foe(n+1)=4.540789e+000; krok(n+1)=2.753897e-002; ng(n+1)=1.436917e+002;
n=90; farx(n+1)=2.057135e+000; foe(n+1)=4.451848e+000; krok(n+1)=4.438543e-002; ng(n+1)=2.786563e+002;
n=91; farx(n+1)=2.018853e+000; foe(n+1)=4.352957e+000; krok(n+1)=2.733118e-002; ng(n+1)=3.309889e+002;
n=92; farx(n+1)=2.026709e+000; foe(n+1)=4.271595e+000; krok(n+1)=3.239658e-002; ng(n+1)=3.801681e+002;
n=93; farx(n+1)=2.020023e+000; foe(n+1)=4.212154e+000; krok(n+1)=1.360507e-002; ng(n+1)=3.451158e+002;
n=94; farx(n+1)=2.031134e+000; foe(n+1)=4.151498e+000; krok(n+1)=3.805310e-002; ng(n+1)=3.175751e+002;
n=95; farx(n+1)=2.029162e+000; foe(n+1)=4.030843e+000; krok(n+1)=1.036878e-001; ng(n+1)=2.187465e+002;
n=96; farx(n+1)=2.023016e+000; foe(n+1)=3.924114e+000; krok(n+1)=6.206136e-002; ng(n+1)=2.427396e+002;
n=97; farx(n+1)=1.861127e+000; foe(n+1)=3.734575e+000; krok(n+1)=4.486659e-002; ng(n+1)=2.280429e+002;
n=98; farx(n+1)=1.702286e+000; foe(n+1)=3.568339e+000; krok(n+1)=3.731451e-002; ng(n+1)=3.367563e+002;
n=99; farx(n+1)=1.651222e+000; foe(n+1)=3.489942e+000; krok(n+1)=4.050274e-002; ng(n+1)=3.285084e+002;
n=100; farx(n+1)=1.623948e+000; foe(n+1)=3.390438e+000; krok(n+1)=6.864555e-002; ng(n+1)=2.095291e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
