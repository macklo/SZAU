%uczenie predyktora arx
clear all;
n=0; farx(n+1)=4.795438e+003; foe(n+1)=4.817033e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.123275e+003; foe(n+1)=3.834763e+003; krok(n+1)=3.697036e-004; ng(n+1)=8.974725e+003;
n=2; farx(n+1)=1.715068e+003; foe(n+1)=1.055426e+004; krok(n+1)=2.247399e-004; ng(n+1)=1.044602e+004;
n=3; farx(n+1)=8.525207e+002; foe(n+1)=1.684471e+004; krok(n+1)=2.625683e-004; ng(n+1)=8.953673e+003;
n=4; farx(n+1)=5.080513e+002; foe(n+1)=1.236924e+004; krok(n+1)=1.846469e-003; ng(n+1)=6.756993e+003;
n=5; farx(n+1)=3.829786e+002; foe(n+1)=9.700805e+003; krok(n+1)=1.203868e-003; ng(n+1)=3.618994e+003;
n=6; farx(n+1)=2.505340e+002; foe(n+1)=9.841274e+003; krok(n+1)=1.863415e-003; ng(n+1)=1.531871e+003;
n=7; farx(n+1)=2.104565e+002; foe(n+1)=9.437062e+003; krok(n+1)=7.564534e-004; ng(n+1)=1.563315e+003;
n=8; farx(n+1)=5.955792e+001; foe(n+1)=5.616021e+003; krok(n+1)=6.130766e-003; ng(n+1)=9.772964e+002;
n=9; farx(n+1)=2.201683e+001; foe(n+1)=8.315002e+002; krok(n+1)=1.983052e-003; ng(n+1)=1.776148e+003;
n=10; farx(n+1)=1.735716e+001; foe(n+1)=5.637240e+002; krok(n+1)=3.376366e-004; ng(n+1)=8.769060e+002;
n=11; farx(n+1)=8.860625e+000; foe(n+1)=1.101169e+002; krok(n+1)=6.399091e-003; ng(n+1)=2.491960e+002;
n=12; farx(n+1)=5.333855e+000; foe(n+1)=1.067982e+002; krok(n+1)=3.937497e-003; ng(n+1)=5.040711e+002;
n=13; farx(n+1)=3.235211e+000; foe(n+1)=1.033428e+002; krok(n+1)=1.165911e-002; ng(n+1)=3.997748e+002;
n=14; farx(n+1)=2.265025e+000; foe(n+1)=7.173040e+001; krok(n+1)=1.750226e-002; ng(n+1)=3.305452e+002;
n=15; farx(n+1)=2.005546e+000; foe(n+1)=7.125502e+001; krok(n+1)=2.001452e-002; ng(n+1)=1.300975e+002;
n=16; farx(n+1)=1.648562e+000; foe(n+1)=7.068362e+001; krok(n+1)=1.032026e-001; ng(n+1)=6.053352e+001;
n=17; farx(n+1)=1.386448e+000; foe(n+1)=6.431253e+001; krok(n+1)=9.685822e-002; ng(n+1)=9.607505e+001;
n=18; farx(n+1)=1.274392e+000; foe(n+1)=5.982131e+001; krok(n+1)=4.633705e-002; ng(n+1)=8.896080e+001;
n=19; farx(n+1)=1.214575e+000; foe(n+1)=5.640298e+001; krok(n+1)=4.828142e-002; ng(n+1)=5.894558e+001;
n=20; farx(n+1)=1.135909e+000; foe(n+1)=6.418918e+001; krok(n+1)=4.285764e-002; ng(n+1)=5.971347e+001;
n=21; farx(n+1)=9.286290e-001; foe(n+1)=8.931315e+001; krok(n+1)=4.572801e-001; ng(n+1)=5.742943e+001;
n=22; farx(n+1)=8.681899e-001; foe(n+1)=9.146443e+001; krok(n+1)=9.963479e-002; ng(n+1)=6.960436e+001;
n=23; farx(n+1)=7.700240e-001; foe(n+1)=1.186019e+002; krok(n+1)=2.243785e-001; ng(n+1)=6.159280e+001;
n=24; farx(n+1)=7.138738e-001; foe(n+1)=6.887579e+001; krok(n+1)=1.447238e-001; ng(n+1)=2.708656e+001;
n=25; farx(n+1)=6.649596e-001; foe(n+1)=5.722363e+001; krok(n+1)=1.228837e-001; ng(n+1)=2.215088e+001;
%odnowa zmiennej metryki
n=26; farx(n+1)=6.552315e-001; foe(n+1)=5.759352e+001; krok(n+1)=1.131625e-004; ng(n+1)=3.635658e+001;
n=27; farx(n+1)=6.476106e-001; foe(n+1)=5.979870e+001; krok(n+1)=2.218843e-003; ng(n+1)=7.977004e+000;
n=28; farx(n+1)=6.411466e-001; foe(n+1)=6.161607e+001; krok(n+1)=1.775310e-004; ng(n+1)=2.328048e+001;
n=29; farx(n+1)=6.187099e-001; foe(n+1)=7.883982e+001; krok(n+1)=5.368339e-003; ng(n+1)=8.739334e+000;
n=30; farx(n+1)=5.778598e-001; foe(n+1)=6.032629e+001; krok(n+1)=1.611092e-002; ng(n+1)=1.013763e+001;
n=31; farx(n+1)=5.693688e-001; foe(n+1)=5.140261e+001; krok(n+1)=2.173419e-002; ng(n+1)=4.843518e+001;
n=32; farx(n+1)=5.417752e-001; foe(n+1)=4.657645e+001; krok(n+1)=7.453607e-002; ng(n+1)=4.858057e+001;
n=33; farx(n+1)=5.006146e-001; foe(n+1)=4.403481e+001; krok(n+1)=6.211076e-002; ng(n+1)=4.129635e+001;
n=34; farx(n+1)=4.774185e-001; foe(n+1)=4.901198e+001; krok(n+1)=1.124683e-001; ng(n+1)=2.735074e+001;
n=35; farx(n+1)=4.546848e-001; foe(n+1)=4.057410e+001; krok(n+1)=1.230771e-001; ng(n+1)=2.172919e+001;
n=36; farx(n+1)=4.418494e-001; foe(n+1)=3.754976e+001; krok(n+1)=6.045742e-002; ng(n+1)=3.054993e+001;
n=37; farx(n+1)=4.193252e-001; foe(n+1)=2.953232e+001; krok(n+1)=3.856870e-001; ng(n+1)=2.291429e+001;
n=38; farx(n+1)=4.044848e-001; foe(n+1)=2.561846e+001; krok(n+1)=1.196783e-001; ng(n+1)=2.823636e+001;
n=39; farx(n+1)=3.851767e-001; foe(n+1)=2.703337e+001; krok(n+1)=4.764226e-001; ng(n+1)=2.096762e+001;
n=40; farx(n+1)=3.716302e-001; foe(n+1)=2.231751e+001; krok(n+1)=1.956250e-001; ng(n+1)=1.938101e+001;
n=41; farx(n+1)=3.615748e-001; foe(n+1)=2.530982e+001; krok(n+1)=1.902702e-001; ng(n+1)=2.498098e+001;
n=42; farx(n+1)=3.441671e-001; foe(n+1)=2.266104e+001; krok(n+1)=5.403833e-001; ng(n+1)=2.001437e+001;
n=43; farx(n+1)=3.343208e-001; foe(n+1)=2.563360e+001; krok(n+1)=3.505976e-001; ng(n+1)=2.269258e+001;
n=44; farx(n+1)=3.173482e-001; foe(n+1)=1.572272e+001; krok(n+1)=3.819012e-001; ng(n+1)=1.003520e+001;
n=45; farx(n+1)=3.132161e-001; foe(n+1)=1.402138e+001; krok(n+1)=1.961845e-001; ng(n+1)=4.666852e+000;
n=46; farx(n+1)=3.083477e-001; foe(n+1)=1.495886e+001; krok(n+1)=2.688700e-001; ng(n+1)=1.803375e+001;
n=47; farx(n+1)=3.040490e-001; foe(n+1)=2.059994e+001; krok(n+1)=2.192961e-001; ng(n+1)=7.882521e+000;
n=48; farx(n+1)=3.009275e-001; foe(n+1)=2.103325e+001; krok(n+1)=1.634459e-001; ng(n+1)=1.106337e+001;
n=49; farx(n+1)=2.991959e-001; foe(n+1)=1.976657e+001; krok(n+1)=1.892969e-001; ng(n+1)=1.205898e+001;
n=50; farx(n+1)=2.930008e-001; foe(n+1)=2.087000e+001; krok(n+1)=8.551844e-001; ng(n+1)=8.754429e+000;
%odnowa zmiennej metryki
n=51; farx(n+1)=2.921907e-001; foe(n+1)=1.852686e+001; krok(n+1)=1.046887e-004; ng(n+1)=1.072803e+001;
n=52; farx(n+1)=2.921112e-001; foe(n+1)=1.846648e+001; krok(n+1)=1.874926e-004; ng(n+1)=2.700054e+000;
n=53; farx(n+1)=2.917265e-001; foe(n+1)=1.718042e+001; krok(n+1)=9.053002e-004; ng(n+1)=3.113430e+000;
n=54; farx(n+1)=2.909161e-001; foe(n+1)=1.548228e+001; krok(n+1)=9.416683e-003; ng(n+1)=1.371299e+000;
n=55; farx(n+1)=2.905481e-001; foe(n+1)=1.585514e+001; krok(n+1)=6.854056e-003; ng(n+1)=1.152785e+000;
n=56; farx(n+1)=2.892541e-001; foe(n+1)=1.739754e+001; krok(n+1)=5.019954e-002; ng(n+1)=7.816933e-001;
n=57; farx(n+1)=2.880673e-001; foe(n+1)=1.477753e+001; krok(n+1)=4.215364e-002; ng(n+1)=8.863117e-001;
n=58; farx(n+1)=2.865294e-001; foe(n+1)=1.194693e+001; krok(n+1)=2.114442e-001; ng(n+1)=2.143078e+000;
n=59; farx(n+1)=2.821898e-001; foe(n+1)=9.465804e+000; krok(n+1)=4.641165e-001; ng(n+1)=3.742224e+000;
n=60; farx(n+1)=2.751727e-001; foe(n+1)=5.793500e+000; krok(n+1)=2.357265e-001; ng(n+1)=1.161056e+001;
n=61; farx(n+1)=2.710823e-001; foe(n+1)=4.397637e+000; krok(n+1)=2.134155e-001; ng(n+1)=7.985344e+000;
n=62; farx(n+1)=2.684170e-001; foe(n+1)=3.529766e+000; krok(n+1)=1.677326e-001; ng(n+1)=1.714633e+001;
n=63; farx(n+1)=2.645340e-001; foe(n+1)=3.257492e+000; krok(n+1)=1.497077e-001; ng(n+1)=8.910762e+000;
n=64; farx(n+1)=2.613073e-001; foe(n+1)=3.007905e+000; krok(n+1)=1.269153e-001; ng(n+1)=1.109823e+001;
n=65; farx(n+1)=2.572316e-001; foe(n+1)=2.413803e+000; krok(n+1)=6.744582e-001; ng(n+1)=1.560551e+001;
n=66; farx(n+1)=2.555105e-001; foe(n+1)=2.205930e+000; krok(n+1)=3.636509e-001; ng(n+1)=6.595236e+000;
n=67; farx(n+1)=2.525911e-001; foe(n+1)=2.276458e+000; krok(n+1)=3.647485e-001; ng(n+1)=4.996600e+000;
n=68; farx(n+1)=2.511686e-001; foe(n+1)=2.361384e+000; krok(n+1)=5.989693e-001; ng(n+1)=5.762120e+000;
n=69; farx(n+1)=2.494526e-001; foe(n+1)=2.493803e+000; krok(n+1)=6.776371e-001; ng(n+1)=9.481153e+000;
n=70; farx(n+1)=2.480576e-001; foe(n+1)=2.499454e+000; krok(n+1)=7.631010e-001; ng(n+1)=3.378051e+000;
n=71; farx(n+1)=2.458841e-001; foe(n+1)=2.486618e+000; krok(n+1)=8.131859e-001; ng(n+1)=4.253260e+000;
n=72; farx(n+1)=2.436131e-001; foe(n+1)=2.661711e+000; krok(n+1)=5.714050e-001; ng(n+1)=2.348457e+000;
n=73; farx(n+1)=2.395474e-001; foe(n+1)=2.583434e+000; krok(n+1)=4.498732e-001; ng(n+1)=9.534549e+000;
n=74; farx(n+1)=2.372497e-001; foe(n+1)=2.620895e+000; krok(n+1)=1.710829e-001; ng(n+1)=9.325020e+000;
n=75; farx(n+1)=2.354851e-001; foe(n+1)=2.503629e+000; krok(n+1)=3.019693e-001; ng(n+1)=7.315450e+000;
%odnowa zmiennej metryki
n=76; farx(n+1)=2.351150e-001; foe(n+1)=2.528113e+000; krok(n+1)=1.168734e-004; ng(n+1)=7.865845e+000;
n=77; farx(n+1)=2.349079e-001; foe(n+1)=2.518065e+000; krok(n+1)=2.707062e-003; ng(n+1)=1.349176e+000;
n=78; farx(n+1)=2.348035e-001; foe(n+1)=2.516525e+000; krok(n+1)=1.185230e-004; ng(n+1)=3.514379e+000;
n=79; farx(n+1)=2.346384e-001; foe(n+1)=2.577736e+000; krok(n+1)=2.774090e-003; ng(n+1)=1.116442e+000;
n=80; farx(n+1)=2.344231e-001; foe(n+1)=2.598954e+000; krok(n+1)=9.564855e-003; ng(n+1)=7.219274e-001;
n=81; farx(n+1)=2.341259e-001; foe(n+1)=2.517149e+000; krok(n+1)=5.286105e-002; ng(n+1)=3.685082e-001;
n=82; farx(n+1)=2.338624e-001; foe(n+1)=2.509096e+000; krok(n+1)=5.421512e-002; ng(n+1)=3.393540e-001;
n=83; farx(n+1)=2.333169e-001; foe(n+1)=2.609120e+000; krok(n+1)=1.191056e-001; ng(n+1)=4.417463e-001;
n=84; farx(n+1)=2.331096e-001; foe(n+1)=2.623233e+000; krok(n+1)=1.365450e-001; ng(n+1)=1.683135e+000;
n=85; farx(n+1)=2.327840e-001; foe(n+1)=2.508182e+000; krok(n+1)=2.183343e-001; ng(n+1)=2.927475e+000;
n=86; farx(n+1)=2.322618e-001; foe(n+1)=2.476234e+000; krok(n+1)=8.847528e-001; ng(n+1)=5.028591e+000;
n=87; farx(n+1)=2.315993e-001; foe(n+1)=2.507106e+000; krok(n+1)=6.193808e-001; ng(n+1)=2.030387e+000;
n=88; farx(n+1)=2.308021e-001; foe(n+1)=2.470058e+000; krok(n+1)=3.805404e-001; ng(n+1)=6.206110e+000;
n=89; farx(n+1)=2.296236e-001; foe(n+1)=2.499692e+000; krok(n+1)=8.155462e-001; ng(n+1)=7.064360e+000;
n=90; farx(n+1)=2.290193e-001; foe(n+1)=2.408183e+000; krok(n+1)=2.343012e-001; ng(n+1)=8.196128e+000;
n=91; farx(n+1)=2.284708e-001; foe(n+1)=2.401192e+000; krok(n+1)=6.115403e-001; ng(n+1)=1.396638e+000;
n=92; farx(n+1)=2.278439e-001; foe(n+1)=2.424811e+000; krok(n+1)=5.450576e-001; ng(n+1)=5.341692e+000;
n=93; farx(n+1)=2.275970e-001; foe(n+1)=2.484434e+000; krok(n+1)=4.556164e-001; ng(n+1)=7.521393e-001;
n=94; farx(n+1)=2.273410e-001; foe(n+1)=2.456300e+000; krok(n+1)=1.106940e+000; ng(n+1)=2.072772e+000;
n=95; farx(n+1)=2.271208e-001; foe(n+1)=2.465425e+000; krok(n+1)=6.986165e-001; ng(n+1)=3.768180e+000;
n=96; farx(n+1)=2.266125e-001; foe(n+1)=2.455032e+000; krok(n+1)=1.197662e+000; ng(n+1)=3.826874e+000;
n=97; farx(n+1)=2.263607e-001; foe(n+1)=2.505224e+000; krok(n+1)=4.787134e-001; ng(n+1)=2.869308e+000;
n=98; farx(n+1)=2.260578e-001; foe(n+1)=2.585958e+000; krok(n+1)=7.461552e-001; ng(n+1)=4.399079e+000;
n=99; farx(n+1)=2.258321e-001; foe(n+1)=2.592255e+000; krok(n+1)=1.507294e+000; ng(n+1)=2.246357e+000;
n=100; farx(n+1)=2.255680e-001; foe(n+1)=2.642451e+000; krok(n+1)=1.879567e+000; ng(n+1)=2.344214e+000;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora ARX');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
