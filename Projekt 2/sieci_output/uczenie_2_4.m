%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.959700e+003; foe(n+1)=4.878933e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.897523e+003; foe(n+1)=3.856995e+003; krok(n+1)=3.505202e-004; ng(n+1)=2.663243e+003;
n=2; farx(n+1)=8.743063e+002; foe(n+1)=7.954593e+002; krok(n+1)=3.804775e-003; ng(n+1)=1.074787e+003;
n=3; farx(n+1)=8.894358e+002; foe(n+1)=7.075620e+002; krok(n+1)=1.992091e-004; ng(n+1)=1.802830e+003;
n=4; farx(n+1)=5.068220e+002; foe(n+1)=5.181733e+002; krok(n+1)=7.251821e-003; ng(n+1)=1.246527e+003;
n=5; farx(n+1)=4.963038e+002; foe(n+1)=4.901774e+002; krok(n+1)=2.032216e-003; ng(n+1)=7.996587e+002;
n=6; farx(n+1)=4.463899e+002; foe(n+1)=4.719406e+002; krok(n+1)=8.861138e-004; ng(n+1)=5.111766e+002;
n=7; farx(n+1)=4.538704e+002; foe(n+1)=4.497821e+002; krok(n+1)=2.687524e-002; ng(n+1)=4.658079e+002;
n=8; farx(n+1)=5.351382e+002; foe(n+1)=4.373503e+002; krok(n+1)=3.882853e-002; ng(n+1)=3.828353e+002;
n=9; farx(n+1)=5.327704e+002; foe(n+1)=4.356534e+002; krok(n+1)=1.810600e-003; ng(n+1)=3.916120e+002;
n=10; farx(n+1)=4.602165e+002; foe(n+1)=4.151428e+002; krok(n+1)=8.452792e-002; ng(n+1)=3.834006e+002;
n=11; farx(n+1)=4.596573e+002; foe(n+1)=4.118782e+002; krok(n+1)=1.241631e-001; ng(n+1)=3.575348e+002;
n=12; farx(n+1)=4.518181e+002; foe(n+1)=4.088451e+002; krok(n+1)=2.128042e-002; ng(n+1)=4.342402e+002;
n=13; farx(n+1)=3.720238e+002; foe(n+1)=3.861323e+002; krok(n+1)=2.701727e-001; ng(n+1)=3.802920e+002;
n=14; farx(n+1)=3.782794e+002; foe(n+1)=3.820927e+002; krok(n+1)=3.780408e-001; ng(n+1)=1.814689e+002;
n=15; farx(n+1)=1.077945e+002; foe(n+1)=3.170474e+002; krok(n+1)=4.703379e+000; ng(n+1)=2.250376e+002;
n=16; farx(n+1)=1.057403e+002; foe(n+1)=3.168629e+002; krok(n+1)=1.064566e-003; ng(n+1)=1.292986e+003;
n=17; farx(n+1)=8.313560e+001; foe(n+1)=2.494201e+002; krok(n+1)=1.034461e+000; ng(n+1)=1.523435e+003;
n=18; farx(n+1)=8.309925e+001; foe(n+1)=2.493901e+002; krok(n+1)=9.677161e-005; ng(n+1)=9.420103e+002;
n=19; farx(n+1)=5.356958e+001; foe(n+1)=2.045293e+002; krok(n+1)=3.380830e-001; ng(n+1)=8.699474e+002;
n=20; farx(n+1)=5.850243e+001; foe(n+1)=1.719909e+002; krok(n+1)=2.869153e-001; ng(n+1)=6.890841e+003;
n=21; farx(n+1)=5.608976e+001; foe(n+1)=1.689417e+002; krok(n+1)=1.058808e-002; ng(n+1)=4.826672e+003;
n=22; farx(n+1)=3.307126e+001; foe(n+1)=1.491723e+002; krok(n+1)=5.695205e-002; ng(n+1)=2.046748e+003;
n=23; farx(n+1)=3.423201e+001; foe(n+1)=1.303179e+002; krok(n+1)=1.068992e-001; ng(n+1)=7.978298e+003;
n=24; farx(n+1)=2.048112e+001; foe(n+1)=9.581567e+001; krok(n+1)=1.617795e+000; ng(n+1)=1.984219e+003;
n=25; farx(n+1)=2.062953e+001; foe(n+1)=9.553515e+001; krok(n+1)=6.699412e-003; ng(n+1)=4.176990e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=2.051505e+001; foe(n+1)=9.518576e+001; krok(n+1)=6.777058e-007; ng(n+1)=1.769896e+003;
n=27; farx(n+1)=1.947230e+001; foe(n+1)=9.204609e+001; krok(n+1)=2.302017e-004; ng(n+1)=2.876997e+002;
n=28; farx(n+1)=1.816979e+001; foe(n+1)=8.968084e+001; krok(n+1)=7.204982e-005; ng(n+1)=1.377256e+003;
n=29; farx(n+1)=1.810203e+001; foe(n+1)=8.671862e+001; krok(n+1)=2.394547e-004; ng(n+1)=3.845418e+002;
n=30; farx(n+1)=1.568076e+001; foe(n+1)=8.153200e+001; krok(n+1)=2.026610e-004; ng(n+1)=1.801242e+003;
n=31; farx(n+1)=1.419195e+001; foe(n+1)=8.020435e+001; krok(n+1)=3.625910e-003; ng(n+1)=4.134864e+003;
n=32; farx(n+1)=1.375855e+001; foe(n+1)=7.702481e+001; krok(n+1)=1.810600e-003; ng(n+1)=3.781233e+003;
n=33; farx(n+1)=1.385374e+001; foe(n+1)=7.227805e+001; krok(n+1)=2.287557e-003; ng(n+1)=4.061678e+003;
n=34; farx(n+1)=1.279465e+001; foe(n+1)=6.905371e+001; krok(n+1)=9.195867e-003; ng(n+1)=2.149587e+003;
n=35; farx(n+1)=1.314651e+001; foe(n+1)=6.712140e+001; krok(n+1)=3.616714e-002; ng(n+1)=1.426598e+003;
n=36; farx(n+1)=1.283695e+001; foe(n+1)=6.576867e+001; krok(n+1)=2.147336e-002; ng(n+1)=2.358849e+003;
n=37; farx(n+1)=1.343659e+001; foe(n+1)=6.233278e+001; krok(n+1)=2.409244e-002; ng(n+1)=3.453966e+003;
n=38; farx(n+1)=1.363812e+001; foe(n+1)=5.124006e+001; krok(n+1)=1.931257e-001; ng(n+1)=1.695252e+003;
n=39; farx(n+1)=9.305929e+000; foe(n+1)=4.230860e+001; krok(n+1)=9.372047e-001; ng(n+1)=2.043735e+003;
n=40; farx(n+1)=9.066509e+000; foe(n+1)=3.922065e+001; krok(n+1)=9.270274e-001; ng(n+1)=1.618647e+003;
n=41; farx(n+1)=9.000807e+000; foe(n+1)=3.692399e+001; krok(n+1)=6.477968e-001; ng(n+1)=3.987288e+002;
n=42; farx(n+1)=8.707993e+000; foe(n+1)=3.550008e+001; krok(n+1)=4.446475e-001; ng(n+1)=6.824622e+002;
n=43; farx(n+1)=9.017857e+000; foe(n+1)=3.444612e+001; krok(n+1)=4.577617e-001; ng(n+1)=4.763919e+002;
n=44; farx(n+1)=8.144164e+000; foe(n+1)=3.248933e+001; krok(n+1)=2.697833e+000; ng(n+1)=5.639130e+002;
n=45; farx(n+1)=8.066643e+000; foe(n+1)=3.186571e+001; krok(n+1)=1.015323e+000; ng(n+1)=6.345936e+002;
n=46; farx(n+1)=7.687517e+000; foe(n+1)=3.100220e+001; krok(n+1)=1.361947e+000; ng(n+1)=2.544477e+002;
n=47; farx(n+1)=7.170891e+000; foe(n+1)=3.059324e+001; krok(n+1)=1.210739e+000; ng(n+1)=6.002850e+002;
n=48; farx(n+1)=6.864380e+000; foe(n+1)=3.010349e+001; krok(n+1)=1.390988e+000; ng(n+1)=2.710022e+002;
n=49; farx(n+1)=6.541915e+000; foe(n+1)=2.983631e+001; krok(n+1)=1.084311e+000; ng(n+1)=1.415909e+002;
n=50; farx(n+1)=6.162353e+000; foe(n+1)=2.968458e+001; krok(n+1)=1.707324e+000; ng(n+1)=1.820533e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=6.152532e+000; foe(n+1)=2.967727e+001; krok(n+1)=7.599740e-006; ng(n+1)=7.251647e+001;
n=52; farx(n+1)=6.128572e+000; foe(n+1)=2.965136e+001; krok(n+1)=1.244425e-005; ng(n+1)=1.206511e+002;
n=53; farx(n+1)=6.118909e+000; foe(n+1)=2.965010e+001; krok(n+1)=4.525103e-005; ng(n+1)=1.473955e+001;
n=54; farx(n+1)=6.085978e+000; foe(n+1)=2.943378e+001; krok(n+1)=2.128042e-002; ng(n+1)=9.063633e+000;
n=55; farx(n+1)=6.038862e+000; foe(n+1)=2.942263e+001; krok(n+1)=4.187548e-004; ng(n+1)=8.340136e+001;
n=56; farx(n+1)=5.846506e+000; foe(n+1)=2.921508e+001; krok(n+1)=3.207144e-002; ng(n+1)=9.990792e+001;
n=57; farx(n+1)=5.607830e+000; foe(n+1)=2.909872e+001; krok(n+1)=3.195265e-003; ng(n+1)=2.421209e+002;
n=58; farx(n+1)=5.024741e+000; foe(n+1)=2.869898e+001; krok(n+1)=4.568872e-002; ng(n+1)=5.501125e+002;
n=59; farx(n+1)=4.676007e+000; foe(n+1)=2.863230e+001; krok(n+1)=2.407667e-002; ng(n+1)=3.760830e+002;
n=60; farx(n+1)=4.495413e+000; foe(n+1)=2.860503e+001; krok(n+1)=4.681906e-002; ng(n+1)=3.255845e+002;
n=61; farx(n+1)=4.583746e+000; foe(n+1)=2.856155e+001; krok(n+1)=3.307424e-001; ng(n+1)=1.101035e+002;
n=62; farx(n+1)=4.520932e+000; foe(n+1)=2.853982e+001; krok(n+1)=1.280929e+000; ng(n+1)=6.365299e+001;
n=63; farx(n+1)=4.345160e+000; foe(n+1)=2.843117e+001; krok(n+1)=5.909229e+000; ng(n+1)=7.963012e+001;
n=64; farx(n+1)=4.015665e+000; foe(n+1)=2.823576e+001; krok(n+1)=1.841069e+000; ng(n+1)=8.969632e+001;
n=65; farx(n+1)=4.040356e+000; foe(n+1)=2.812695e+001; krok(n+1)=5.281352e-001; ng(n+1)=4.637843e+002;
n=66; farx(n+1)=3.460895e+000; foe(n+1)=2.778584e+001; krok(n+1)=1.202875e+000; ng(n+1)=7.306239e+002;
n=67; farx(n+1)=3.470345e+000; foe(n+1)=2.741716e+001; krok(n+1)=9.073805e-001; ng(n+1)=1.407307e+002;
n=68; farx(n+1)=3.480716e+000; foe(n+1)=2.692186e+001; krok(n+1)=7.238598e-001; ng(n+1)=9.818320e+002;
n=69; farx(n+1)=3.172595e+000; foe(n+1)=2.635793e+001; krok(n+1)=1.646105e-001; ng(n+1)=9.963498e+002;
n=70; farx(n+1)=3.172595e+000; foe(n+1)=2.635793e+001; krok(n+1)=7.390764e-016; ng(n+1)=3.809171e+003;
n=71; farx(n+1)=3.172809e+000; foe(n+1)=2.634939e+001; krok(n+1)=8.103568e-006; ng(n+1)=3.809171e+003;
n=72; farx(n+1)=3.083951e+000; foe(n+1)=2.569067e+001; krok(n+1)=3.544829e-001; ng(n+1)=1.061454e+004;
n=73; farx(n+1)=3.083938e+000; foe(n+1)=2.569051e+001; krok(n+1)=8.096016e-007; ng(n+1)=2.094883e+003;
n=74; farx(n+1)=2.798750e+000; foe(n+1)=2.523554e+001; krok(n+1)=3.068907e-001; ng(n+1)=5.504483e+002;
n=75; farx(n+1)=2.800418e+000; foe(n+1)=2.514550e+001; krok(n+1)=3.869575e-001; ng(n+1)=2.995971e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=2.789602e+000; foe(n+1)=2.514276e+001; krok(n+1)=1.570969e-005; ng(n+1)=1.171967e+002;
n=77; farx(n+1)=2.775435e+000; foe(n+1)=2.506088e+001; krok(n+1)=4.503114e-006; ng(n+1)=1.413863e+002;
n=78; farx(n+1)=2.755974e+000; foe(n+1)=2.500954e+001; krok(n+1)=5.807097e-006; ng(n+1)=1.912302e+002;
n=79; farx(n+1)=2.729488e+000; foe(n+1)=2.498179e+001; krok(n+1)=1.284475e-005; ng(n+1)=7.492656e+002;
n=80; farx(n+1)=2.646878e+000; foe(n+1)=2.470810e+001; krok(n+1)=5.216764e-004; ng(n+1)=5.541008e+002;
n=81; farx(n+1)=2.697352e+000; foe(n+1)=2.469292e+001; krok(n+1)=1.308616e-002; ng(n+1)=8.938783e+002;
n=82; farx(n+1)=2.686585e+000; foe(n+1)=2.468316e+001; krok(n+1)=1.281894e-002; ng(n+1)=8.051034e+002;
n=83; farx(n+1)=2.673448e+000; foe(n+1)=2.466990e+001; krok(n+1)=1.199493e-002; ng(n+1)=7.334698e+002;
n=84; farx(n+1)=2.686234e+000; foe(n+1)=2.463512e+001; krok(n+1)=3.281725e-002; ng(n+1)=7.084259e+002;
n=85; farx(n+1)=2.655138e+000; foe(n+1)=2.462344e+001; krok(n+1)=2.967304e-002; ng(n+1)=3.040292e+002;
n=86; farx(n+1)=2.457371e+000; foe(n+1)=2.449591e+001; krok(n+1)=1.402391e+000; ng(n+1)=8.732975e+002;
n=87; farx(n+1)=2.421497e+000; foe(n+1)=2.445282e+001; krok(n+1)=6.400584e-001; ng(n+1)=2.075715e+002;
n=88; farx(n+1)=2.364606e+000; foe(n+1)=2.442756e+001; krok(n+1)=3.346481e-001; ng(n+1)=7.460329e+001;
n=89; farx(n+1)=2.363061e+000; foe(n+1)=2.440535e+001; krok(n+1)=1.315598e+000; ng(n+1)=1.322289e+002;
n=90; farx(n+1)=2.351797e+000; foe(n+1)=2.439797e+001; krok(n+1)=9.073805e-001; ng(n+1)=2.899548e+002;
n=91; farx(n+1)=2.331815e+000; foe(n+1)=2.439218e+001; krok(n+1)=1.549732e+000; ng(n+1)=7.609969e+001;
n=92; farx(n+1)=2.301581e+000; foe(n+1)=2.438688e+001; krok(n+1)=1.914853e+000; ng(n+1)=6.527206e+001;
n=93; farx(n+1)=2.295955e+000; foe(n+1)=2.437492e+001; krok(n+1)=1.905690e+000; ng(n+1)=2.571009e+001;
n=94; farx(n+1)=2.282773e+000; foe(n+1)=2.435590e+001; krok(n+1)=2.054927e+000; ng(n+1)=1.167838e+002;
n=95; farx(n+1)=2.254116e+000; foe(n+1)=2.434657e+001; krok(n+1)=2.067461e+000; ng(n+1)=2.632408e+002;
n=96; farx(n+1)=2.251298e+000; foe(n+1)=2.434172e+001; krok(n+1)=7.927530e-001; ng(n+1)=2.046497e+002;
n=97; farx(n+1)=2.250488e+000; foe(n+1)=2.433529e+001; krok(n+1)=2.674687e-002; ng(n+1)=9.009335e+001;
n=98; farx(n+1)=2.235267e+000; foe(n+1)=2.432455e+001; krok(n+1)=1.103953e+000; ng(n+1)=1.195374e+002;
n=99; farx(n+1)=2.233761e+000; foe(n+1)=2.432321e+001; krok(n+1)=9.961132e-001; ng(n+1)=2.786541e+001;
n=100; farx(n+1)=2.241984e+000; foe(n+1)=2.432263e+001; krok(n+1)=2.190890e+000; ng(n+1)=4.338215e+001;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
