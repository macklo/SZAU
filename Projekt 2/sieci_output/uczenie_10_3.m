%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.164222e+003; foe(n+1)=4.446029e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.322886e+003; foe(n+1)=2.569300e+003; krok(n+1)=4.395065e-004; ng(n+1)=7.099593e+003;
n=2; farx(n+1)=9.321115e+002; foe(n+1)=9.830706e+002; krok(n+1)=1.233983e-003; ng(n+1)=4.516340e+003;
n=3; farx(n+1)=1.148335e+003; foe(n+1)=5.474490e+002; krok(n+1)=2.376454e-004; ng(n+1)=9.071841e+003;
n=4; farx(n+1)=1.288275e+003; foe(n+1)=5.303592e+002; krok(n+1)=5.608324e-003; ng(n+1)=1.664435e+003;
n=5; farx(n+1)=8.987538e+002; foe(n+1)=4.633463e+002; krok(n+1)=3.456371e-003; ng(n+1)=5.241889e+002;
n=6; farx(n+1)=8.304986e+002; foe(n+1)=4.578831e+002; krok(n+1)=3.569431e-004; ng(n+1)=7.529463e+002;
n=7; farx(n+1)=5.896745e+002; foe(n+1)=4.150742e+002; krok(n+1)=1.063821e-003; ng(n+1)=3.219997e+003;
n=8; farx(n+1)=3.514287e+002; foe(n+1)=3.834146e+002; krok(n+1)=9.884278e-004; ng(n+1)=2.249214e+003;
n=9; farx(n+1)=2.853309e+002; foe(n+1)=3.569854e+002; krok(n+1)=4.460746e-004; ng(n+1)=3.588988e+003;
n=10; farx(n+1)=2.553047e+002; foe(n+1)=3.443944e+002; krok(n+1)=7.123994e-004; ng(n+1)=5.379155e+003;
n=11; farx(n+1)=2.140087e+002; foe(n+1)=3.344092e+002; krok(n+1)=2.456077e-003; ng(n+1)=5.637114e+003;
n=12; farx(n+1)=1.611122e+002; foe(n+1)=3.245920e+002; krok(n+1)=8.989513e-004; ng(n+1)=5.345143e+003;
n=13; farx(n+1)=1.216703e+002; foe(n+1)=3.152311e+002; krok(n+1)=1.869974e-003; ng(n+1)=4.270670e+003;
n=14; farx(n+1)=1.105815e+002; foe(n+1)=2.974049e+002; krok(n+1)=2.513970e-003; ng(n+1)=3.091398e+003;
n=15; farx(n+1)=1.112582e+002; foe(n+1)=2.931056e+002; krok(n+1)=2.808967e-004; ng(n+1)=4.482430e+003;
n=16; farx(n+1)=1.038586e+002; foe(n+1)=2.699413e+002; krok(n+1)=8.624636e-003; ng(n+1)=4.220565e+003;
n=17; farx(n+1)=9.214197e+001; foe(n+1)=2.575948e+002; krok(n+1)=1.389400e-003; ng(n+1)=2.765547e+003;
n=18; farx(n+1)=7.798524e+001; foe(n+1)=2.438142e+002; krok(n+1)=1.120856e-003; ng(n+1)=1.982620e+003;
n=19; farx(n+1)=6.829193e+001; foe(n+1)=2.393268e+002; krok(n+1)=6.130864e-004; ng(n+1)=2.836318e+003;
n=20; farx(n+1)=6.608026e+001; foe(n+1)=2.352896e+002; krok(n+1)=9.696269e-004; ng(n+1)=3.996232e+003;
n=21; farx(n+1)=6.891817e+001; foe(n+1)=2.127207e+002; krok(n+1)=1.877073e-003; ng(n+1)=6.901251e+003;
n=22; farx(n+1)=6.892285e+001; foe(n+1)=2.127150e+002; krok(n+1)=6.368120e-008; ng(n+1)=3.965167e+003;
n=23; farx(n+1)=6.909250e+001; foe(n+1)=2.073135e+002; krok(n+1)=1.186919e-003; ng(n+1)=1.656090e+003;
n=24; farx(n+1)=6.931040e+001; foe(n+1)=2.046264e+002; krok(n+1)=4.156331e-005; ng(n+1)=4.811583e+004;
n=25; farx(n+1)=8.043548e+001; foe(n+1)=1.951123e+002; krok(n+1)=4.975647e-003; ng(n+1)=1.396930e+004;
%odnowa zmiennej metryki
n=26; farx(n+1)=8.048434e+001; foe(n+1)=1.936758e+002; krok(n+1)=3.440769e-007; ng(n+1)=5.362709e+003;
n=27; farx(n+1)=6.797268e+001; foe(n+1)=1.876687e+002; krok(n+1)=3.012945e-005; ng(n+1)=2.107596e+003;
n=28; farx(n+1)=4.054314e+001; foe(n+1)=1.596830e+002; krok(n+1)=2.004501e-004; ng(n+1)=5.113557e+003;
n=29; farx(n+1)=2.266782e+001; foe(n+1)=1.366633e+002; krok(n+1)=1.687093e-003; ng(n+1)=6.975157e+003;
n=30; farx(n+1)=2.266595e+001; foe(n+1)=1.366538e+002; krok(n+1)=3.962771e-007; ng(n+1)=3.855233e+003;
n=31; farx(n+1)=1.715223e+001; foe(n+1)=1.275895e+002; krok(n+1)=7.101242e-004; ng(n+1)=4.905454e+003;
n=32; farx(n+1)=1.709507e+001; foe(n+1)=1.269211e+002; krok(n+1)=2.201687e-005; ng(n+1)=3.241835e+004;
n=33; farx(n+1)=1.769220e+001; foe(n+1)=1.168930e+002; krok(n+1)=2.345612e-004; ng(n+1)=2.848300e+004;
n=34; farx(n+1)=2.934313e+001; foe(n+1)=4.598416e+001; krok(n+1)=3.159755e-003; ng(n+1)=1.982014e+004;
n=35; farx(n+1)=3.012981e+001; foe(n+1)=3.649510e+001; krok(n+1)=1.427772e-003; ng(n+1)=9.502321e+002;
n=36; farx(n+1)=2.720377e+001; foe(n+1)=3.420517e+001; krok(n+1)=4.395065e-004; ng(n+1)=1.336563e+003;
n=37; farx(n+1)=2.468995e+001; foe(n+1)=3.250735e+001; krok(n+1)=7.663457e-004; ng(n+1)=1.008088e+003;
n=38; farx(n+1)=2.066715e+001; foe(n+1)=2.593288e+001; krok(n+1)=1.142218e-002; ng(n+1)=4.459767e+002;
n=39; farx(n+1)=2.018395e+001; foe(n+1)=2.389463e+001; krok(n+1)=4.312903e-003; ng(n+1)=1.462896e+003;
n=40; farx(n+1)=1.808318e+001; foe(n+1)=2.177435e+001; krok(n+1)=6.700076e-003; ng(n+1)=4.782251e+002;
n=41; farx(n+1)=1.797082e+001; foe(n+1)=2.141649e+001; krok(n+1)=3.752980e-003; ng(n+1)=6.591133e+002;
n=42; farx(n+1)=1.680440e+001; foe(n+1)=1.985910e+001; krok(n+1)=6.681129e-003; ng(n+1)=1.142906e+003;
n=43; farx(n+1)=1.499885e+001; foe(n+1)=1.876029e+001; krok(n+1)=7.852675e-003; ng(n+1)=2.535768e+002;
n=44; farx(n+1)=1.465701e+001; foe(n+1)=1.766180e+001; krok(n+1)=1.479563e-002; ng(n+1)=5.999305e+002;
n=45; farx(n+1)=1.307298e+001; foe(n+1)=1.609692e+001; krok(n+1)=9.367765e-003; ng(n+1)=5.867725e+002;
n=46; farx(n+1)=1.081057e+001; foe(n+1)=1.500058e+001; krok(n+1)=6.216787e-003; ng(n+1)=4.644750e+002;
n=47; farx(n+1)=9.909740e+000; foe(n+1)=1.397937e+001; krok(n+1)=6.349841e-003; ng(n+1)=4.321948e+002;
n=48; farx(n+1)=8.273696e+000; foe(n+1)=1.232907e+001; krok(n+1)=8.029022e-003; ng(n+1)=8.037957e+002;
n=49; farx(n+1)=7.952689e+000; foe(n+1)=1.170941e+001; krok(n+1)=9.600724e-003; ng(n+1)=4.821669e+002;
n=50; farx(n+1)=7.249319e+000; foe(n+1)=1.124340e+001; krok(n+1)=6.192369e-003; ng(n+1)=3.586300e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=7.094759e+000; foe(n+1)=1.112935e+001; krok(n+1)=4.789661e-005; ng(n+1)=2.385673e+002;
n=52; farx(n+1)=7.054858e+000; foe(n+1)=1.104845e+001; krok(n+1)=1.737277e-005; ng(n+1)=3.117939e+002;
n=53; farx(n+1)=6.258092e+000; foe(n+1)=9.937229e+000; krok(n+1)=1.835515e-003; ng(n+1)=1.459083e+002;
n=54; farx(n+1)=5.955816e+000; foe(n+1)=9.040072e+000; krok(n+1)=1.841613e-003; ng(n+1)=3.338265e+002;
n=55; farx(n+1)=5.807423e+000; foe(n+1)=8.746610e+000; krok(n+1)=7.663457e-004; ng(n+1)=6.755418e+002;
n=56; farx(n+1)=5.999974e+000; foe(n+1)=8.458518e+000; krok(n+1)=2.326282e-004; ng(n+1)=6.445501e+002;
n=57; farx(n+1)=5.447486e+000; foe(n+1)=7.892671e+000; krok(n+1)=3.363630e-003; ng(n+1)=6.215983e+002;
n=58; farx(n+1)=5.126320e+000; foe(n+1)=7.372993e+000; krok(n+1)=8.130648e-003; ng(n+1)=4.756204e+002;
n=59; farx(n+1)=4.390650e+000; foe(n+1)=6.698167e+000; krok(n+1)=1.626130e-002; ng(n+1)=3.978535e+002;
n=60; farx(n+1)=4.078645e+000; foe(n+1)=6.431357e+000; krok(n+1)=1.965086e-003; ng(n+1)=4.781188e+002;
n=61; farx(n+1)=3.475622e+000; foe(n+1)=5.963729e+000; krok(n+1)=4.830641e-003; ng(n+1)=4.230966e+002;
n=62; farx(n+1)=2.867895e+000; foe(n+1)=5.340470e+000; krok(n+1)=6.545122e-003; ng(n+1)=1.827214e+002;
n=63; farx(n+1)=2.555901e+000; foe(n+1)=5.031291e+000; krok(n+1)=1.008597e-002; ng(n+1)=1.741674e+002;
n=64; farx(n+1)=2.470682e+000; foe(n+1)=4.811119e+000; krok(n+1)=4.150488e-003; ng(n+1)=4.662655e+002;
n=65; farx(n+1)=2.289816e+000; foe(n+1)=4.606610e+000; krok(n+1)=1.059021e-002; ng(n+1)=2.512292e+002;
n=66; farx(n+1)=2.038156e+000; foe(n+1)=4.193595e+000; krok(n+1)=2.334578e-002; ng(n+1)=2.544522e+002;
n=67; farx(n+1)=1.919672e+000; foe(n+1)=3.853040e+000; krok(n+1)=1.806042e-002; ng(n+1)=1.220354e+002;
n=68; farx(n+1)=1.769783e+000; foe(n+1)=3.488413e+000; krok(n+1)=1.532189e-002; ng(n+1)=1.470759e+002;
n=69; farx(n+1)=1.723598e+000; foe(n+1)=3.167820e+000; krok(n+1)=1.019447e-002; ng(n+1)=3.749039e+002;
n=70; farx(n+1)=1.714184e+000; foe(n+1)=3.128986e+000; krok(n+1)=3.340564e-003; ng(n+1)=2.453060e+002;
n=71; farx(n+1)=1.604269e+000; foe(n+1)=2.950579e+000; krok(n+1)=1.958773e-002; ng(n+1)=7.032546e+001;
n=72; farx(n+1)=1.560963e+000; foe(n+1)=2.880435e+000; krok(n+1)=5.478088e-003; ng(n+1)=2.361495e+002;
n=73; farx(n+1)=1.454210e+000; foe(n+1)=2.578801e+000; krok(n+1)=3.072094e-002; ng(n+1)=9.791490e+001;
n=74; farx(n+1)=1.412554e+000; foe(n+1)=2.364503e+000; krok(n+1)=2.642647e-002; ng(n+1)=2.886647e+002;
n=75; farx(n+1)=1.411547e+000; foe(n+1)=2.301490e+000; krok(n+1)=7.632150e-003; ng(n+1)=3.161170e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.407181e+000; foe(n+1)=2.275999e+000; krok(n+1)=1.404624e-005; ng(n+1)=2.049404e+002;
n=77; farx(n+1)=1.404442e+000; foe(n+1)=2.267446e+000; krok(n+1)=8.333245e-006; ng(n+1)=1.359716e+002;
n=78; farx(n+1)=1.397582e+000; foe(n+1)=2.250424e+000; krok(n+1)=4.687315e-005; ng(n+1)=8.546244e+001;
n=79; farx(n+1)=1.386508e+000; foe(n+1)=2.156704e+000; krok(n+1)=7.625735e-004; ng(n+1)=5.583637e+001;
n=80; farx(n+1)=1.386150e+000; foe(n+1)=2.091466e+000; krok(n+1)=5.067687e-004; ng(n+1)=5.524993e+001;
n=81; farx(n+1)=1.374625e+000; foe(n+1)=2.007999e+000; krok(n+1)=1.897913e-003; ng(n+1)=5.376621e+001;
n=82; farx(n+1)=1.381072e+000; foe(n+1)=1.971624e+000; krok(n+1)=1.427772e-003; ng(n+1)=8.526981e+001;
n=83; farx(n+1)=1.381415e+000; foe(n+1)=1.921643e+000; krok(n+1)=9.533414e-003; ng(n+1)=1.034302e+002;
n=84; farx(n+1)=1.366833e+000; foe(n+1)=1.875260e+000; krok(n+1)=1.932256e-002; ng(n+1)=1.633277e+002;
n=85; farx(n+1)=1.340102e+000; foe(n+1)=1.815298e+000; krok(n+1)=1.301252e-002; ng(n+1)=2.301674e+002;
n=86; farx(n+1)=1.290377e+000; foe(n+1)=1.733714e+000; krok(n+1)=1.477175e-002; ng(n+1)=3.742610e+002;
n=87; farx(n+1)=1.238236e+000; foe(n+1)=1.670724e+000; krok(n+1)=1.390246e-002; ng(n+1)=1.060813e+002;
n=88; farx(n+1)=1.189969e+000; foe(n+1)=1.615576e+000; krok(n+1)=1.183052e-002; ng(n+1)=1.589842e+002;
n=89; farx(n+1)=1.153985e+000; foe(n+1)=1.572794e+000; krok(n+1)=2.439969e-002; ng(n+1)=9.303685e+001;
n=90; farx(n+1)=1.120485e+000; foe(n+1)=1.530071e+000; krok(n+1)=1.394089e-002; ng(n+1)=2.045331e+002;
n=91; farx(n+1)=1.077473e+000; foe(n+1)=1.483665e+000; krok(n+1)=2.243330e-002; ng(n+1)=1.495449e+002;
n=92; farx(n+1)=1.062478e+000; foe(n+1)=1.459791e+000; krok(n+1)=1.778972e-002; ng(n+1)=6.049508e+001;
n=93; farx(n+1)=1.045485e+000; foe(n+1)=1.427794e+000; krok(n+1)=1.778267e-002; ng(n+1)=1.181372e+002;
n=94; farx(n+1)=1.019656e+000; foe(n+1)=1.369129e+000; krok(n+1)=1.729712e-002; ng(n+1)=1.940370e+002;
n=95; farx(n+1)=1.006978e+000; foe(n+1)=1.348485e+000; krok(n+1)=1.865271e-002; ng(n+1)=1.145362e+002;
n=96; farx(n+1)=9.918202e-001; foe(n+1)=1.311513e+000; krok(n+1)=1.856024e-002; ng(n+1)=6.425518e+001;
n=97; farx(n+1)=9.923944e-001; foe(n+1)=1.278559e+000; krok(n+1)=4.256083e-002; ng(n+1)=1.615131e+002;
n=98; farx(n+1)=9.987003e-001; foe(n+1)=1.191057e+000; krok(n+1)=1.204927e-001; ng(n+1)=2.142753e+002;
n=99; farx(n+1)=9.978894e-001; foe(n+1)=1.155012e+000; krok(n+1)=1.448480e-002; ng(n+1)=2.627868e+002;
n=100; farx(n+1)=9.959712e-001; foe(n+1)=1.115824e+000; krok(n+1)=5.625683e-002; ng(n+1)=2.936232e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
