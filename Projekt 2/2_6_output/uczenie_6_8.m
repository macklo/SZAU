%uczenie predyktora arx
clear all;
n=0; farx(n+1)=4.312535e+003; foe(n+1)=4.265578e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.170761e+003; foe(n+1)=3.103171e+003; krok(n+1)=4.168272e-004; ng(n+1)=2.756117e+003;
n=2; farx(n+1)=2.859799e+003; foe(n+1)=3.075260e+003; krok(n+1)=3.325065e-004; ng(n+1)=2.189792e+003;
n=3; farx(n+1)=1.812442e+003; foe(n+1)=7.288121e+003; krok(n+1)=1.078226e-003; ng(n+1)=1.885904e+003;
n=4; farx(n+1)=1.459222e+003; foe(n+1)=1.047200e+004; krok(n+1)=6.856719e-004; ng(n+1)=3.052866e+003;
n=5; farx(n+1)=1.399529e+003; foe(n+1)=1.110394e+004; krok(n+1)=1.688698e-002; ng(n+1)=2.400428e+003;
n=6; farx(n+1)=2.924860e+002; foe(n+1)=5.886417e+003; krok(n+1)=1.670810e+000; ng(n+1)=2.229702e+003;
n=7; farx(n+1)=2.287370e+002; foe(n+1)=5.745493e+003; krok(n+1)=1.071175e+000; ng(n+1)=5.313575e+002;
n=8; farx(n+1)=1.471639e+002; foe(n+1)=6.794830e+003; krok(n+1)=5.690980e+000; ng(n+1)=2.641971e+002;
n=9; farx(n+1)=1.019906e+002; foe(n+1)=4.582209e+003; krok(n+1)=1.447720e+000; ng(n+1)=4.396245e+002;
n=10; farx(n+1)=7.487987e+001; foe(n+1)=5.151800e+003; krok(n+1)=8.551938e-001; ng(n+1)=2.118689e+002;
n=11; farx(n+1)=7.361280e+001; foe(n+1)=5.052584e+003; krok(n+1)=1.707324e+000; ng(n+1)=6.862832e+001;
n=12; farx(n+1)=7.083009e+001; foe(n+1)=1.067855e+004; krok(n+1)=7.622761e+000; ng(n+1)=6.279415e+001;
n=13; farx(n+1)=6.654940e+001; foe(n+1)=3.815365e+003; krok(n+1)=3.781569e+000; ng(n+1)=9.529756e+001;
n=14; farx(n+1)=6.163512e+001; foe(n+1)=6.418380e+003; krok(n+1)=5.934102e-001; ng(n+1)=1.187279e+002;
n=15; farx(n+1)=5.488663e+001; foe(n+1)=6.685775e+003; krok(n+1)=1.295594e+000; ng(n+1)=1.236506e+002;
n=16; farx(n+1)=4.747240e+001; foe(n+1)=2.064104e+002; krok(n+1)=4.170483e-001; ng(n+1)=1.095693e+002;
n=17; farx(n+1)=4.119686e+001; foe(n+1)=2.029383e+002; krok(n+1)=7.214488e-001; ng(n+1)=1.291805e+002;
n=18; farx(n+1)=2.988531e+001; foe(n+1)=2.042747e+003; krok(n+1)=1.893841e+000; ng(n+1)=1.589211e+002;
n=19; farx(n+1)=2.675252e+001; foe(n+1)=1.083824e+003; krok(n+1)=2.230612e-001; ng(n+1)=2.139395e+002;
n=20; farx(n+1)=2.573553e+001; foe(n+1)=7.453079e+002; krok(n+1)=2.489869e-001; ng(n+1)=1.078029e+002;
n=21; farx(n+1)=2.399363e+001; foe(n+1)=2.042051e+002; krok(n+1)=7.560816e-001; ng(n+1)=4.430254e+001;
n=22; farx(n+1)=2.316849e+001; foe(n+1)=6.863131e+002; krok(n+1)=1.386016e+000; ng(n+1)=4.208935e+001;
n=23; farx(n+1)=2.294465e+001; foe(n+1)=6.025660e+002; krok(n+1)=1.006124e+000; ng(n+1)=4.526563e+001;
n=24; farx(n+1)=2.287213e+001; foe(n+1)=6.467548e+002; krok(n+1)=1.001473e+000; ng(n+1)=1.457446e+001;
n=25; farx(n+1)=2.284341e+001; foe(n+1)=5.343580e+002; krok(n+1)=9.282331e-001; ng(n+1)=1.039829e+001;
%odnowa zmiennej metryki
n=26; farx(n+1)=2.283874e+001; foe(n+1)=4.546710e+002; krok(n+1)=2.594055e-004; ng(n+1)=1.023560e+001;
n=27; farx(n+1)=2.283570e+001; foe(n+1)=4.670908e+002; krok(n+1)=5.723719e-004; ng(n+1)=4.107403e+000;
n=28; farx(n+1)=2.282544e+001; foe(n+1)=4.210436e+002; krok(n+1)=2.660052e-003; ng(n+1)=3.360326e+000;
n=29; farx(n+1)=2.282453e+001; foe(n+1)=3.883268e+002; krok(n+1)=5.852382e-003; ng(n+1)=9.900027e-001;
n=30; farx(n+1)=2.282361e+001; foe(n+1)=3.984175e+002; krok(n+1)=3.126305e-002; ng(n+1)=8.048445e-001;
n=31; farx(n+1)=2.277026e+001; foe(n+1)=2.962379e+002; krok(n+1)=4.250100e+000; ng(n+1)=9.543349e-001;
n=32; farx(n+1)=7.435381e+000; foe(n+1)=4.443213e+002; krok(n+1)=1.136267e+001; ng(n+1)=4.076202e+000;
n=33; farx(n+1)=7.277153e+000; foe(n+1)=4.418675e+002; krok(n+1)=1.154146e-002; ng(n+1)=6.312437e+001;
n=34; farx(n+1)=6.400524e+000; foe(n+1)=1.923376e+002; krok(n+1)=5.467184e-001; ng(n+1)=6.392211e+001;
n=35; farx(n+1)=5.513559e+000; foe(n+1)=1.419414e+002; krok(n+1)=1.822466e+000; ng(n+1)=8.958919e+001;
n=36; farx(n+1)=4.741004e+000; foe(n+1)=1.931451e+002; krok(n+1)=6.710813e-001; ng(n+1)=1.402515e+002;
n=37; farx(n+1)=4.471607e+000; foe(n+1)=2.391335e+002; krok(n+1)=4.686024e-001; ng(n+1)=6.898371e+001;
n=38; farx(n+1)=4.133387e+000; foe(n+1)=2.064335e+002; krok(n+1)=6.063266e-001; ng(n+1)=5.498371e+001;
n=39; farx(n+1)=3.880157e+000; foe(n+1)=1.508308e+002; krok(n+1)=7.461831e-001; ng(n+1)=4.394393e+001;
n=40; farx(n+1)=3.658651e+000; foe(n+1)=1.464649e+002; krok(n+1)=9.112329e-001; ng(n+1)=1.022529e+002;
n=41; farx(n+1)=3.443626e+000; foe(n+1)=1.663294e+002; krok(n+1)=1.025225e+000; ng(n+1)=2.868694e+001;
n=42; farx(n+1)=3.298869e+000; foe(n+1)=1.457189e+002; krok(n+1)=6.930078e-001; ng(n+1)=4.712507e+001;
n=43; farx(n+1)=3.209947e+000; foe(n+1)=1.530042e+002; krok(n+1)=5.939275e-001; ng(n+1)=8.071310e+001;
n=44; farx(n+1)=3.032590e+000; foe(n+1)=1.360280e+002; krok(n+1)=1.386016e+000; ng(n+1)=6.424030e+000;
n=45; farx(n+1)=2.964886e+000; foe(n+1)=1.448566e+002; krok(n+1)=1.447720e+000; ng(n+1)=5.740914e+001;
n=46; farx(n+1)=2.894617e+000; foe(n+1)=1.291920e+002; krok(n+1)=4.411592e-001; ng(n+1)=2.163847e+001;
n=47; farx(n+1)=2.842930e+000; foe(n+1)=1.279679e+002; krok(n+1)=9.870098e-001; ng(n+1)=3.894456e+001;
n=48; farx(n+1)=2.827282e+000; foe(n+1)=1.242001e+002; krok(n+1)=3.614383e-001; ng(n+1)=1.544994e+001;
n=49; farx(n+1)=2.792351e+000; foe(n+1)=1.291574e+002; krok(n+1)=1.104103e+000; ng(n+1)=3.894892e+001;
n=50; farx(n+1)=2.760381e+000; foe(n+1)=1.302367e+002; krok(n+1)=3.862514e-001; ng(n+1)=4.649769e+000;
%odnowa zmiennej metryki
n=51; farx(n+1)=2.752544e+000; foe(n+1)=1.327763e+002; krok(n+1)=2.114892e-005; ng(n+1)=4.312959e+001;
n=52; farx(n+1)=2.750728e+000; foe(n+1)=1.301125e+002; krok(n+1)=4.129136e-004; ng(n+1)=4.728239e+000;
n=53; farx(n+1)=2.749071e+000; foe(n+1)=1.257033e+002; krok(n+1)=2.972899e-004; ng(n+1)=4.456699e+000;
n=54; farx(n+1)=2.739043e+000; foe(n+1)=1.281016e+002; krok(n+1)=5.962927e-002; ng(n+1)=9.781958e-001;
n=55; farx(n+1)=2.724551e+000; foe(n+1)=1.319797e+002; krok(n+1)=1.068981e-001; ng(n+1)=1.374434e+000;
n=56; farx(n+1)=2.714390e+000; foe(n+1)=1.242832e+002; krok(n+1)=1.139041e-001; ng(n+1)=1.160581e+001;
n=57; farx(n+1)=2.710838e+000; foe(n+1)=1.211557e+002; krok(n+1)=1.480267e-001; ng(n+1)=3.272710e+001;
n=58; farx(n+1)=2.694124e+000; foe(n+1)=1.187886e+002; krok(n+1)=6.428475e-001; ng(n+1)=3.626121e+001;
n=59; farx(n+1)=2.672662e+000; foe(n+1)=1.235274e+002; krok(n+1)=1.663205e+000; ng(n+1)=3.359382e+001;
n=60; farx(n+1)=2.664555e+000; foe(n+1)=1.178382e+002; krok(n+1)=8.612093e-001; ng(n+1)=1.173719e+001;
n=61; farx(n+1)=2.656867e+000; foe(n+1)=1.170349e+002; krok(n+1)=4.310227e-001; ng(n+1)=2.746572e+001;
n=62; farx(n+1)=2.647938e+000; foe(n+1)=1.145268e+002; krok(n+1)=1.201404e+000; ng(n+1)=1.561522e+001;
n=63; farx(n+1)=2.644080e+000; foe(n+1)=1.177300e+002; krok(n+1)=6.954939e-001; ng(n+1)=2.310721e+001;
n=64; farx(n+1)=2.639649e+000; foe(n+1)=1.154383e+002; krok(n+1)=1.075645e+000; ng(n+1)=7.865022e+000;
n=65; farx(n+1)=2.635593e+000; foe(n+1)=1.138612e+002; krok(n+1)=1.552709e+000; ng(n+1)=1.735976e+001;
n=66; farx(n+1)=2.633484e+000; foe(n+1)=1.133087e+002; krok(n+1)=1.218939e+000; ng(n+1)=1.637676e+001;
n=67; farx(n+1)=2.632331e+000; foe(n+1)=1.145812e+002; krok(n+1)=2.885219e-001; ng(n+1)=1.030771e+001;
n=68; farx(n+1)=2.631137e+000; foe(n+1)=1.132155e+002; krok(n+1)=6.128735e-001; ng(n+1)=6.521303e+000;
n=69; farx(n+1)=2.630263e+000; foe(n+1)=1.128357e+002; krok(n+1)=1.273593e+000; ng(n+1)=7.941522e+000;
n=70; farx(n+1)=2.629735e+000; foe(n+1)=1.132448e+002; krok(n+1)=8.295027e-001; ng(n+1)=1.101126e+001;
n=71; farx(n+1)=2.629241e+000; foe(n+1)=1.115861e+002; krok(n+1)=1.236406e+000; ng(n+1)=2.528246e+000;
n=72; farx(n+1)=2.628981e+000; foe(n+1)=1.125998e+002; krok(n+1)=1.499788e+000; ng(n+1)=4.081556e+000;
n=73; farx(n+1)=2.628876e+000; foe(n+1)=1.119838e+002; krok(n+1)=1.420334e+000; ng(n+1)=2.468118e+000;
n=74; farx(n+1)=2.628843e+000; foe(n+1)=1.116641e+002; krok(n+1)=5.988308e-001; ng(n+1)=1.213734e+000;
n=75; farx(n+1)=2.628829e+000; foe(n+1)=1.119116e+002; krok(n+1)=7.310195e-001; ng(n+1)=2.028075e+000;
%odnowa zmiennej metryki
n=76; farx(n+1)=2.628829e+000; foe(n+1)=1.119165e+002; krok(n+1)=1.259177e-005; ng(n+1)=1.331397e-001;
n=77; farx(n+1)=2.628828e+000; foe(n+1)=1.118649e+002; krok(n+1)=1.176191e-003; ng(n+1)=4.157800e-002;
n=78; farx(n+1)=2.628825e+000; foe(n+1)=1.116734e+002; krok(n+1)=6.023036e-004; ng(n+1)=1.246844e-001;
n=79; farx(n+1)=2.628808e+000; foe(n+1)=1.118794e+002; krok(n+1)=3.258189e-003; ng(n+1)=1.541418e-001;
n=80; farx(n+1)=2.628787e+000; foe(n+1)=1.118287e+002; krok(n+1)=4.226396e-002; ng(n+1)=4.317868e-002;
n=81; farx(n+1)=2.628782e+000; foe(n+1)=1.117657e+002; krok(n+1)=5.411858e-002; ng(n+1)=2.036511e-002;
n=82; farx(n+1)=2.628781e+000; foe(n+1)=1.117471e+002; krok(n+1)=7.166328e+000; ng(n+1)=1.700869e-002;
n=83; farx(n+1)=2.628781e+000; foe(n+1)=1.117337e+002; krok(n+1)=1.157348e+000; ng(n+1)=8.544193e-002;
n=84; farx(n+1)=2.628781e+000; foe(n+1)=1.117337e+002; krok(n+1)=3.844490e-006; ng(n+1)=1.504602e-003;
n=85; farx(n+1)=2.628781e+000; foe(n+1)=1.117337e+002; krok(n+1)=2.921410e-005; ng(n+1)=1.504596e-003;
n=86; farx(n+1)=2.628781e+000; foe(n+1)=1.117337e+002; krok(n+1)=1.069939e-005; ng(n+1)=1.504552e-003;
n=87; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=9.929818e-001; ng(n+1)=1.504536e-003;
n=88; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=5.866408e-007; ng(n+1)=1.093396e-005;
n=89; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.075770e-005; ng(n+1)=1.093396e-005;
n=90; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=8.140526e-006; ng(n+1)=1.093384e-005;
n=91; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.210359e-007; ng(n+1)=1.093376e-005;
 % z�y kierunek w metodzie zm - odnowa 
n=92; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=2.530219e-011; ng(n+1)=1.093375e-005;
 % z�y kierunek w metodzie zm - odnowa 
n=93; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=2.143402e-006; ng(n+1)=1.093373e-005;
n=94; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=8.600086e-006; ng(n+1)=9.046285e-006;
n=95; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.337424e-006; ng(n+1)=9.156205e-006;
n=96; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=3.012167e-006; ng(n+1)=9.156474e-006;
n=97; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=2.043113e-008; ng(n+1)=9.156439e-006;
n=98; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=4.472653e-011; ng(n+1)=9.156439e-006;
 % z�y kierunek w metodzie zm - odnowa 
n=99; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=9.888364e-006; ng(n+1)=9.156439e-006;
n=100; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.337424e-006; ng(n+1)=3.221133e-006;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora ARX');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)