%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.832581e+003; foe(n+1)=4.628626e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.096716e+003; foe(n+1)=3.005791e+003; krok(n+1)=5.113081e-004; ng(n+1)=3.545477e+003;
n=2; farx(n+1)=1.352737e+003; foe(n+1)=9.297171e+002; krok(n+1)=5.829556e-003; ng(n+1)=1.182554e+003;
n=3; farx(n+1)=1.633982e+003; foe(n+1)=7.047689e+002; krok(n+1)=1.347782e-004; ng(n+1)=5.581636e+003;
n=4; farx(n+1)=2.181082e+003; foe(n+1)=5.980750e+002; krok(n+1)=2.285820e-004; ng(n+1)=5.988922e+003;
n=5; farx(n+1)=1.704391e+003; foe(n+1)=5.211384e+002; krok(n+1)=3.254371e-003; ng(n+1)=9.448089e+002;
n=6; farx(n+1)=4.725259e+002; foe(n+1)=3.840377e+002; krok(n+1)=1.189189e-002; ng(n+1)=9.044452e+002;
n=7; farx(n+1)=3.195492e+002; foe(n+1)=3.673491e+002; krok(n+1)=3.355212e-004; ng(n+1)=2.717753e+003;
n=8; farx(n+1)=2.511883e+002; foe(n+1)=3.579714e+002; krok(n+1)=5.104038e-004; ng(n+1)=2.218881e+003;
n=9; farx(n+1)=1.892059e+002; foe(n+1)=3.395239e+002; krok(n+1)=2.213894e-003; ng(n+1)=2.930456e+003;
n=10; farx(n+1)=1.590122e+002; foe(n+1)=3.324352e+002; krok(n+1)=3.179535e-004; ng(n+1)=5.418608e+003;
n=11; farx(n+1)=1.265091e+002; foe(n+1)=3.055642e+002; krok(n+1)=8.855578e-003; ng(n+1)=6.139495e+003;
n=12; farx(n+1)=1.248062e+002; foe(n+1)=3.044412e+002; krok(n+1)=6.619770e-005; ng(n+1)=2.761583e+003;
n=13; farx(n+1)=1.234307e+002; foe(n+1)=2.834439e+002; krok(n+1)=1.251384e-002; ng(n+1)=3.264423e+003;
n=14; farx(n+1)=1.330121e+002; foe(n+1)=2.747401e+002; krok(n+1)=6.727176e-004; ng(n+1)=2.654761e+003;
n=15; farx(n+1)=1.207918e+002; foe(n+1)=2.467650e+002; krok(n+1)=2.031930e-003; ng(n+1)=3.925756e+003;
n=16; farx(n+1)=1.207310e+002; foe(n+1)=2.429659e+002; krok(n+1)=8.170121e-005; ng(n+1)=7.077025e+003;
n=17; farx(n+1)=1.237228e+002; foe(n+1)=2.319061e+002; krok(n+1)=1.915637e-003; ng(n+1)=4.477280e+003;
n=18; farx(n+1)=1.258253e+002; foe(n+1)=2.259752e+002; krok(n+1)=1.593673e-003; ng(n+1)=4.428400e+003;
n=19; farx(n+1)=1.291942e+002; foe(n+1)=2.203533e+002; krok(n+1)=5.440190e-004; ng(n+1)=2.404838e+003;
n=20; farx(n+1)=1.376536e+002; foe(n+1)=2.051545e+002; krok(n+1)=1.250908e-003; ng(n+1)=4.833901e+003;
n=21; farx(n+1)=1.400399e+002; foe(n+1)=2.022633e+002; krok(n+1)=6.037139e-004; ng(n+1)=7.775248e+003;
n=22; farx(n+1)=1.508274e+002; foe(n+1)=1.930886e+002; krok(n+1)=3.271540e-003; ng(n+1)=7.162908e+003;
n=23; farx(n+1)=1.447467e+002; foe(n+1)=1.788600e+002; krok(n+1)=4.702595e-003; ng(n+1)=5.195420e+003;
n=24; farx(n+1)=1.474787e+002; foe(n+1)=1.722267e+002; krok(n+1)=2.162140e-003; ng(n+1)=5.131611e+003;
n=25; farx(n+1)=1.335036e+002; foe(n+1)=1.620988e+002; krok(n+1)=3.237896e-003; ng(n+1)=5.293846e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=1.162318e+002; foe(n+1)=1.532641e+002; krok(n+1)=5.126555e-006; ng(n+1)=5.610880e+003;
n=27; farx(n+1)=1.129579e+002; foe(n+1)=1.416719e+002; krok(n+1)=1.153689e-005; ng(n+1)=3.634235e+003;
n=28; farx(n+1)=1.089415e+002; foe(n+1)=1.354966e+002; krok(n+1)=1.319295e-004; ng(n+1)=1.163929e+003;
n=29; farx(n+1)=1.066473e+002; foe(n+1)=1.339095e+002; krok(n+1)=2.769106e-005; ng(n+1)=1.285408e+003;
n=30; farx(n+1)=8.890828e+001; foe(n+1)=1.153620e+002; krok(n+1)=1.600738e-003; ng(n+1)=1.382876e+003;
n=31; farx(n+1)=3.184728e+001; foe(n+1)=5.630647e+001; krok(n+1)=6.518159e-004; ng(n+1)=6.542964e+003;
n=32; farx(n+1)=2.512402e+001; foe(n+1)=4.386582e+001; krok(n+1)=9.053002e-004; ng(n+1)=5.207818e+003;
n=33; farx(n+1)=2.223618e+001; foe(n+1)=3.812822e+001; krok(n+1)=3.081204e-004; ng(n+1)=6.247069e+003;
n=34; farx(n+1)=1.414840e+001; foe(n+1)=3.148611e+001; krok(n+1)=8.031637e-004; ng(n+1)=2.893545e+003;
n=35; farx(n+1)=1.327005e+001; foe(n+1)=2.862517e+001; krok(n+1)=1.532691e-003; ng(n+1)=4.268420e+003;
n=36; farx(n+1)=7.815796e+000; foe(n+1)=2.128652e+001; krok(n+1)=5.915258e-003; ng(n+1)=4.624013e+003;
n=37; farx(n+1)=5.827745e+000; foe(n+1)=1.878972e+001; krok(n+1)=4.403244e-003; ng(n+1)=1.977683e+003;
n=38; farx(n+1)=5.403850e+000; foe(n+1)=1.691393e+001; krok(n+1)=1.308315e-003; ng(n+1)=4.425313e+003;
n=39; farx(n+1)=5.142260e+000; foe(n+1)=1.607472e+001; krok(n+1)=2.641497e-003; ng(n+1)=5.601053e+002;
n=40; farx(n+1)=5.164612e+000; foe(n+1)=1.513571e+001; krok(n+1)=5.414124e-003; ng(n+1)=1.301702e+003;
n=41; farx(n+1)=4.684205e+000; foe(n+1)=1.355179e+001; krok(n+1)=3.303815e-003; ng(n+1)=1.645855e+003;
n=42; farx(n+1)=4.436562e+000; foe(n+1)=1.282253e+001; krok(n+1)=3.783371e-003; ng(n+1)=2.870267e+003;
n=43; farx(n+1)=3.989729e+000; foe(n+1)=1.173188e+001; krok(n+1)=1.226153e-002; ng(n+1)=9.201379e+002;
n=44; farx(n+1)=3.479815e+000; foe(n+1)=1.009518e+001; krok(n+1)=1.082825e-002; ng(n+1)=1.643171e+003;
n=45; farx(n+1)=3.394296e+000; foe(n+1)=9.820001e+000; krok(n+1)=5.858547e-003; ng(n+1)=2.548221e+003;
n=46; farx(n+1)=3.264441e+000; foe(n+1)=9.223189e+000; krok(n+1)=8.518466e-003; ng(n+1)=9.462916e+002;
n=47; farx(n+1)=3.198106e+000; foe(n+1)=9.065450e+000; krok(n+1)=1.575380e-002; ng(n+1)=6.007947e+002;
n=48; farx(n+1)=3.194658e+000; foe(n+1)=8.907197e+000; krok(n+1)=1.226153e-002; ng(n+1)=7.298964e+002;
n=49; farx(n+1)=3.181310e+000; foe(n+1)=8.276233e+000; krok(n+1)=3.369650e-002; ng(n+1)=6.252555e+002;
n=50; farx(n+1)=3.193234e+000; foe(n+1)=8.171287e+000; krok(n+1)=8.705227e-003; ng(n+1)=9.185923e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=3.194557e+000; foe(n+1)=8.084501e+000; krok(n+1)=8.613172e-007; ng(n+1)=1.291793e+003;
n=52; farx(n+1)=3.190571e+000; foe(n+1)=8.060991e+000; krok(n+1)=3.652351e-006; ng(n+1)=3.374219e+002;
n=53; farx(n+1)=3.185990e+000; foe(n+1)=8.015643e+000; krok(n+1)=2.394830e-005; ng(n+1)=2.002761e+002;
n=54; farx(n+1)=3.223450e+000; foe(n+1)=7.555564e+000; krok(n+1)=2.954923e-004; ng(n+1)=1.763900e+002;
n=55; farx(n+1)=3.262631e+000; foe(n+1)=7.327772e+000; krok(n+1)=6.483498e-004; ng(n+1)=1.740036e+002;
n=56; farx(n+1)=3.259648e+000; foe(n+1)=7.101313e+000; krok(n+1)=6.647718e-004; ng(n+1)=1.643125e+002;
n=57; farx(n+1)=3.265786e+000; foe(n+1)=6.259305e+000; krok(n+1)=1.012183e-002; ng(n+1)=4.459658e+002;
n=58; farx(n+1)=3.133133e+000; foe(n+1)=5.982339e+000; krok(n+1)=3.783371e-003; ng(n+1)=1.053865e+003;
n=59; farx(n+1)=2.796768e+000; foe(n+1)=5.351470e+000; krok(n+1)=1.019133e-002; ng(n+1)=6.675546e+002;
n=60; farx(n+1)=2.374388e+000; foe(n+1)=4.833940e+000; krok(n+1)=3.452691e-003; ng(n+1)=1.512690e+003;
n=61; farx(n+1)=2.239350e+000; foe(n+1)=4.579164e+000; krok(n+1)=1.636169e-003; ng(n+1)=5.549415e+002;
n=62; farx(n+1)=2.106440e+000; foe(n+1)=4.198799e+000; krok(n+1)=6.905382e-003; ng(n+1)=7.476999e+002;
n=63; farx(n+1)=2.040705e+000; foe(n+1)=4.092623e+000; krok(n+1)=3.343359e-003; ng(n+1)=9.950869e+002;
n=64; farx(n+1)=1.966579e+000; foe(n+1)=3.917540e+000; krok(n+1)=1.071441e-002; ng(n+1)=2.154155e+002;
n=65; farx(n+1)=1.721787e+000; foe(n+1)=3.584925e+000; krok(n+1)=1.633292e-002; ng(n+1)=4.713492e+002;
n=66; farx(n+1)=1.608662e+000; foe(n+1)=3.425711e+000; krok(n+1)=9.348086e-003; ng(n+1)=7.547970e+002;
n=67; farx(n+1)=1.559908e+000; foe(n+1)=3.300932e+000; krok(n+1)=3.350038e-003; ng(n+1)=1.789077e+003;
n=68; farx(n+1)=1.505161e+000; foe(n+1)=3.228231e+000; krok(n+1)=9.808332e-003; ng(n+1)=8.575805e+002;
n=69; farx(n+1)=1.587744e+000; foe(n+1)=3.097851e+000; krok(n+1)=1.479856e-002; ng(n+1)=5.606423e+002;
n=70; farx(n+1)=1.507991e+000; foe(n+1)=2.980370e+000; krok(n+1)=2.399905e-002; ng(n+1)=3.893407e+002;
n=71; farx(n+1)=1.501089e+000; foe(n+1)=2.912028e+000; krok(n+1)=1.865725e-002; ng(n+1)=6.071754e+002;
n=72; farx(n+1)=1.480607e+000; foe(n+1)=2.837379e+000; krok(n+1)=1.574999e-002; ng(n+1)=5.324922e+002;
n=73; farx(n+1)=1.430161e+000; foe(n+1)=2.652140e+000; krok(n+1)=1.269968e-002; ng(n+1)=5.729745e+002;
n=74; farx(n+1)=1.395961e+000; foe(n+1)=2.472598e+000; krok(n+1)=1.434453e-002; ng(n+1)=1.082851e+003;
n=75; farx(n+1)=1.343184e+000; foe(n+1)=2.331421e+000; krok(n+1)=5.213103e-002; ng(n+1)=3.654836e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.340190e+000; foe(n+1)=2.305724e+000; krok(n+1)=6.432245e-007; ng(n+1)=7.800769e+002;
n=77; farx(n+1)=1.340299e+000; foe(n+1)=2.300850e+000; krok(n+1)=4.915571e-007; ng(n+1)=3.163091e+002;
n=78; farx(n+1)=1.340859e+000; foe(n+1)=2.295459e+000; krok(n+1)=9.392863e-006; ng(n+1)=9.621382e+001;
n=79; farx(n+1)=1.307577e+000; foe(n+1)=2.135417e+000; krok(n+1)=3.831728e-004; ng(n+1)=7.951950e+001;
n=80; farx(n+1)=1.307254e+000; foe(n+1)=2.129127e+000; krok(n+1)=2.337468e-004; ng(n+1)=2.726161e+002;
n=81; farx(n+1)=1.315139e+000; foe(n+1)=2.112140e+000; krok(n+1)=1.195607e-003; ng(n+1)=2.825732e+002;
n=82; farx(n+1)=1.319078e+000; foe(n+1)=2.048860e+000; krok(n+1)=1.736904e-003; ng(n+1)=3.140381e+002;
n=83; farx(n+1)=1.306922e+000; foe(n+1)=1.983814e+000; krok(n+1)=2.962391e-003; ng(n+1)=1.965785e+002;
n=84; farx(n+1)=1.369874e+000; foe(n+1)=1.943551e+000; krok(n+1)=7.932207e-003; ng(n+1)=4.533476e+002;
n=85; farx(n+1)=1.343717e+000; foe(n+1)=1.914894e+000; krok(n+1)=7.947193e-003; ng(n+1)=1.826596e+002;
n=86; farx(n+1)=1.298301e+000; foe(n+1)=1.838376e+000; krok(n+1)=1.697325e-002; ng(n+1)=6.454735e+002;
n=87; farx(n+1)=1.269224e+000; foe(n+1)=1.791121e+000; krok(n+1)=9.661282e-003; ng(n+1)=2.829907e+002;
n=88; farx(n+1)=1.237556e+000; foe(n+1)=1.734036e+000; krok(n+1)=1.776017e-002; ng(n+1)=3.774992e+002;
n=89; farx(n+1)=1.224207e+000; foe(n+1)=1.714119e+000; krok(n+1)=3.416397e-003; ng(n+1)=6.902165e+002;
n=90; farx(n+1)=1.232305e+000; foe(n+1)=1.680028e+000; krok(n+1)=2.954351e-002; ng(n+1)=1.043709e+002;
n=91; farx(n+1)=1.247924e+000; foe(n+1)=1.639160e+000; krok(n+1)=4.138006e-002; ng(n+1)=2.471053e+002;
n=92; farx(n+1)=1.271497e+000; foe(n+1)=1.594029e+000; krok(n+1)=2.667694e-002; ng(n+1)=2.581120e+002;
n=93; farx(n+1)=1.272142e+000; foe(n+1)=1.582512e+000; krok(n+1)=1.333847e-002; ng(n+1)=4.922199e+002;
n=94; farx(n+1)=1.271571e+000; foe(n+1)=1.559967e+000; krok(n+1)=2.194150e-002; ng(n+1)=6.247499e+002;
n=95; farx(n+1)=1.280486e+000; foe(n+1)=1.494722e+000; krok(n+1)=6.307385e-002; ng(n+1)=1.039816e+003;
n=96; farx(n+1)=1.298036e+000; foe(n+1)=1.444329e+000; krok(n+1)=3.481722e-002; ng(n+1)=4.694062e+002;
n=97; farx(n+1)=1.297945e+000; foe(n+1)=1.425896e+000; krok(n+1)=3.871130e-002; ng(n+1)=4.790630e+002;
n=98; farx(n+1)=1.291145e+000; foe(n+1)=1.404198e+000; krok(n+1)=6.601690e-002; ng(n+1)=2.520842e+002;
n=99; farx(n+1)=1.280691e+000; foe(n+1)=1.368075e+000; krok(n+1)=1.079418e-001; ng(n+1)=4.408801e+002;
n=100; farx(n+1)=1.274675e+000; foe(n+1)=1.350241e+000; krok(n+1)=8.973423e-002; ng(n+1)=7.491636e+001;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
