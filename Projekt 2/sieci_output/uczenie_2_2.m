%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.786317e+003; foe(n+1)=4.832948e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.573689e+003; foe(n+1)=3.767576e+003; krok(n+1)=3.395966e-004; ng(n+1)=2.719428e+003;
n=2; farx(n+1)=7.715071e+002; foe(n+1)=8.568911e+002; krok(n+1)=2.456077e-003; ng(n+1)=1.545245e+003;
n=3; farx(n+1)=7.526988e+002; foe(n+1)=6.420186e+002; krok(n+1)=1.915864e-004; ng(n+1)=3.061519e+003;
n=4; farx(n+1)=6.700479e+002; foe(n+1)=5.231399e+002; krok(n+1)=1.945210e-003; ng(n+1)=1.714512e+003;
n=5; farx(n+1)=4.938021e+002; foe(n+1)=4.862752e+002; krok(n+1)=3.485222e-003; ng(n+1)=7.014976e+002;
n=6; farx(n+1)=4.010599e+002; foe(n+1)=4.173075e+002; krok(n+1)=4.241720e-003; ng(n+1)=1.220046e+003;
n=7; farx(n+1)=3.895660e+002; foe(n+1)=4.143323e+002; krok(n+1)=1.532691e-003; ng(n+1)=3.944480e+002;
n=8; farx(n+1)=3.882987e+002; foe(n+1)=4.089361e+002; krok(n+1)=4.776032e-003; ng(n+1)=4.679000e+002;
n=9; farx(n+1)=3.812483e+002; foe(n+1)=4.023550e+002; krok(n+1)=3.722051e-003; ng(n+1)=1.719376e+002;
n=10; farx(n+1)=3.444648e+002; foe(n+1)=3.915030e+002; krok(n+1)=7.142563e-002; ng(n+1)=1.769059e+002;
n=11; farx(n+1)=3.707182e+002; foe(n+1)=3.860985e+002; krok(n+1)=7.903310e-002; ng(n+1)=3.007508e+002;
n=12; farx(n+1)=3.515978e+002; foe(n+1)=3.777323e+002; krok(n+1)=2.723382e-001; ng(n+1)=2.122085e+002;
n=13; farx(n+1)=3.357461e+002; foe(n+1)=3.671929e+002; krok(n+1)=1.752988e-001; ng(n+1)=2.716479e+002;
n=14; farx(n+1)=3.035137e+002; foe(n+1)=3.528411e+002; krok(n+1)=1.382959e+000; ng(n+1)=2.150975e+002;
n=15; farx(n+1)=2.708561e+002; foe(n+1)=3.445098e+002; krok(n+1)=6.995615e-002; ng(n+1)=3.078150e+002;
n=16; farx(n+1)=2.161847e+002; foe(n+1)=3.306867e+002; krok(n+1)=1.261373e-001; ng(n+1)=7.182370e+002;
n=17; farx(n+1)=1.430575e+002; foe(n+1)=3.086695e+002; krok(n+1)=1.439513e-001; ng(n+1)=1.262490e+003;
n=18; farx(n+1)=8.038997e+001; foe(n+1)=2.660435e+002; krok(n+1)=5.005537e-002; ng(n+1)=2.682954e+003;
n=19; farx(n+1)=7.981638e+001; foe(n+1)=2.646669e+002; krok(n+1)=4.880194e-003; ng(n+1)=2.994743e+003;
n=20; farx(n+1)=7.281984e+001; foe(n+1)=2.464442e+002; krok(n+1)=8.268560e-002; ng(n+1)=2.993234e+003;
n=21; farx(n+1)=7.409767e+001; foe(n+1)=2.130693e+002; krok(n+1)=3.655097e-001; ng(n+1)=1.702643e+003;
n=22; farx(n+1)=5.641823e+001; foe(n+1)=1.878322e+002; krok(n+1)=2.663672e-001; ng(n+1)=2.398617e+003;
n=23; farx(n+1)=4.086882e+001; foe(n+1)=1.649028e+002; krok(n+1)=1.827549e-001; ng(n+1)=8.773011e+002;
n=24; farx(n+1)=3.535802e+001; foe(n+1)=1.500764e+002; krok(n+1)=2.429329e-001; ng(n+1)=8.550499e+002;
n=25; farx(n+1)=3.189959e+001; foe(n+1)=1.307156e+002; krok(n+1)=9.725754e-001; ng(n+1)=1.972202e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=2.871170e+001; foe(n+1)=1.282365e+002; krok(n+1)=1.486954e-005; ng(n+1)=9.502654e+002;
n=27; farx(n+1)=2.632189e+001; foe(n+1)=1.248541e+002; krok(n+1)=3.498787e-005; ng(n+1)=7.323115e+002;
n=28; farx(n+1)=2.459903e+001; foe(n+1)=1.168727e+002; krok(n+1)=2.263251e-004; ng(n+1)=4.336572e+002;
n=29; farx(n+1)=1.734116e+001; foe(n+1)=1.007389e+002; krok(n+1)=2.681146e-004; ng(n+1)=5.411776e+002;
n=30; farx(n+1)=1.470296e+001; foe(n+1)=9.490422e+001; krok(n+1)=1.842068e-004; ng(n+1)=1.541704e+003;
n=31; farx(n+1)=1.418418e+001; foe(n+1)=8.690473e+001; krok(n+1)=8.351411e-004; ng(n+1)=3.428254e+003;
n=32; farx(n+1)=1.380030e+001; foe(n+1)=8.595023e+001; krok(n+1)=9.247269e-004; ng(n+1)=3.293760e+003;
n=33; farx(n+1)=1.400406e+001; foe(n+1)=8.458527e+001; krok(n+1)=2.567947e-002; ng(n+1)=3.590475e+003;
n=34; farx(n+1)=1.364669e+001; foe(n+1)=8.224466e+001; krok(n+1)=1.250522e-001; ng(n+1)=2.268108e+003;
n=35; farx(n+1)=9.886802e+000; foe(n+1)=7.655386e+001; krok(n+1)=1.525346e-001; ng(n+1)=2.691234e+003;
n=36; farx(n+1)=6.524076e+000; foe(n+1)=6.677923e+001; krok(n+1)=2.990259e-001; ng(n+1)=1.764208e+003;
n=37; farx(n+1)=6.483492e+000; foe(n+1)=6.390891e+001; krok(n+1)=6.774452e-002; ng(n+1)=5.961698e+002;
n=38; farx(n+1)=5.972982e+000; foe(n+1)=5.955517e+001; krok(n+1)=1.775417e-001; ng(n+1)=7.165091e+002;
n=39; farx(n+1)=5.488335e+000; foe(n+1)=5.596864e+001; krok(n+1)=2.002215e-001; ng(n+1)=8.391227e+002;
n=40; farx(n+1)=4.443363e+000; foe(n+1)=4.924432e+001; krok(n+1)=9.906032e-001; ng(n+1)=2.079342e+003;
n=41; farx(n+1)=4.191873e+000; foe(n+1)=4.777346e+001; krok(n+1)=3.785938e-001; ng(n+1)=1.133016e+003;
n=42; farx(n+1)=3.909916e+000; foe(n+1)=4.665642e+001; krok(n+1)=4.228235e-001; ng(n+1)=1.737390e+002;
n=43; farx(n+1)=3.768436e+000; foe(n+1)=4.632714e+001; krok(n+1)=3.106282e-001; ng(n+1)=3.348400e+002;
n=44; farx(n+1)=3.734457e+000; foe(n+1)=4.620670e+001; krok(n+1)=1.288588e-001; ng(n+1)=6.145039e+002;
n=45; farx(n+1)=3.730017e+000; foe(n+1)=4.600702e+001; krok(n+1)=6.480438e-001; ng(n+1)=2.942436e+002;
n=46; farx(n+1)=3.753839e+000; foe(n+1)=4.574424e+001; krok(n+1)=2.608981e+000; ng(n+1)=2.298851e+002;
n=47; farx(n+1)=3.900663e+000; foe(n+1)=4.554203e+001; krok(n+1)=8.316025e-001; ng(n+1)=1.751299e+002;
n=48; farx(n+1)=3.802883e+000; foe(n+1)=4.537361e+001; krok(n+1)=2.137416e+000; ng(n+1)=2.362834e+002;
n=49; farx(n+1)=3.924870e+000; foe(n+1)=4.489934e+001; krok(n+1)=3.629522e+000; ng(n+1)=2.227004e+002;
n=50; farx(n+1)=3.826757e+000; foe(n+1)=4.389265e+001; krok(n+1)=1.707324e+000; ng(n+1)=4.272839e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=3.827753e+000; foe(n+1)=4.382633e+001; krok(n+1)=7.172631e-006; ng(n+1)=2.835693e+002;
n=52; farx(n+1)=3.794399e+000; foe(n+1)=4.378894e+001; krok(n+1)=2.466737e-005; ng(n+1)=1.069573e+002;
n=53; farx(n+1)=3.756656e+000; foe(n+1)=4.371580e+001; krok(n+1)=6.019984e-006; ng(n+1)=2.657792e+002;
n=54; farx(n+1)=3.637114e+000; foe(n+1)=4.345650e+001; krok(n+1)=6.767654e-004; ng(n+1)=5.606485e+001;
n=55; farx(n+1)=3.584115e+000; foe(n+1)=4.336139e+001; krok(n+1)=4.651747e-004; ng(n+1)=5.347034e+001;
n=56; farx(n+1)=3.570052e+000; foe(n+1)=4.307655e+001; krok(n+1)=3.135616e-003; ng(n+1)=6.274719e+001;
n=57; farx(n+1)=2.985997e+000; foe(n+1)=4.199113e+001; krok(n+1)=3.203828e-002; ng(n+1)=7.468637e+001;
n=58; farx(n+1)=3.550810e+000; foe(n+1)=4.145491e+001; krok(n+1)=5.059367e-002; ng(n+1)=3.906035e+002;
n=59; farx(n+1)=3.580133e+000; foe(n+1)=4.125402e+001; krok(n+1)=2.320653e-002; ng(n+1)=6.658969e+002;
n=60; farx(n+1)=3.382619e+000; foe(n+1)=4.103805e+001; krok(n+1)=3.150760e-002; ng(n+1)=4.577452e+002;
n=61; farx(n+1)=3.573344e+000; foe(n+1)=4.091634e+001; krok(n+1)=3.650930e-002; ng(n+1)=8.066135e+002;
n=62; farx(n+1)=3.120238e+000; foe(n+1)=4.012661e+001; krok(n+1)=8.612093e-001; ng(n+1)=5.828103e+002;
n=63; farx(n+1)=2.518957e+000; foe(n+1)=3.961925e+001; krok(n+1)=6.522451e-001; ng(n+1)=1.979554e+002;
n=64; farx(n+1)=2.284520e+000; foe(n+1)=3.912508e+001; krok(n+1)=5.154350e-001; ng(n+1)=2.988735e+002;
n=65; farx(n+1)=2.228357e+000; foe(n+1)=3.881622e+001; krok(n+1)=4.128102e-001; ng(n+1)=4.140131e+002;
n=66; farx(n+1)=2.331317e+000; foe(n+1)=3.810371e+001; krok(n+1)=8.536620e-001; ng(n+1)=3.505491e+002;
n=67; farx(n+1)=2.196601e+000; foe(n+1)=3.724541e+001; krok(n+1)=2.427388e-001; ng(n+1)=8.633169e+002;
n=68; farx(n+1)=1.789504e+000; foe(n+1)=3.590136e+001; krok(n+1)=3.262230e-001; ng(n+1)=5.143773e+002;
n=69; farx(n+1)=1.650819e+000; foe(n+1)=3.497735e+001; krok(n+1)=7.841683e-002; ng(n+1)=2.126477e+003;
n=70; farx(n+1)=1.501066e+000; foe(n+1)=3.425819e+001; krok(n+1)=6.153853e-002; ng(n+1)=1.538721e+003;
n=71; farx(n+1)=1.400015e+000; foe(n+1)=3.388739e+001; krok(n+1)=9.547529e-002; ng(n+1)=2.197859e+003;
n=72; farx(n+1)=1.319710e+000; foe(n+1)=3.328714e+001; krok(n+1)=1.245141e-001; ng(n+1)=2.589458e+003;
n=73; farx(n+1)=1.337178e+000; foe(n+1)=3.300526e+001; krok(n+1)=2.097829e-001; ng(n+1)=9.586638e+002;
n=74; farx(n+1)=1.332859e+000; foe(n+1)=3.276938e+001; krok(n+1)=3.427774e-001; ng(n+1)=2.490614e+002;
n=75; farx(n+1)=1.275428e+000; foe(n+1)=3.245119e+001; krok(n+1)=6.182030e-001; ng(n+1)=7.243944e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.273661e+000; foe(n+1)=3.244020e+001; krok(n+1)=2.618307e-007; ng(n+1)=4.994938e+002;
n=77; farx(n+1)=1.274830e+000; foe(n+1)=3.243361e+001; krok(n+1)=3.430845e-006; ng(n+1)=1.039416e+002;
n=78; farx(n+1)=1.272443e+000; foe(n+1)=3.242552e+001; krok(n+1)=1.527265e-006; ng(n+1)=1.832304e+002;
n=79; farx(n+1)=1.278725e+000; foe(n+1)=3.239811e+001; krok(n+1)=1.043926e-004; ng(n+1)=4.599171e+001;
n=80; farx(n+1)=1.269737e+000; foe(n+1)=3.230521e+001; krok(n+1)=3.968650e-004; ng(n+1)=4.953610e+001;
n=81; farx(n+1)=1.299461e+000; foe(n+1)=3.216841e+001; krok(n+1)=3.378712e-003; ng(n+1)=8.464608e+001;
n=82; farx(n+1)=1.315026e+000; foe(n+1)=3.203348e+001; krok(n+1)=2.032662e-003; ng(n+1)=1.915638e+002;
n=83; farx(n+1)=1.281303e+000; foe(n+1)=3.193973e+001; krok(n+1)=2.883948e-003; ng(n+1)=7.348194e+002;
n=84; farx(n+1)=1.291834e+000; foe(n+1)=3.192974e+001; krok(n+1)=6.053639e-003; ng(n+1)=1.183860e+003;
n=85; farx(n+1)=1.303558e+000; foe(n+1)=3.191456e+001; krok(n+1)=1.427439e-002; ng(n+1)=1.054098e+003;
n=86; farx(n+1)=1.335392e+000; foe(n+1)=3.178300e+001; krok(n+1)=2.134155e-001; ng(n+1)=9.735755e+002;
n=87; farx(n+1)=1.175800e+000; foe(n+1)=3.080182e+001; krok(n+1)=8.005809e-002; ng(n+1)=6.818126e+002;
n=88; farx(n+1)=1.084688e+000; foe(n+1)=3.040520e+001; krok(n+1)=5.015999e-001; ng(n+1)=1.828315e+003;
n=89; farx(n+1)=1.078515e+000; foe(n+1)=3.000411e+001; krok(n+1)=3.903951e-001; ng(n+1)=8.757687e+002;
n=90; farx(n+1)=1.113744e+000; foe(n+1)=2.966877e+001; krok(n+1)=5.271762e-001; ng(n+1)=1.982923e+003;
n=91; farx(n+1)=1.122584e+000; foe(n+1)=2.957617e+001; krok(n+1)=4.747686e-001; ng(n+1)=7.592489e+002;
n=92; farx(n+1)=1.115875e+000; foe(n+1)=2.946254e+001; krok(n+1)=3.856870e-001; ng(n+1)=4.194378e+002;
n=93; farx(n+1)=1.095629e+000; foe(n+1)=2.939861e+001; krok(n+1)=3.403631e-001; ng(n+1)=3.335523e+002;
n=94; farx(n+1)=1.071191e+000; foe(n+1)=2.923245e+001; krok(n+1)=7.679697e-001; ng(n+1)=7.354975e+002;
n=95; farx(n+1)=1.093396e+000; foe(n+1)=2.907073e+001; krok(n+1)=1.092255e+000; ng(n+1)=1.174500e+003;
n=96; farx(n+1)=1.109006e+000; foe(n+1)=2.895310e+001; krok(n+1)=3.238984e-001; ng(n+1)=7.144876e+002;
n=97; farx(n+1)=1.087406e+000; foe(n+1)=2.883078e+001; krok(n+1)=1.199074e+000; ng(n+1)=6.151099e+002;
n=98; farx(n+1)=1.093772e+000; foe(n+1)=2.877100e+001; krok(n+1)=8.295027e-001; ng(n+1)=1.215893e+003;
n=99; farx(n+1)=1.084373e+000; foe(n+1)=2.871004e+001; krok(n+1)=3.080168e-001; ng(n+1)=4.106880e+002;
n=100; farx(n+1)=1.061356e+000; foe(n+1)=2.864563e+001; krok(n+1)=3.427774e-001; ng(n+1)=4.012928e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
