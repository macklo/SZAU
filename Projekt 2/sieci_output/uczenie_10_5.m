%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.567760e+003; foe(n+1)=4.737427e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.418975e+003; foe(n+1)=3.583345e+003; krok(n+1)=3.467612e-004; ng(n+1)=5.903875e+003;
n=2; farx(n+1)=9.535570e+002; foe(n+1)=1.191850e+003; krok(n+1)=9.722412e-004; ng(n+1)=4.947510e+003;
n=3; farx(n+1)=1.076151e+003; foe(n+1)=8.963003e+002; krok(n+1)=1.185248e-004; ng(n+1)=1.003279e+004;
n=4; farx(n+1)=1.416430e+003; foe(n+1)=8.459432e+002; krok(n+1)=4.898436e-004; ng(n+1)=6.407198e+003;
n=5; farx(n+1)=5.720142e+002; foe(n+1)=5.590328e+002; krok(n+1)=1.244697e-002; ng(n+1)=1.352489e+003;
n=6; farx(n+1)=5.447245e+002; foe(n+1)=5.536930e+002; krok(n+1)=1.106017e-005; ng(n+1)=4.598917e+003;
n=7; farx(n+1)=1.449349e+002; foe(n+1)=3.083329e+002; krok(n+1)=6.518015e-004; ng(n+1)=5.226299e+003;
n=8; farx(n+1)=1.392382e+002; foe(n+1)=3.029891e+002; krok(n+1)=2.418894e-005; ng(n+1)=4.446280e+003;
n=9; farx(n+1)=1.204463e+002; foe(n+1)=2.981740e+002; krok(n+1)=2.902986e-004; ng(n+1)=2.960713e+003;
n=10; farx(n+1)=1.121226e+002; foe(n+1)=2.930758e+002; krok(n+1)=7.138862e-004; ng(n+1)=3.549058e+003;
n=11; farx(n+1)=1.186806e+002; foe(n+1)=2.844084e+002; krok(n+1)=2.293768e-003; ng(n+1)=3.754337e+003;
n=12; farx(n+1)=1.077630e+002; foe(n+1)=2.516115e+002; krok(n+1)=4.440041e-003; ng(n+1)=2.946989e+003;
n=13; farx(n+1)=1.014048e+002; foe(n+1)=2.474476e+002; krok(n+1)=4.129769e-004; ng(n+1)=1.283110e+003;
n=14; farx(n+1)=8.478962e+001; foe(n+1)=2.328936e+002; krok(n+1)=5.637357e-003; ng(n+1)=3.008832e+003;
n=15; farx(n+1)=7.228080e+001; foe(n+1)=2.257098e+002; krok(n+1)=1.250114e-003; ng(n+1)=1.701823e+003;
n=16; farx(n+1)=6.865046e+001; foe(n+1)=2.190597e+002; krok(n+1)=1.389400e-003; ng(n+1)=3.850763e+003;
n=17; farx(n+1)=6.560898e+001; foe(n+1)=2.156968e+002; krok(n+1)=1.070672e-003; ng(n+1)=1.118703e+003;
n=18; farx(n+1)=6.807387e+001; foe(n+1)=2.112185e+002; krok(n+1)=1.812955e-003; ng(n+1)=2.682227e+003;
n=19; farx(n+1)=6.457815e+001; foe(n+1)=2.042801e+002; krok(n+1)=1.505777e-003; ng(n+1)=9.691322e+002;
n=20; farx(n+1)=5.855324e+001; foe(n+1)=2.023022e+002; krok(n+1)=3.910575e-004; ng(n+1)=3.074564e+003;
n=21; farx(n+1)=3.615202e+001; foe(n+1)=1.789068e+002; krok(n+1)=3.625910e-003; ng(n+1)=3.019249e+003;
n=22; farx(n+1)=3.371108e+001; foe(n+1)=1.767103e+002; krok(n+1)=1.954789e-004; ng(n+1)=2.183148e+003;
n=23; farx(n+1)=3.353100e+001; foe(n+1)=1.755338e+002; krok(n+1)=1.127391e-003; ng(n+1)=3.253230e+003;
n=24; farx(n+1)=3.284412e+001; foe(n+1)=1.711108e+002; krok(n+1)=1.150834e-003; ng(n+1)=2.212703e+003;
n=25; farx(n+1)=2.773881e+001; foe(n+1)=1.591751e+002; krok(n+1)=7.740461e-004; ng(n+1)=5.197139e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=2.795568e+001; foe(n+1)=1.579560e+002; krok(n+1)=1.130174e-006; ng(n+1)=5.283583e+003;
n=27; farx(n+1)=3.216759e+001; foe(n+1)=1.475820e+002; krok(n+1)=2.286087e-005; ng(n+1)=3.322987e+003;
n=28; farx(n+1)=3.405767e+001; foe(n+1)=1.415618e+002; krok(n+1)=2.857275e-005; ng(n+1)=2.663613e+003;
n=29; farx(n+1)=2.441975e+001; foe(n+1)=1.312867e+002; krok(n+1)=3.375449e-004; ng(n+1)=1.292959e+003;
n=30; farx(n+1)=1.818967e+001; foe(n+1)=1.148754e+002; krok(n+1)=1.172806e-004; ng(n+1)=4.492486e+003;
n=31; farx(n+1)=1.284508e+001; foe(n+1)=8.879025e+001; krok(n+1)=1.059163e-003; ng(n+1)=2.422506e+003;
n=32; farx(n+1)=1.105389e+001; foe(n+1)=7.941142e+001; krok(n+1)=1.491544e-004; ng(n+1)=7.522394e+003;
n=33; farx(n+1)=1.045984e+001; foe(n+1)=6.855605e+001; krok(n+1)=4.532388e-004; ng(n+1)=6.150095e+003;
n=34; farx(n+1)=1.010134e+001; foe(n+1)=4.983723e+001; krok(n+1)=8.336543e-004; ng(n+1)=1.091964e+004;
n=35; farx(n+1)=1.020873e+001; foe(n+1)=4.457707e+001; krok(n+1)=4.168272e-004; ng(n+1)=3.610352e+003;
n=36; farx(n+1)=1.007489e+001; foe(n+1)=4.115360e+001; krok(n+1)=5.945797e-004; ng(n+1)=7.682091e+003;
n=37; farx(n+1)=9.756481e+000; foe(n+1)=3.513345e+001; krok(n+1)=1.736904e-003; ng(n+1)=3.899943e+003;
n=38; farx(n+1)=8.368204e+000; foe(n+1)=2.833651e+001; krok(n+1)=3.516052e-003; ng(n+1)=3.343973e+003;
n=39; farx(n+1)=9.018134e+000; foe(n+1)=2.426534e+001; krok(n+1)=9.914730e-004; ng(n+1)=5.323628e+003;
n=40; farx(n+1)=9.525720e+000; foe(n+1)=2.186528e+001; krok(n+1)=1.387045e-003; ng(n+1)=1.426128e+003;
n=41; farx(n+1)=1.012540e+001; foe(n+1)=1.977961e+001; krok(n+1)=2.586254e-003; ng(n+1)=1.034394e+003;
n=42; farx(n+1)=9.198388e+000; foe(n+1)=1.770159e+001; krok(n+1)=2.651602e-003; ng(n+1)=2.140434e+003;
n=43; farx(n+1)=8.234175e+000; foe(n+1)=1.667862e+001; krok(n+1)=3.127439e-003; ng(n+1)=7.329718e+002;
n=44; farx(n+1)=6.384512e+000; foe(n+1)=1.477798e+001; krok(n+1)=2.548617e-003; ng(n+1)=2.312210e+003;
n=45; farx(n+1)=5.588715e+000; foe(n+1)=1.411807e+001; krok(n+1)=1.629650e-003; ng(n+1)=1.734338e+003;
n=46; farx(n+1)=5.340460e+000; foe(n+1)=1.368596e+001; krok(n+1)=3.056641e-003; ng(n+1)=5.813553e+002;
n=47; farx(n+1)=5.265683e+000; foe(n+1)=1.312838e+001; krok(n+1)=4.777659e-003; ng(n+1)=4.880640e+002;
n=48; farx(n+1)=4.357859e+000; foe(n+1)=1.205707e+001; krok(n+1)=9.888795e-003; ng(n+1)=4.515935e+002;
n=49; farx(n+1)=4.068321e+000; foe(n+1)=1.157413e+001; krok(n+1)=1.742611e-003; ng(n+1)=1.910814e+003;
n=50; farx(n+1)=3.879783e+000; foe(n+1)=1.110533e+001; krok(n+1)=5.400718e-003; ng(n+1)=1.678446e+003;
%odnowa zmiennej metryki
n=51; farx(n+1)=3.871718e+000; foe(n+1)=1.065325e+001; krok(n+1)=3.021959e-006; ng(n+1)=2.339130e+003;
n=52; farx(n+1)=3.493459e+000; foe(n+1)=1.010395e+001; krok(n+1)=2.867992e-005; ng(n+1)=8.185443e+002;
n=53; farx(n+1)=3.411891e+000; foe(n+1)=9.899591e+000; krok(n+1)=1.526516e-005; ng(n+1)=6.322671e+002;
n=54; farx(n+1)=3.258551e+000; foe(n+1)=9.550283e+000; krok(n+1)=5.291906e-005; ng(n+1)=5.190998e+002;
n=55; farx(n+1)=3.004698e+000; foe(n+1)=9.046172e+000; krok(n+1)=1.317301e-003; ng(n+1)=1.543719e+002;
n=56; farx(n+1)=2.225421e+000; foe(n+1)=7.436663e+000; krok(n+1)=3.475168e-003; ng(n+1)=1.530062e+002;
n=57; farx(n+1)=2.056287e+000; foe(n+1)=7.066824e+000; krok(n+1)=9.843744e-004; ng(n+1)=5.296386e+002;
n=58; farx(n+1)=1.808907e+000; foe(n+1)=6.386476e+000; krok(n+1)=1.721185e-003; ng(n+1)=5.708443e+002;
n=59; farx(n+1)=1.581158e+000; foe(n+1)=5.188125e+000; krok(n+1)=3.752980e-003; ng(n+1)=9.244617e+002;
n=60; farx(n+1)=1.529009e+000; foe(n+1)=4.824858e+000; krok(n+1)=5.842554e-004; ng(n+1)=1.782088e+003;
n=61; farx(n+1)=1.530678e+000; foe(n+1)=4.649782e+000; krok(n+1)=3.456065e-003; ng(n+1)=8.800538e+002;
n=62; farx(n+1)=1.496711e+000; foe(n+1)=4.348777e+000; krok(n+1)=6.274942e-003; ng(n+1)=1.199557e+003;
n=63; farx(n+1)=1.458247e+000; foe(n+1)=4.076792e+000; krok(n+1)=2.343845e-003; ng(n+1)=5.785216e+002;
n=64; farx(n+1)=1.437870e+000; foe(n+1)=3.906464e+000; krok(n+1)=3.096691e-003; ng(n+1)=8.474035e+002;
n=65; farx(n+1)=1.412532e+000; foe(n+1)=3.823473e+000; krok(n+1)=2.364412e-003; ng(n+1)=5.553664e+002;
n=66; farx(n+1)=1.310127e+000; foe(n+1)=3.391096e+000; krok(n+1)=2.362755e-002; ng(n+1)=7.577913e+002;
n=67; farx(n+1)=1.258097e+000; foe(n+1)=3.285445e+000; krok(n+1)=4.055547e-003; ng(n+1)=6.852584e+002;
n=68; farx(n+1)=1.194505e+000; foe(n+1)=3.136930e+000; krok(n+1)=4.853566e-003; ng(n+1)=9.379592e+002;
n=69; farx(n+1)=1.085582e+000; foe(n+1)=2.948648e+000; krok(n+1)=1.136199e-002; ng(n+1)=3.314045e+002;
n=70; farx(n+1)=1.044143e+000; foe(n+1)=2.862679e+000; krok(n+1)=3.516052e-003; ng(n+1)=4.530946e+002;
n=71; farx(n+1)=1.045013e+000; foe(n+1)=2.675619e+000; krok(n+1)=1.226153e-002; ng(n+1)=3.660777e+002;
n=72; farx(n+1)=1.067776e+000; foe(n+1)=2.534254e+000; krok(n+1)=6.727259e-003; ng(n+1)=3.649416e+002;
n=73; farx(n+1)=1.066126e+000; foe(n+1)=2.448493e+000; krok(n+1)=5.906887e-003; ng(n+1)=9.079634e+002;
n=74; farx(n+1)=1.006431e+000; foe(n+1)=2.346562e+000; krok(n+1)=1.865271e-002; ng(n+1)=2.600371e+002;
n=75; farx(n+1)=9.854466e-001; foe(n+1)=2.292839e+000; krok(n+1)=7.557177e-003; ng(n+1)=6.701257e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=9.787625e-001; foe(n+1)=2.226906e+000; krok(n+1)=1.508897e-006; ng(n+1)=1.111664e+003;
n=77; farx(n+1)=9.775671e-001; foe(n+1)=2.222011e+000; krok(n+1)=6.394006e-006; ng(n+1)=1.514781e+002;
n=78; farx(n+1)=9.791388e-001; foe(n+1)=2.211914e+000; krok(n+1)=2.354047e-005; ng(n+1)=1.111512e+002;
n=79; farx(n+1)=9.759870e-001; foe(n+1)=2.205744e+000; krok(n+1)=4.176829e-005; ng(n+1)=7.187370e+001;
n=80; farx(n+1)=9.698109e-001; foe(n+1)=2.170285e+000; krok(n+1)=3.748417e-004; ng(n+1)=5.772751e+001;
n=81; farx(n+1)=9.759383e-001; foe(n+1)=2.127292e+000; krok(n+1)=7.068748e-004; ng(n+1)=4.407723e+001;
n=82; farx(n+1)=9.893347e-001; foe(n+1)=2.020375e+000; krok(n+1)=1.573435e-003; ng(n+1)=4.921110e+001;
n=83; farx(n+1)=9.887208e-001; foe(n+1)=1.920719e+000; krok(n+1)=1.632192e-003; ng(n+1)=1.198525e+002;
n=84; farx(n+1)=9.834794e-001; foe(n+1)=1.830663e+000; krok(n+1)=3.576264e-003; ng(n+1)=4.091960e+002;
n=85; farx(n+1)=9.538622e-001; foe(n+1)=1.761587e+000; krok(n+1)=3.941672e-003; ng(n+1)=4.740075e+002;
n=86; farx(n+1)=9.340211e-001; foe(n+1)=1.689273e+000; krok(n+1)=4.738412e-003; ng(n+1)=8.274121e+002;
n=87; farx(n+1)=9.295191e-001; foe(n+1)=1.666389e+000; krok(n+1)=2.778799e-003; ng(n+1)=3.657162e+002;
n=88; farx(n+1)=8.509107e-001; foe(n+1)=1.584808e+000; krok(n+1)=1.877194e-002; ng(n+1)=4.711349e+002;
n=89; farx(n+1)=7.981461e-001; foe(n+1)=1.520189e+000; krok(n+1)=6.283621e-003; ng(n+1)=7.547815e+002;
n=90; farx(n+1)=7.632404e-001; foe(n+1)=1.421876e+000; krok(n+1)=1.375044e-002; ng(n+1)=2.677004e+002;
n=91; farx(n+1)=7.413329e-001; foe(n+1)=1.393664e+000; krok(n+1)=4.483426e-003; ng(n+1)=2.321958e+002;
n=92; farx(n+1)=6.478495e-001; foe(n+1)=1.257310e+000; krok(n+1)=1.080437e-002; ng(n+1)=3.791943e+002;
n=93; farx(n+1)=6.067288e-001; foe(n+1)=1.175842e+000; krok(n+1)=3.032635e-003; ng(n+1)=8.025839e+002;
n=94; farx(n+1)=5.955852e-001; foe(n+1)=1.143631e+000; krok(n+1)=2.051078e-003; ng(n+1)=2.184681e+002;
n=95; farx(n+1)=5.857277e-001; foe(n+1)=1.100014e+000; krok(n+1)=4.221448e-003; ng(n+1)=4.201179e+002;
n=96; farx(n+1)=5.547197e-001; foe(n+1)=1.044562e+000; krok(n+1)=1.181377e-002; ng(n+1)=2.880660e+002;
n=97; farx(n+1)=5.446918e-001; foe(n+1)=9.768674e-001; krok(n+1)=1.522640e-002; ng(n+1)=3.933180e+002;
n=98; farx(n+1)=5.487805e-001; foe(n+1)=9.418945e-001; krok(n+1)=7.191676e-003; ng(n+1)=3.410109e+002;
n=99; farx(n+1)=5.516166e-001; foe(n+1)=9.355932e-001; krok(n+1)=3.225080e-003; ng(n+1)=2.362767e+002;
n=100; farx(n+1)=5.569658e-001; foe(n+1)=8.881826e-001; krok(n+1)=4.171530e-002; ng(n+1)=2.434237e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
