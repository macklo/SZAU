%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.596397e+003; foe(n+1)=4.726998e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.671236e+003; foe(n+1)=3.696688e+003; krok(n+1)=3.467612e-004; ng(n+1)=5.012421e+003;
n=2; farx(n+1)=1.050486e+003; foe(n+1)=1.124338e+003; krok(n+1)=1.441974e-003; ng(n+1)=3.787900e+003;
n=3; farx(n+1)=1.209255e+003; foe(n+1)=8.563897e+002; krok(n+1)=1.711903e-004; ng(n+1)=7.445451e+003;
n=4; farx(n+1)=1.482411e+003; foe(n+1)=7.819908e+002; krok(n+1)=7.433745e-004; ng(n+1)=4.649800e+003;
n=5; farx(n+1)=7.881424e+002; foe(n+1)=6.109956e+002; krok(n+1)=3.401267e-003; ng(n+1)=8.363506e+002;
n=6; farx(n+1)=3.335024e+002; foe(n+1)=4.356396e+002; krok(n+1)=8.500296e-004; ng(n+1)=3.465518e+003;
n=7; farx(n+1)=2.718046e+002; foe(n+1)=4.011266e+002; krok(n+1)=1.673718e-004; ng(n+1)=3.985083e+003;
n=8; farx(n+1)=2.716694e+002; foe(n+1)=3.860809e+002; krok(n+1)=3.721397e-003; ng(n+1)=3.947250e+003;
n=9; farx(n+1)=2.778739e+002; foe(n+1)=3.763487e+002; krok(n+1)=8.456028e-004; ng(n+1)=1.350278e+003;
n=10; farx(n+1)=2.760800e+002; foe(n+1)=3.736097e+002; krok(n+1)=9.249103e-004; ng(n+1)=1.240595e+003;
n=11; farx(n+1)=2.824322e+002; foe(n+1)=3.616831e+002; krok(n+1)=2.032216e-003; ng(n+1)=1.004150e+003;
n=12; farx(n+1)=2.256359e+002; foe(n+1)=3.405172e+002; krok(n+1)=4.777659e-003; ng(n+1)=8.762112e+002;
n=13; farx(n+1)=1.990873e+002; foe(n+1)=3.307915e+002; krok(n+1)=7.209870e-004; ng(n+1)=2.207505e+003;
n=14; farx(n+1)=1.534345e+002; foe(n+1)=3.065434e+002; krok(n+1)=9.854180e-004; ng(n+1)=2.865903e+003;
n=15; farx(n+1)=8.791497e+001; foe(n+1)=2.671293e+002; krok(n+1)=8.594024e-004; ng(n+1)=2.573575e+003;
n=16; farx(n+1)=6.849267e+001; foe(n+1)=2.559991e+002; krok(n+1)=3.534374e-004; ng(n+1)=2.769978e+003;
n=17; farx(n+1)=5.676400e+001; foe(n+1)=2.452591e+002; krok(n+1)=4.988661e-004; ng(n+1)=3.038041e+003;
n=18; farx(n+1)=5.162223e+001; foe(n+1)=2.264148e+002; krok(n+1)=6.946998e-004; ng(n+1)=2.420703e+003;
n=19; farx(n+1)=4.353404e+001; foe(n+1)=2.046579e+002; krok(n+1)=1.303632e-003; ng(n+1)=2.969295e+003;
n=20; farx(n+1)=4.507586e+001; foe(n+1)=1.841400e+002; krok(n+1)=1.219912e-003; ng(n+1)=4.273642e+003;
n=21; farx(n+1)=4.095955e+001; foe(n+1)=1.627387e+002; krok(n+1)=4.740993e-004; ng(n+1)=5.903441e+003;
n=22; farx(n+1)=4.342787e+001; foe(n+1)=1.454830e+002; krok(n+1)=4.465626e-004; ng(n+1)=5.676060e+003;
n=23; farx(n+1)=4.118159e+001; foe(n+1)=1.378811e+002; krok(n+1)=8.605590e-004; ng(n+1)=6.088477e+003;
n=24; farx(n+1)=4.149441e+001; foe(n+1)=1.328808e+002; krok(n+1)=9.232346e-004; ng(n+1)=7.036833e+003;
n=25; farx(n+1)=4.450977e+001; foe(n+1)=1.314126e+002; krok(n+1)=1.181843e-003; ng(n+1)=6.117508e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=3.859132e+001; foe(n+1)=1.207219e+002; krok(n+1)=7.048309e-006; ng(n+1)=5.029588e+003;
n=27; farx(n+1)=3.481504e+001; foe(n+1)=1.172422e+002; krok(n+1)=2.444110e-005; ng(n+1)=1.392982e+003;
n=28; farx(n+1)=2.769185e+001; foe(n+1)=1.087667e+002; krok(n+1)=1.592886e-004; ng(n+1)=1.096506e+003;
n=29; farx(n+1)=1.830242e+001; foe(n+1)=9.209949e+001; krok(n+1)=3.403614e-004; ng(n+1)=1.181002e+003;
n=30; farx(n+1)=1.617722e+001; foe(n+1)=8.177346e+001; krok(n+1)=7.138862e-004; ng(n+1)=1.792186e+003;
n=31; farx(n+1)=1.600111e+001; foe(n+1)=6.286969e+001; krok(n+1)=1.244450e-003; ng(n+1)=7.197155e+003;
n=32; farx(n+1)=1.630091e+001; foe(n+1)=5.257135e+001; krok(n+1)=4.216476e-004; ng(n+1)=3.512090e+003;
n=33; farx(n+1)=1.484518e+001; foe(n+1)=4.673126e+001; krok(n+1)=4.630782e-004; ng(n+1)=3.550540e+003;
n=34; farx(n+1)=1.422050e+001; foe(n+1)=4.055742e+001; krok(n+1)=7.394073e-004; ng(n+1)=1.442889e+003;
n=35; farx(n+1)=1.439680e+001; foe(n+1)=3.765820e+001; krok(n+1)=1.291963e-003; ng(n+1)=2.155536e+003;
n=36; farx(n+1)=1.468842e+001; foe(n+1)=3.193761e+001; krok(n+1)=2.027075e-003; ng(n+1)=1.197579e+003;
n=37; farx(n+1)=1.382667e+001; foe(n+1)=2.978604e+001; krok(n+1)=2.883948e-003; ng(n+1)=1.322165e+003;
n=38; farx(n+1)=1.250194e+001; foe(n+1)=2.647865e+001; krok(n+1)=1.882699e-003; ng(n+1)=1.691842e+003;
n=39; farx(n+1)=1.147804e+001; foe(n+1)=2.388043e+001; krok(n+1)=4.375565e-003; ng(n+1)=2.773189e+003;
n=40; farx(n+1)=1.095120e+001; foe(n+1)=2.279309e+001; krok(n+1)=3.868579e-003; ng(n+1)=1.007043e+003;
n=41; farx(n+1)=9.187971e+000; foe(n+1)=2.097401e+001; krok(n+1)=1.333847e-002; ng(n+1)=5.706830e+002;
n=42; farx(n+1)=6.996585e+000; foe(n+1)=1.776857e+001; krok(n+1)=2.165649e-002; ng(n+1)=7.654790e+002;
n=43; farx(n+1)=6.279249e+000; foe(n+1)=1.651554e+001; krok(n+1)=8.948804e-003; ng(n+1)=3.414463e+002;
n=44; farx(n+1)=4.918690e+000; foe(n+1)=1.377759e+001; krok(n+1)=1.140023e-002; ng(n+1)=1.281394e+003;
n=45; farx(n+1)=4.554698e+000; foe(n+1)=1.295310e+001; krok(n+1)=3.572423e-003; ng(n+1)=2.219738e+003;
n=46; farx(n+1)=4.199662e+000; foe(n+1)=1.216403e+001; krok(n+1)=1.095618e-002; ng(n+1)=1.187338e+003;
n=47; farx(n+1)=3.851295e+000; foe(n+1)=1.130512e+001; krok(n+1)=1.432529e-002; ng(n+1)=1.079920e+003;
n=48; farx(n+1)=3.695445e+000; foe(n+1)=1.091840e+001; krok(n+1)=3.544455e-003; ng(n+1)=7.036847e+002;
n=49; farx(n+1)=3.108716e+000; foe(n+1)=8.477767e+000; krok(n+1)=4.305393e-002; ng(n+1)=6.149426e+002;
n=50; farx(n+1)=3.049669e+000; foe(n+1)=8.190788e+000; krok(n+1)=3.881663e-003; ng(n+1)=4.902981e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=3.025255e+000; foe(n+1)=8.108002e+000; krok(n+1)=4.022035e-006; ng(n+1)=5.857776e+002;
n=52; farx(n+1)=3.020086e+000; foe(n+1)=8.084407e+000; krok(n+1)=4.264017e-005; ng(n+1)=1.125660e+002;
n=53; farx(n+1)=2.969255e+000; foe(n+1)=7.991784e+000; krok(n+1)=7.454994e-005; ng(n+1)=1.622779e+002;
n=54; farx(n+1)=2.878209e+000; foe(n+1)=7.742155e+000; krok(n+1)=1.413318e-004; ng(n+1)=2.196506e+002;
n=55; farx(n+1)=2.818017e+000; foe(n+1)=7.550996e+000; krok(n+1)=9.988735e-004; ng(n+1)=8.205341e+001;
n=56; farx(n+1)=2.672097e+000; foe(n+1)=7.034743e+000; krok(n+1)=3.113587e-003; ng(n+1)=1.032257e+002;
n=57; farx(n+1)=2.665428e+000; foe(n+1)=6.888763e+000; krok(n+1)=1.338974e-003; ng(n+1)=5.726434e+002;
n=58; farx(n+1)=2.816159e+000; foe(n+1)=6.400317e+000; krok(n+1)=1.095618e-002; ng(n+1)=7.239318e+002;
n=59; farx(n+1)=2.803264e+000; foe(n+1)=6.117706e+000; krok(n+1)=6.956269e-003; ng(n+1)=1.328233e+002;
n=60; farx(n+1)=2.742111e+000; foe(n+1)=6.019270e+000; krok(n+1)=2.284436e-002; ng(n+1)=1.317402e+002;
n=61; farx(n+1)=2.682006e+000; foe(n+1)=5.962916e+000; krok(n+1)=1.270603e-002; ng(n+1)=2.422334e+002;
n=62; farx(n+1)=2.592839e+000; foe(n+1)=5.846439e+000; krok(n+1)=1.736512e-002; ng(n+1)=4.838314e+002;
n=63; farx(n+1)=2.383822e+000; foe(n+1)=5.677548e+000; krok(n+1)=1.269968e-002; ng(n+1)=4.576775e+002;
n=64; farx(n+1)=2.258400e+000; foe(n+1)=5.506496e+000; krok(n+1)=5.139053e-003; ng(n+1)=5.723387e+002;
n=65; farx(n+1)=2.185277e+000; foe(n+1)=5.431359e+000; krok(n+1)=5.945944e-003; ng(n+1)=4.871765e+002;
n=66; farx(n+1)=2.029697e+000; foe(n+1)=5.166494e+000; krok(n+1)=2.107682e-002; ng(n+1)=4.594366e+002;
n=67; farx(n+1)=2.009382e+000; foe(n+1)=5.105559e+000; krok(n+1)=1.219985e-002; ng(n+1)=3.032112e+002;
n=68; farx(n+1)=1.877765e+000; foe(n+1)=4.904177e+000; krok(n+1)=2.173419e-002; ng(n+1)=4.957066e+002;
n=69; farx(n+1)=1.757415e+000; foe(n+1)=4.663438e+000; krok(n+1)=1.417782e-002; ng(n+1)=7.162600e+002;
n=70; farx(n+1)=1.671948e+000; foe(n+1)=4.525708e+000; krok(n+1)=1.417782e-002; ng(n+1)=5.403506e+002;
n=71; farx(n+1)=1.597187e+000; foe(n+1)=4.345191e+000; krok(n+1)=4.681906e-002; ng(n+1)=7.033949e+002;
n=72; farx(n+1)=1.499440e+000; foe(n+1)=4.136166e+000; krok(n+1)=3.864513e-002; ng(n+1)=4.990418e+002;
n=73; farx(n+1)=1.484638e+000; foe(n+1)=4.074143e+000; krok(n+1)=3.500452e-002; ng(n+1)=3.099058e+002;
n=74; farx(n+1)=1.427530e+000; foe(n+1)=3.900717e+000; krok(n+1)=5.715876e-002; ng(n+1)=1.630684e+002;
n=75; farx(n+1)=1.409144e+000; foe(n+1)=3.640270e+000; krok(n+1)=5.424665e-002; ng(n+1)=2.515561e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.401323e+000; foe(n+1)=3.579904e+000; krok(n+1)=5.674684e-006; ng(n+1)=4.886257e+002;
n=77; farx(n+1)=1.397768e+000; foe(n+1)=3.559233e+000; krok(n+1)=5.288329e-006; ng(n+1)=3.074611e+002;
n=78; farx(n+1)=1.371710e+000; foe(n+1)=3.456769e+000; krok(n+1)=3.308775e-004; ng(n+1)=8.023832e+001;
n=79; farx(n+1)=1.372516e+000; foe(n+1)=3.433192e+000; krok(n+1)=5.234435e-005; ng(n+1)=1.130032e+002;
n=80; farx(n+1)=1.351227e+000; foe(n+1)=3.368017e+000; krok(n+1)=1.463096e-003; ng(n+1)=3.766557e+001;
n=81; farx(n+1)=1.297022e+000; foe(n+1)=3.256993e+000; krok(n+1)=2.007255e-003; ng(n+1)=4.448526e+001;
n=82; farx(n+1)=1.294306e+000; foe(n+1)=3.221925e+000; krok(n+1)=3.740913e-003; ng(n+1)=3.402811e+002;
n=83; farx(n+1)=1.286012e+000; foe(n+1)=3.168370e+000; krok(n+1)=3.096185e-003; ng(n+1)=4.669971e+002;
n=84; farx(n+1)=1.279104e+000; foe(n+1)=3.098290e+000; krok(n+1)=6.283621e-003; ng(n+1)=6.377817e+002;
n=85; farx(n+1)=1.285414e+000; foe(n+1)=2.939082e+000; krok(n+1)=1.251384e-002; ng(n+1)=7.462854e+002;
n=86; farx(n+1)=1.307499e+000; foe(n+1)=2.810571e+000; krok(n+1)=7.321912e-003; ng(n+1)=1.670565e+002;
n=87; farx(n+1)=1.279572e+000; foe(n+1)=2.733073e+000; krok(n+1)=3.280985e-002; ng(n+1)=2.173190e+002;
n=88; farx(n+1)=1.217943e+000; foe(n+1)=2.669432e+000; krok(n+1)=1.688579e-002; ng(n+1)=3.053726e+002;
n=89; farx(n+1)=1.191208e+000; foe(n+1)=2.592090e+000; krok(n+1)=9.939766e-003; ng(n+1)=2.213970e+002;
n=90; farx(n+1)=1.108619e+000; foe(n+1)=2.400983e+000; krok(n+1)=1.086709e-002; ng(n+1)=3.567018e+002;
n=91; farx(n+1)=1.075741e+000; foe(n+1)=2.341403e+000; krok(n+1)=1.086709e-002; ng(n+1)=4.486174e+002;
n=92; farx(n+1)=1.046488e+000; foe(n+1)=2.255825e+000; krok(n+1)=1.859267e-002; ng(n+1)=4.153927e+002;
n=93; farx(n+1)=1.033885e+000; foe(n+1)=2.199851e+000; krok(n+1)=1.464382e-002; ng(n+1)=4.088002e+002;
n=94; farx(n+1)=1.018601e+000; foe(n+1)=2.095869e+000; krok(n+1)=2.414071e-002; ng(n+1)=3.178323e+002;
n=95; farx(n+1)=9.612592e-001; foe(n+1)=2.022314e+000; krok(n+1)=4.000365e-002; ng(n+1)=9.268125e+001;
n=96; farx(n+1)=9.356815e-001; foe(n+1)=1.971272e+000; krok(n+1)=1.214360e-002; ng(n+1)=5.224748e+002;
n=97; farx(n+1)=9.276783e-001; foe(n+1)=1.942107e+000; krok(n+1)=1.021537e-002; ng(n+1)=3.313288e+002;
n=98; farx(n+1)=8.899401e-001; foe(n+1)=1.828526e+000; krok(n+1)=6.402713e-002; ng(n+1)=3.163337e+002;
n=99; farx(n+1)=8.316473e-001; foe(n+1)=1.740862e+000; krok(n+1)=4.356626e-002; ng(n+1)=4.294100e+002;
n=100; farx(n+1)=7.668242e-001; foe(n+1)=1.685934e+000; krok(n+1)=3.039298e-002; ng(n+1)=2.551305e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
