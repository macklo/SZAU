%uczenie predyktora oe
clear all;
n=0; farx(n+1)=3.981474e+003; foe(n+1)=3.945250e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.951694e+003; foe(n+1)=2.889532e+003; krok(n+1)=5.192137e-004; ng(n+1)=2.212335e+003;
n=2; farx(n+1)=1.044100e+003; foe(n+1)=7.642522e+002; krok(n+1)=3.427028e-003; ng(n+1)=1.292605e+003;
n=3; farx(n+1)=1.029442e+003; foe(n+1)=4.754651e+002; krok(n+1)=4.347668e-004; ng(n+1)=2.685698e+003;
n=4; farx(n+1)=1.148897e+003; foe(n+1)=4.488201e+002; krok(n+1)=2.957764e-003; ng(n+1)=7.796854e+002;
n=5; farx(n+1)=8.023481e+002; foe(n+1)=4.265033e+002; krok(n+1)=3.272337e-003; ng(n+1)=2.308081e+002;
n=6; farx(n+1)=6.749776e+002; foe(n+1)=3.782434e+002; krok(n+1)=5.946996e-003; ng(n+1)=5.251940e+002;
n=7; farx(n+1)=6.321875e+002; foe(n+1)=3.702148e+002; krok(n+1)=6.767654e-004; ng(n+1)=6.629783e+002;
n=8; farx(n+1)=6.173960e+002; foe(n+1)=3.677950e+002; krok(n+1)=5.849663e-003; ng(n+1)=7.055454e+002;
n=9; farx(n+1)=4.715484e+002; foe(n+1)=3.501847e+002; krok(n+1)=1.142218e-002; ng(n+1)=7.321371e+002;
n=10; farx(n+1)=4.050008e+002; foe(n+1)=3.410063e+002; krok(n+1)=1.121665e-002; ng(n+1)=8.103180e+002;
n=11; farx(n+1)=4.721471e+002; foe(n+1)=3.344048e+002; krok(n+1)=2.876670e-002; ng(n+1)=4.321707e+002;
n=12; farx(n+1)=4.606763e+002; foe(n+1)=3.319665e+002; krok(n+1)=5.536703e-002; ng(n+1)=5.181797e+002;
n=13; farx(n+1)=4.110968e+002; foe(n+1)=3.225961e+002; krok(n+1)=1.591991e-001; ng(n+1)=4.897065e+002;
n=14; farx(n+1)=2.833308e+002; foe(n+1)=2.873505e+002; krok(n+1)=1.289896e+000; ng(n+1)=2.373353e+002;
n=15; farx(n+1)=2.618500e+002; foe(n+1)=2.835578e+002; krok(n+1)=1.312690e-001; ng(n+1)=3.464540e+002;
n=16; farx(n+1)=1.171373e+002; foe(n+1)=2.490150e+002; krok(n+1)=9.859540e-001; ng(n+1)=3.796344e+002;
n=17; farx(n+1)=1.105729e+002; foe(n+1)=2.454950e+002; krok(n+1)=6.533169e-002; ng(n+1)=7.913107e+002;
n=18; farx(n+1)=1.071531e+002; foe(n+1)=2.287626e+002; krok(n+1)=2.161212e-001; ng(n+1)=8.789424e+002;
n=19; farx(n+1)=8.389713e+001; foe(n+1)=2.062704e+002; krok(n+1)=2.278082e-001; ng(n+1)=7.061661e+002;
n=20; farx(n+1)=7.443830e+001; foe(n+1)=1.969077e+002; krok(n+1)=9.048247e-002; ng(n+1)=1.018556e+003;
n=21; farx(n+1)=2.481238e+001; foe(n+1)=1.493317e+002; krok(n+1)=3.550835e-001; ng(n+1)=1.527764e+003;
n=22; farx(n+1)=2.021679e+001; foe(n+1)=1.439492e+002; krok(n+1)=1.254519e-001; ng(n+1)=1.022001e+003;
n=23; farx(n+1)=2.300979e+001; foe(n+1)=1.314823e+002; krok(n+1)=1.082372e-001; ng(n+1)=8.960143e+002;
n=24; farx(n+1)=2.437898e+001; foe(n+1)=1.186862e+002; krok(n+1)=3.556863e-001; ng(n+1)=1.822092e+003;
n=25; farx(n+1)=2.040597e+001; foe(n+1)=1.046124e+002; krok(n+1)=1.601772e+000; ng(n+1)=3.392068e+002;
%odnowa zmiennej metryki
n=26; farx(n+1)=1.951287e+001; foe(n+1)=9.992388e+001; krok(n+1)=5.921156e-006; ng(n+1)=1.901831e+003;
n=27; farx(n+1)=1.840013e+001; foe(n+1)=9.933300e+001; krok(n+1)=4.279756e-005; ng(n+1)=3.096484e+002;
n=28; farx(n+1)=1.856212e+001; foe(n+1)=9.689618e+001; krok(n+1)=7.315478e-004; ng(n+1)=1.936614e+002;
n=29; farx(n+1)=1.790114e+001; foe(n+1)=9.076361e+001; krok(n+1)=9.869792e-004; ng(n+1)=4.148689e+002;
n=30; farx(n+1)=1.877127e+001; foe(n+1)=8.903937e+001; krok(n+1)=1.336585e-003; ng(n+1)=8.281470e+002;
n=31; farx(n+1)=1.784814e+001; foe(n+1)=8.821752e+001; krok(n+1)=1.651908e-003; ng(n+1)=9.940741e+002;
n=32; farx(n+1)=1.723901e+001; foe(n+1)=8.629945e+001; krok(n+1)=1.859958e-002; ng(n+1)=1.394314e+003;
n=33; farx(n+1)=1.292649e+001; foe(n+1)=7.703413e+001; krok(n+1)=7.162643e-003; ng(n+1)=1.329420e+003;
n=34; farx(n+1)=1.222690e+001; foe(n+1)=7.635023e+001; krok(n+1)=1.526430e-002; ng(n+1)=1.837289e+003;
n=35; farx(n+1)=1.070097e+001; foe(n+1)=7.434986e+001; krok(n+1)=3.979978e-002; ng(n+1)=1.684062e+003;
n=36; farx(n+1)=8.544593e+000; foe(n+1)=7.124744e+001; krok(n+1)=1.684825e-002; ng(n+1)=1.158102e+003;
n=37; farx(n+1)=6.033507e+000; foe(n+1)=5.571119e+001; krok(n+1)=7.847380e-001; ng(n+1)=1.050068e+003;
n=38; farx(n+1)=5.653623e+000; foe(n+1)=3.627059e+001; krok(n+1)=2.758349e-001; ng(n+1)=1.165086e+003;
n=39; farx(n+1)=5.619466e+000; foe(n+1)=3.613170e+001; krok(n+1)=3.187345e-003; ng(n+1)=7.457586e+002;
n=40; farx(n+1)=5.622073e+000; foe(n+1)=3.220047e+001; krok(n+1)=5.851700e-002; ng(n+1)=7.256958e+002;
n=41; farx(n+1)=5.613699e+000; foe(n+1)=2.553399e+001; krok(n+1)=5.154350e-001; ng(n+1)=8.606117e+002;
n=42; farx(n+1)=5.649936e+000; foe(n+1)=2.476543e+001; krok(n+1)=5.734876e-002; ng(n+1)=7.240369e+002;
n=43; farx(n+1)=5.931096e+000; foe(n+1)=2.286136e+001; krok(n+1)=5.558094e-002; ng(n+1)=6.534078e+002;
n=44; farx(n+1)=5.885610e+000; foe(n+1)=2.213071e+001; krok(n+1)=2.301336e-001; ng(n+1)=4.853252e+002;
n=45; farx(n+1)=5.381663e+000; foe(n+1)=2.061951e+001; krok(n+1)=4.641165e-001; ng(n+1)=7.676831e+002;
n=46; farx(n+1)=5.246448e+000; foe(n+1)=1.950330e+001; krok(n+1)=1.210207e+000; ng(n+1)=2.840413e+002;
n=47; farx(n+1)=4.472995e+000; foe(n+1)=1.847169e+001; krok(n+1)=8.812124e-001; ng(n+1)=4.333287e+002;
n=48; farx(n+1)=4.138629e+000; foe(n+1)=1.754578e+001; krok(n+1)=1.635976e+000; ng(n+1)=3.093173e+002;
n=49; farx(n+1)=3.988698e+000; foe(n+1)=1.696457e+001; krok(n+1)=4.063898e-001; ng(n+1)=1.005872e+003;
n=50; farx(n+1)=3.853063e+000; foe(n+1)=1.654286e+001; krok(n+1)=2.196657e+000; ng(n+1)=3.482293e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=3.854548e+000; foe(n+1)=1.653285e+001; krok(n+1)=1.460365e-006; ng(n+1)=1.720886e+002;
n=52; farx(n+1)=3.839770e+000; foe(n+1)=1.648037e+001; krok(n+1)=1.069939e-005; ng(n+1)=1.524350e+002;
n=53; farx(n+1)=3.810593e+000; foe(n+1)=1.643146e+001; krok(n+1)=6.745433e-005; ng(n+1)=7.176286e+001;
n=54; farx(n+1)=3.792071e+000; foe(n+1)=1.627221e+001; krok(n+1)=1.433605e-004; ng(n+1)=6.083324e+001;
n=55; farx(n+1)=3.808118e+000; foe(n+1)=1.598361e+001; krok(n+1)=1.270603e-002; ng(n+1)=5.464604e+001;
n=56; farx(n+1)=3.873106e+000; foe(n+1)=1.591742e+001; krok(n+1)=9.751354e-004; ng(n+1)=1.433115e+002;
n=57; farx(n+1)=3.873333e+000; foe(n+1)=1.591457e+001; krok(n+1)=2.286608e-003; ng(n+1)=1.702788e+002;
n=58; farx(n+1)=3.882059e+000; foe(n+1)=1.589377e+001; krok(n+1)=1.620122e-003; ng(n+1)=1.809610e+002;
n=59; farx(n+1)=3.895010e+000; foe(n+1)=1.586167e+001; krok(n+1)=6.852302e-002; ng(n+1)=2.061816e+002;
n=60; farx(n+1)=3.953023e+000; foe(n+1)=1.573065e+001; krok(n+1)=1.248868e-001; ng(n+1)=2.156973e+002;
n=61; farx(n+1)=3.980981e+000; foe(n+1)=1.563052e+001; krok(n+1)=1.571994e-002; ng(n+1)=2.193960e+002;
n=62; farx(n+1)=3.945548e+000; foe(n+1)=1.552961e+001; krok(n+1)=7.485385e-002; ng(n+1)=1.993058e+002;
n=63; farx(n+1)=4.031227e+000; foe(n+1)=1.548764e+001; krok(n+1)=2.640812e-001; ng(n+1)=7.918593e+001;
n=64; farx(n+1)=4.088289e+000; foe(n+1)=1.545727e+001; krok(n+1)=1.072007e+000; ng(n+1)=1.068737e+002;
n=65; farx(n+1)=4.134033e+000; foe(n+1)=1.543461e+001; krok(n+1)=1.893841e+000; ng(n+1)=1.073150e+002;
n=66; farx(n+1)=4.245926e+000; foe(n+1)=1.540926e+001; krok(n+1)=3.138952e+000; ng(n+1)=1.082830e+001;
n=67; farx(n+1)=4.318005e+000; foe(n+1)=1.539468e+001; krok(n+1)=1.334889e+000; ng(n+1)=4.743873e+001;
n=68; farx(n+1)=4.376395e+000; foe(n+1)=1.538651e+001; krok(n+1)=1.127230e+000; ng(n+1)=7.434075e+001;
n=69; farx(n+1)=4.403935e+000; foe(n+1)=1.538310e+001; krok(n+1)=4.635137e-001; ng(n+1)=5.150850e+001;
n=70; farx(n+1)=4.431798e+000; foe(n+1)=1.538063e+001; krok(n+1)=4.353621e-001; ng(n+1)=7.020472e+001;
n=71; farx(n+1)=4.427556e+000; foe(n+1)=1.537936e+001; krok(n+1)=1.003200e+000; ng(n+1)=5.109998e+001;
n=72; farx(n+1)=4.438692e+000; foe(n+1)=1.537898e+001; krok(n+1)=8.067975e-001; ng(n+1)=2.687454e+000;
n=73; farx(n+1)=4.448048e+000; foe(n+1)=1.537884e+001; krok(n+1)=1.441640e+000; ng(n+1)=9.401440e+000;
n=74; farx(n+1)=4.451889e+000; foe(n+1)=1.537873e+001; krok(n+1)=2.994570e+000; ng(n+1)=3.073252e+000;
n=75; farx(n+1)=4.454138e+000; foe(n+1)=1.537848e+001; krok(n+1)=3.600437e+000; ng(n+1)=6.462505e+000;
%odnowa zmiennej metryki
n=76; farx(n+1)=4.454111e+000; foe(n+1)=1.537846e+001; krok(n+1)=2.455049e-006; ng(n+1)=7.327617e+000;
n=77; farx(n+1)=4.454020e+000; foe(n+1)=1.537840e+001; krok(n+1)=9.165986e-006; ng(n+1)=5.494103e+000;
n=78; farx(n+1)=4.453912e+000; foe(n+1)=1.537840e+001; krok(n+1)=5.210339e-005; ng(n+1)=7.817432e-001;
n=79; farx(n+1)=4.454449e+000; foe(n+1)=1.537839e+001; krok(n+1)=1.769627e-004; ng(n+1)=5.423686e-001;
n=80; farx(n+1)=4.455572e+000; foe(n+1)=1.537837e+001; krok(n+1)=2.749912e-004; ng(n+1)=7.899609e-001;
n=81; farx(n+1)=4.455729e+000; foe(n+1)=1.537834e+001; krok(n+1)=3.176507e-003; ng(n+1)=1.987814e-001;
n=82; farx(n+1)=4.453552e+000; foe(n+1)=1.537828e+001; krok(n+1)=2.621411e-003; ng(n+1)=3.707693e-001;
n=83; farx(n+1)=4.452967e+000; foe(n+1)=1.537826e+001; krok(n+1)=1.154146e-002; ng(n+1)=1.930262e-001;
n=84; farx(n+1)=4.452859e+000; foe(n+1)=1.537826e+001; krok(n+1)=1.750226e-002; ng(n+1)=1.864785e-001;
n=85; farx(n+1)=4.452173e+000; foe(n+1)=1.537826e+001; krok(n+1)=2.835564e-002; ng(n+1)=1.780329e-001;
n=86; farx(n+1)=4.453005e+000; foe(n+1)=1.537818e+001; krok(n+1)=3.071879e+000; ng(n+1)=1.689367e-001;
n=87; farx(n+1)=4.468571e+000; foe(n+1)=1.537729e+001; krok(n+1)=3.691454e-001; ng(n+1)=6.513825e-001;
n=88; farx(n+1)=4.473940e+000; foe(n+1)=1.537641e+001; krok(n+1)=2.421479e+000; ng(n+1)=1.978342e+001;
n=89; farx(n+1)=4.467553e+000; foe(n+1)=1.537608e+001; krok(n+1)=7.101669e-001; ng(n+1)=1.135884e+001;
n=90; farx(n+1)=4.460184e+000; foe(n+1)=1.537556e+001; krok(n+1)=1.157348e+000; ng(n+1)=9.819422e+000;
n=91; farx(n+1)=4.448781e+000; foe(n+1)=1.537531e+001; krok(n+1)=1.728528e+000; ng(n+1)=2.877379e+000;
n=92; farx(n+1)=4.445425e+000; foe(n+1)=1.537525e+001; krok(n+1)=1.001473e+000; ng(n+1)=7.725651e+000;
n=93; farx(n+1)=4.443529e+000; foe(n+1)=1.537524e+001; krok(n+1)=8.922169e-001; ng(n+1)=2.540272e+000;
n=94; farx(n+1)=4.442099e+000; foe(n+1)=1.537523e+001; krok(n+1)=1.119298e+000; ng(n+1)=7.513011e-001;
n=95; farx(n+1)=4.441787e+000; foe(n+1)=1.537523e+001; krok(n+1)=2.196657e+000; ng(n+1)=1.182504e+000;
n=96; farx(n+1)=4.441705e+000; foe(n+1)=1.537523e+001; krok(n+1)=1.043878e+000; ng(n+1)=9.396045e-001;
n=97; farx(n+1)=4.441757e+000; foe(n+1)=1.537523e+001; krok(n+1)=9.362720e-001; ng(n+1)=5.483563e-002;
n=98; farx(n+1)=4.441757e+000; foe(n+1)=1.537523e+001; krok(n+1)=3.287163e-005; ng(n+1)=2.129107e-002;
n=99; farx(n+1)=4.441757e+000; foe(n+1)=1.537523e+001; krok(n+1)=1.180110e-008; ng(n+1)=2.129039e-002;
n=100; farx(n+1)=4.441757e+000; foe(n+1)=1.537523e+001; krok(n+1)=4.300043e-006; ng(n+1)=2.129039e-002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
