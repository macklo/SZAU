%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.771534e+003; foe(n+1)=4.859559e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.708108e+003; foe(n+1)=3.803684e+003; krok(n+1)=3.383827e-004; ng(n+1)=2.677922e+003;
n=2; farx(n+1)=8.105071e+002; foe(n+1)=9.531440e+002; krok(n+1)=2.531421e-003; ng(n+1)=1.308654e+003;
n=3; farx(n+1)=8.326089e+002; foe(n+1)=7.559624e+002; krok(n+1)=1.046887e-004; ng(n+1)=3.920855e+003;
n=4; farx(n+1)=9.647280e+002; foe(n+1)=6.920280e+002; krok(n+1)=1.889294e-003; ng(n+1)=3.137455e+003;
n=5; farx(n+1)=6.751012e+002; foe(n+1)=5.472716e+002; krok(n+1)=4.161987e-003; ng(n+1)=1.447412e+003;
n=6; farx(n+1)=5.678909e+002; foe(n+1)=5.282837e+002; krok(n+1)=5.631478e-004; ng(n+1)=7.159706e+002;
n=7; farx(n+1)=5.186645e+002; foe(n+1)=5.123018e+002; krok(n+1)=1.308616e-002; ng(n+1)=9.594584e+002;
n=8; farx(n+1)=4.840851e+002; foe(n+1)=5.015669e+002; krok(n+1)=5.704577e-004; ng(n+1)=8.908938e+002;
n=9; farx(n+1)=4.311621e+002; foe(n+1)=4.844513e+002; krok(n+1)=1.033570e-002; ng(n+1)=7.403156e+002;
n=10; farx(n+1)=3.762038e+002; foe(n+1)=4.623301e+002; krok(n+1)=1.256724e-002; ng(n+1)=6.316617e+002;
n=11; farx(n+1)=3.048248e+002; foe(n+1)=4.163543e+002; krok(n+1)=8.776600e-002; ng(n+1)=8.105793e+002;
n=12; farx(n+1)=2.738254e+002; foe(n+1)=3.958324e+002; krok(n+1)=1.046893e-001; ng(n+1)=1.040392e+003;
n=13; farx(n+1)=2.007793e+002; foe(n+1)=3.328152e+002; krok(n+1)=2.595550e-001; ng(n+1)=1.461696e+003;
n=14; farx(n+1)=2.159334e+002; foe(n+1)=3.250947e+002; krok(n+1)=7.596282e-002; ng(n+1)=1.322929e+003;
n=15; farx(n+1)=1.260421e+002; foe(n+1)=2.969484e+002; krok(n+1)=9.270274e-001; ng(n+1)=5.680179e+002;
n=16; farx(n+1)=8.091907e+001; foe(n+1)=2.765080e+002; krok(n+1)=1.448235e-001; ng(n+1)=2.320871e+003;
n=17; farx(n+1)=5.078481e+001; foe(n+1)=2.545660e+002; krok(n+1)=1.455442e-001; ng(n+1)=3.972263e+003;
n=18; farx(n+1)=4.598695e+001; foe(n+1)=2.500014e+002; krok(n+1)=8.685192e-003; ng(n+1)=5.982507e+003;
n=19; farx(n+1)=6.086995e+001; foe(n+1)=1.870500e+002; krok(n+1)=1.572879e-001; ng(n+1)=6.278815e+003;
n=20; farx(n+1)=5.933344e+001; foe(n+1)=1.722592e+002; krok(n+1)=7.113068e-002; ng(n+1)=2.168776e+003;
n=21; farx(n+1)=5.337896e+001; foe(n+1)=1.667916e+002; krok(n+1)=1.601162e-001; ng(n+1)=6.292179e+002;
n=22; farx(n+1)=3.197861e+001; foe(n+1)=1.277725e+002; krok(n+1)=2.382113e-001; ng(n+1)=5.730414e+003;
n=23; farx(n+1)=3.426867e+001; foe(n+1)=1.179649e+002; krok(n+1)=1.308949e-001; ng(n+1)=1.413678e+003;
n=24; farx(n+1)=4.104968e+001; foe(n+1)=1.042477e+002; krok(n+1)=1.981558e-001; ng(n+1)=2.653873e+003;
n=25; farx(n+1)=3.995945e+001; foe(n+1)=9.139751e+001; krok(n+1)=6.034212e-001; ng(n+1)=2.293963e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=3.987871e+001; foe(n+1)=9.069464e+001; krok(n+1)=1.028644e-006; ng(n+1)=2.125610e+003;
n=27; farx(n+1)=3.978042e+001; foe(n+1)=8.979064e+001; krok(n+1)=7.063894e-005; ng(n+1)=3.105077e+002;
n=28; farx(n+1)=3.779378e+001; foe(n+1)=8.887331e+001; krok(n+1)=1.625548e-004; ng(n+1)=2.466936e+002;
n=29; farx(n+1)=3.158129e+001; foe(n+1)=8.609267e+001; krok(n+1)=2.469223e-003; ng(n+1)=1.533305e+002;
n=30; farx(n+1)=2.618200e+001; foe(n+1)=8.088855e+001; krok(n+1)=1.718805e-003; ng(n+1)=6.255259e+002;
n=31; farx(n+1)=2.378420e+001; foe(n+1)=7.991755e+001; krok(n+1)=1.300439e-003; ng(n+1)=4.512176e+002;
n=32; farx(n+1)=2.237016e+001; foe(n+1)=7.594093e+001; krok(n+1)=9.824309e-003; ng(n+1)=5.424590e+002;
n=33; farx(n+1)=2.039820e+001; foe(n+1)=7.427337e+001; krok(n+1)=1.524981e-003; ng(n+1)=1.300679e+003;
n=34; farx(n+1)=2.009312e+001; foe(n+1)=7.105309e+001; krok(n+1)=2.629929e-002; ng(n+1)=8.295339e+002;
n=35; farx(n+1)=2.057330e+001; foe(n+1)=6.924078e+001; krok(n+1)=7.833534e-003; ng(n+1)=8.650927e+002;
n=36; farx(n+1)=2.150010e+001; foe(n+1)=6.765149e+001; krok(n+1)=1.803262e-002; ng(n+1)=1.032344e+003;
n=37; farx(n+1)=1.926440e+001; foe(n+1)=5.947582e+001; krok(n+1)=3.404866e-001; ng(n+1)=8.689920e+002;
n=38; farx(n+1)=2.357047e+001; foe(n+1)=5.438018e+001; krok(n+1)=1.600146e-001; ng(n+1)=8.968804e+002;
n=39; farx(n+1)=1.767139e+001; foe(n+1)=4.825773e+001; krok(n+1)=9.792329e-001; ng(n+1)=8.231238e+002;
n=40; farx(n+1)=1.517877e+001; foe(n+1)=4.381617e+001; krok(n+1)=1.715220e+000; ng(n+1)=7.116052e+002;
n=41; farx(n+1)=1.530565e+001; foe(n+1)=3.931639e+001; krok(n+1)=1.355274e+000; ng(n+1)=1.521265e+003;
n=42; farx(n+1)=1.260738e+001; foe(n+1)=3.509979e+001; krok(n+1)=1.240335e+000; ng(n+1)=4.841432e+002;
n=43; farx(n+1)=1.044438e+001; foe(n+1)=3.241689e+001; krok(n+1)=9.073805e-001; ng(n+1)=2.271181e+002;
n=44; farx(n+1)=9.123807e+000; foe(n+1)=3.049136e+001; krok(n+1)=8.295027e-001; ng(n+1)=2.422552e+002;
n=45; farx(n+1)=8.217491e+000; foe(n+1)=3.009840e+001; krok(n+1)=3.730916e-001; ng(n+1)=1.266814e+002;
n=46; farx(n+1)=7.844434e+000; foe(n+1)=2.994654e+001; krok(n+1)=4.714530e-001; ng(n+1)=2.724869e+002;
n=47; farx(n+1)=7.213868e+000; foe(n+1)=2.952391e+001; krok(n+1)=2.984343e+000; ng(n+1)=1.672572e+002;
n=48; farx(n+1)=6.793081e+000; foe(n+1)=2.918539e+001; krok(n+1)=6.871474e-001; ng(n+1)=8.263207e+001;
n=49; farx(n+1)=6.632399e+000; foe(n+1)=2.891147e+001; krok(n+1)=1.098329e+000; ng(n+1)=3.912311e+002;
n=50; farx(n+1)=6.065600e+000; foe(n+1)=2.875287e+001; krok(n+1)=1.155850e+000; ng(n+1)=1.372856e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=6.055932e+000; foe(n+1)=2.874476e+001; krok(n+1)=2.288495e-005; ng(n+1)=4.657609e+001;
n=52; farx(n+1)=6.050414e+000; foe(n+1)=2.873935e+001; krok(n+1)=4.654987e-006; ng(n+1)=9.451014e+001;
n=53; farx(n+1)=6.052250e+000; foe(n+1)=2.873351e+001; krok(n+1)=1.223990e-004; ng(n+1)=2.037446e+001;
n=54; farx(n+1)=5.955595e+000; foe(n+1)=2.870301e+001; krok(n+1)=3.018569e-004; ng(n+1)=2.872364e+001;
n=55; farx(n+1)=6.061953e+000; foe(n+1)=2.853624e+001; krok(n+1)=2.144917e-003; ng(n+1)=2.378381e+001;
n=56; farx(n+1)=5.956660e+000; foe(n+1)=2.831293e+001; krok(n+1)=2.564143e-003; ng(n+1)=1.388009e+002;
n=57; farx(n+1)=5.817218e+000; foe(n+1)=2.826537e+001; krok(n+1)=9.600724e-003; ng(n+1)=3.107410e+002;
n=58; farx(n+1)=5.847057e+000; foe(n+1)=2.810719e+001; krok(n+1)=1.047148e-001; ng(n+1)=3.997992e+002;
n=59; farx(n+1)=6.085390e+000; foe(n+1)=2.802756e+001; krok(n+1)=4.756755e-002; ng(n+1)=6.970007e+002;
n=60; farx(n+1)=5.496260e+000; foe(n+1)=2.772905e+001; krok(n+1)=9.205345e-001; ng(n+1)=6.322737e+002;
n=61; farx(n+1)=5.606515e+000; foe(n+1)=2.769654e+001; krok(n+1)=8.276012e-002; ng(n+1)=6.600232e+001;
n=62; farx(n+1)=5.369711e+000; foe(n+1)=2.746303e+001; krok(n+1)=9.634903e-001; ng(n+1)=6.517039e+001;
n=63; farx(n+1)=5.326849e+000; foe(n+1)=2.725678e+001; krok(n+1)=9.744898e-001; ng(n+1)=1.533950e+002;
n=64; farx(n+1)=4.910567e+000; foe(n+1)=2.701918e+001; krok(n+1)=1.006124e+000; ng(n+1)=3.683951e+002;
n=65; farx(n+1)=4.776716e+000; foe(n+1)=2.690830e+001; krok(n+1)=4.641165e-001; ng(n+1)=3.523373e+002;
n=66; farx(n+1)=4.762889e+000; foe(n+1)=2.673968e+001; krok(n+1)=3.924010e-001; ng(n+1)=1.867883e+002;
n=67; farx(n+1)=4.471033e+000; foe(n+1)=2.657944e+001; krok(n+1)=1.492172e+000; ng(n+1)=1.204846e+002;
n=68; farx(n+1)=4.195415e+000; foe(n+1)=2.646700e+001; krok(n+1)=5.692710e-001; ng(n+1)=5.942913e+002;
n=69; farx(n+1)=4.040729e+000; foe(n+1)=2.638039e+001; krok(n+1)=6.002260e-001; ng(n+1)=4.730497e+002;
n=70; farx(n+1)=4.037294e+000; foe(n+1)=2.621159e+001; krok(n+1)=7.807902e-001; ng(n+1)=3.597081e+002;
n=71; farx(n+1)=3.634951e+000; foe(n+1)=2.584131e+001; krok(n+1)=9.270274e-001; ng(n+1)=1.499567e+002;
n=72; farx(n+1)=3.123567e+000; foe(n+1)=2.548292e+001; krok(n+1)=7.113726e-001; ng(n+1)=8.251410e+002;
n=73; farx(n+1)=3.205224e+000; foe(n+1)=2.522226e+001; krok(n+1)=9.606268e-001; ng(n+1)=6.290587e+002;
n=74; farx(n+1)=3.113201e+000; foe(n+1)=2.507151e+001; krok(n+1)=1.268929e+000; ng(n+1)=9.436015e+001;
n=75; farx(n+1)=3.017952e+000; foe(n+1)=2.501893e+001; krok(n+1)=1.961845e-001; ng(n+1)=1.127020e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=3.019397e+000; foe(n+1)=2.500697e+001; krok(n+1)=9.682870e-007; ng(n+1)=2.494925e+002;
n=77; farx(n+1)=3.031654e+000; foe(n+1)=2.497941e+001; krok(n+1)=1.178746e-005; ng(n+1)=1.303132e+002;
n=78; farx(n+1)=3.033798e+000; foe(n+1)=2.497847e+001; krok(n+1)=7.961977e-005; ng(n+1)=7.316890e+000;
n=79; farx(n+1)=3.035754e+000; foe(n+1)=2.495967e+001; krok(n+1)=2.238326e-004; ng(n+1)=2.389293e+001;
n=80; farx(n+1)=3.046050e+000; foe(n+1)=2.495788e+001; krok(n+1)=3.229906e-004; ng(n+1)=6.866720e+000;
n=81; farx(n+1)=3.008447e+000; foe(n+1)=2.492724e+001; krok(n+1)=4.953534e-003; ng(n+1)=6.381294e+000;
n=82; farx(n+1)=3.052261e+000; foe(n+1)=2.491496e+001; krok(n+1)=1.086709e-002; ng(n+1)=1.245637e+001;
n=83; farx(n+1)=2.996619e+000; foe(n+1)=2.490065e+001; krok(n+1)=5.983917e-002; ng(n+1)=4.340512e+001;
n=84; farx(n+1)=2.983248e+000; foe(n+1)=2.489681e+001; krok(n+1)=5.018717e-002; ng(n+1)=2.203716e+002;
n=85; farx(n+1)=2.998732e+000; foe(n+1)=2.489545e+001; krok(n+1)=2.051314e-002; ng(n+1)=1.993435e+002;
n=86; farx(n+1)=2.890021e+000; foe(n+1)=2.486392e+001; krok(n+1)=1.991105e+000; ng(n+1)=1.778597e+002;
n=87; farx(n+1)=2.742561e+000; foe(n+1)=2.480548e+001; krok(n+1)=3.629522e+000; ng(n+1)=1.597181e+002;
n=88; farx(n+1)=2.713053e+000; foe(n+1)=2.478125e+001; krok(n+1)=8.006816e-001; ng(n+1)=7.553729e+001;
n=89; farx(n+1)=2.654955e+000; foe(n+1)=2.476062e+001; krok(n+1)=2.621075e-001; ng(n+1)=2.572296e+002;
n=90; farx(n+1)=2.646484e+000; foe(n+1)=2.473266e+001; krok(n+1)=1.310660e+000; ng(n+1)=3.507268e+002;
n=91; farx(n+1)=2.643021e+000; foe(n+1)=2.469645e+001; krok(n+1)=1.841069e+000; ng(n+1)=4.575615e+001;
n=92; farx(n+1)=2.505521e+000; foe(n+1)=2.465285e+001; krok(n+1)=1.009068e+000; ng(n+1)=9.544262e+001;
n=93; farx(n+1)=2.442352e+000; foe(n+1)=2.461525e+001; krok(n+1)=1.961613e+000; ng(n+1)=1.093235e+002;
n=94; farx(n+1)=2.413767e+000; foe(n+1)=2.457332e+001; krok(n+1)=1.585247e+000; ng(n+1)=4.109265e+002;
n=95; farx(n+1)=2.474300e+000; foe(n+1)=2.449125e+001; krok(n+1)=1.378694e+000; ng(n+1)=2.168054e+002;
n=96; farx(n+1)=2.354642e+000; foe(n+1)=2.440690e+001; krok(n+1)=1.304490e+000; ng(n+1)=2.039927e+002;
n=97; farx(n+1)=2.327456e+000; foe(n+1)=2.437843e+001; krok(n+1)=1.498210e+000; ng(n+1)=1.329657e+002;
n=98; farx(n+1)=2.355554e+000; foe(n+1)=2.436183e+001; krok(n+1)=1.420334e+000; ng(n+1)=4.530137e+001;
n=99; farx(n+1)=2.377742e+000; foe(n+1)=2.434958e+001; krok(n+1)=1.192577e+000; ng(n+1)=1.530096e+002;
n=100; farx(n+1)=2.352614e+000; foe(n+1)=2.434697e+001; krok(n+1)=1.134786e+000; ng(n+1)=2.486212e+001;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
