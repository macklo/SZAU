%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.363777e+003; foe(n+1)=4.384953e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
n=1; farx(n+1)=2.917408e+003; foe(n+1)=2.940115e+003; krok(n+1)=5.081655e-004; ng(n+1)=2.153559e+003;
n=2; farx(n+1)=4.794214e+002; foe(n+1)=4.769271e+002; krok(n+1)=4.210497e-003; ng(n+1)=7.728291e+002;
n=3; farx(n+1)=4.894001e+002; foe(n+1)=4.572310e+002; krok(n+1)=3.729504e-004; ng(n+1)=4.162755e+002;
n=4; farx(n+1)=4.357184e+002; foe(n+1)=4.401216e+002; krok(n+1)=2.322394e-004; ng(n+1)=5.472186e+002;
n=5; farx(n+1)=4.421172e+002; foe(n+1)=4.251073e+002; krok(n+1)=3.800744e-004; ng(n+1)=3.842304e+002;
n=6; farx(n+1)=3.972623e+002; foe(n+1)=4.112135e+002; krok(n+1)=2.224690e-004; ng(n+1)=4.851887e+002;
n=7; farx(n+1)=3.993365e+002; foe(n+1)=3.983340e+002; krok(n+1)=3.831728e-004; ng(n+1)=3.543108e+002;
n=8; farx(n+1)=3.608560e+002; foe(n+1)=3.862614e+002; krok(n+1)=2.154098e-004; ng(n+1)=4.504007e+002;
n=9; farx(n+1)=3.602505e+002; foe(n+1)=3.748729e+002; krok(n+1)=3.870231e-004; ng(n+1)=3.274355e+002;
n=10; farx(n+1)=3.269128e+002; foe(n+1)=3.639900e+002; krok(n+1)=2.084136e-004; ng(n+1)=4.389143e+002;
n=11; farx(n+1)=3.242883e+002; foe(n+1)=3.538528e+002; krok(n+1)=3.812452e-004; ng(n+1)=3.031655e+002;
n=12; farx(n+1)=2.962445e+002; foe(n+1)=3.442359e+002; krok(n+1)=2.045351e-004; ng(n+1)=4.133139e+002;
n=13; farx(n+1)=2.920373e+002; foe(n+1)=3.351457e+002; krok(n+1)=3.831728e-004; ng(n+1)=2.879441e+002;
n=14; farx(n+1)=2.684129e+002; foe(n+1)=3.266248e+002; krok(n+1)=2.015675e-004; ng(n+1)=3.913592e+002;
n=15; farx(n+1)=2.632200e+002; foe(n+1)=3.185577e+002; krok(n+1)=3.870864e-004; ng(n+1)=2.786131e+002;
n=16; farx(n+1)=2.432999e+002; foe(n+1)=3.109902e+002; krok(n+1)=1.985317e-004; ng(n+1)=3.696804e+002;
n=17; farx(n+1)=2.374187e+002; foe(n+1)=3.038523e+002; krok(n+1)=3.920873e-004; ng(n+1)=2.750258e+002;
n=18; farx(n+1)=2.205932e+002; foe(n+1)=2.971545e+002; krok(n+1)=1.969720e-004; ng(n+1)=3.461015e+002;
n=19; farx(n+1)=2.144701e+002; foe(n+1)=2.908992e+002; krok(n+1)=3.889679e-004; ng(n+1)=2.704165e+002;
n=20; farx(n+1)=2.003198e+002; foe(n+1)=2.850356e+002; krok(n+1)=1.966794e-004; ng(n+1)=3.251769e+002;
n=21; farx(n+1)=1.940942e+002; foe(n+1)=2.795158e+002; krok(n+1)=3.889679e-004; ng(n+1)=2.640238e+002;
n=22; farx(n+1)=1.820487e+002; foe(n+1)=2.743119e+002; krok(n+1)=1.960919e-004; ng(n+1)=3.161786e+002;
n=23; farx(n+1)=1.759480e+002; foe(n+1)=2.694259e+002; krok(n+1)=3.831728e-004; ng(n+1)=2.577506e+002;
n=24; farx(n+1)=1.656310e+002; foe(n+1)=2.647897e+002; krok(n+1)=1.963632e-004; ng(n+1)=3.087659e+002;
n=25; farx(n+1)=1.597843e+002; foe(n+1)=2.604362e+002; krok(n+1)=3.749852e-004; ng(n+1)=2.523147e+002;
n=26; farx(n+1)=1.508709e+002; foe(n+1)=2.562319e+002; krok(n+1)=1.954649e-004; ng(n+1)=3.042207e+002;
n=27; farx(n+1)=1.453071e+002; foe(n+1)=2.522959e+002; krok(n+1)=3.643336e-004; ng(n+1)=2.472796e+002;
n=28; farx(n+1)=1.376379e+002; foe(n+1)=2.484349e+002; krok(n+1)=1.915864e-004; ng(n+1)=3.009415e+002;
n=29; farx(n+1)=1.320902e+002; foe(n+1)=2.447486e+002; krok(n+1)=3.674633e-004; ng(n+1)=2.403137e+002;
n=30; farx(n+1)=1.254048e+002; foe(n+1)=2.410773e+002; krok(n+1)=1.858436e-004; ng(n+1)=3.055903e+002;
n=31; farx(n+1)=1.201329e+002; foe(n+1)=2.376479e+002; krok(n+1)=3.551278e-004; ng(n+1)=2.361701e+002;
n=32; farx(n+1)=1.143552e+002; foe(n+1)=2.341804e+002; krok(n+1)=1.803354e-004; ng(n+1)=3.071024e+002;
n=33; farx(n+1)=1.094008e+002; foe(n+1)=2.309804e+002; krok(n+1)=3.395966e-004; ng(n+1)=2.333539e+002;
n=34; farx(n+1)=1.044640e+002; foe(n+1)=2.276707e+002; krok(n+1)=1.711903e-004; ng(n+1)=3.109245e+002;
n=35; farx(n+1)=9.964105e+001; foe(n+1)=2.246071e+002; krok(n+1)=3.347436e-004; ng(n+1)=2.289752e+002;
n=36; farx(n+1)=9.539938e+001; foe(n+1)=2.213731e+002; krok(n+1)=1.614953e-004; ng(n+1)=3.210520e+002;
n=37; farx(n+1)=9.110902e+001; foe(n+1)=2.185539e+002; krok(n+1)=3.019151e-004; ng(n+1)=2.287987e+002;
n=38; farx(n+1)=8.754896e+001; foe(n+1)=2.155388e+002; krok(n+1)=1.522476e-004; ng(n+1)=3.205868e+002;
n=39; farx(n+1)=8.371655e+001; foe(n+1)=2.129411e+002; krok(n+1)=2.734728e-004; ng(n+1)=2.292874e+002;
n=40; farx(n+1)=8.073988e+001; foe(n+1)=2.101156e+002; krok(n+1)=1.410971e-004; ng(n+1)=3.227231e+002;
n=41; farx(n+1)=7.733823e+001; foe(n+1)=2.077352e+002; krok(n+1)=2.463545e-004; ng(n+1)=2.339517e+002;
n=42; farx(n+1)=7.485322e+001; foe(n+1)=2.051006e+002; krok(n+1)=1.297027e-004; ng(n+1)=3.252887e+002;
n=43; farx(n+1)=7.191755e+001; foe(n+1)=2.029675e+002; krok(n+1)=2.162775e-004; ng(n+1)=2.395595e+002;
n=44; farx(n+1)=6.985047e+001; foe(n+1)=2.005442e+002; krok(n+1)=1.180809e-004; ng(n+1)=3.260692e+002;
n=45; farx(n+1)=6.732053e+001; foe(n+1)=1.986371e+002; krok(n+1)=1.897438e-004; ng(n+1)=2.439315e+002;
n=46; farx(n+1)=6.560461e+001; foe(n+1)=1.964254e+002; krok(n+1)=1.062896e-004; ng(n+1)=3.263731e+002;
n=47; farx(n+1)=6.336953e+001; foe(n+1)=1.947078e+002; krok(n+1)=1.714180e-004; ng(n+1)=2.460682e+002;
n=48; farx(n+1)=6.193712e+001; foe(n+1)=1.926602e+002; krok(n+1)=9.501859e-005; ng(n+1)=3.308959e+002;
n=49; farx(n+1)=6.000234e+001; foe(n+1)=1.911275e+002; krok(n+1)=1.516739e-004; ng(n+1)=2.480870e+002;
n=50; farx(n+1)=5.880835e+001; foe(n+1)=1.892719e+002; krok(n+1)=8.459568e-005; ng(n+1)=3.314627e+002;
n=51; farx(n+1)=5.708262e+001; foe(n+1)=1.878868e+002; krok(n+1)=1.389271e-004; ng(n+1)=2.479118e+002;
n=52; farx(n+1)=5.607836e+001; foe(n+1)=1.861725e+002; krok(n+1)=7.528795e-005; ng(n+1)=3.362549e+002;
n=53; farx(n+1)=5.455710e+001; foe(n+1)=1.849290e+002; krok(n+1)=1.258386e-004; ng(n+1)=2.478956e+002;
n=54; farx(n+1)=5.370926e+001; foe(n+1)=1.833613e+002; krok(n+1)=6.703667e-005; ng(n+1)=3.388518e+002;
n=55; farx(n+1)=5.235102e+001; foe(n+1)=1.822327e+002; krok(n+1)=1.156138e-004; ng(n+1)=2.465857e+002;
n=56; farx(n+1)=5.162869e+001; foe(n+1)=1.807936e+002; krok(n+1)=6.000183e-005; ng(n+1)=3.417779e+002;
n=57; farx(n+1)=5.042795e+001; foe(n+1)=1.797756e+002; krok(n+1)=1.052302e-004; ng(n+1)=2.457517e+002;
n=58; farx(n+1)=4.980481e+001; foe(n+1)=1.784629e+002; krok(n+1)=5.427827e-005; ng(n+1)=3.429008e+002;
n=59; farx(n+1)=4.878340e+001; foe(n+1)=1.775649e+002; krok(n+1)=9.199795e-005; ng(n+1)=2.470943e+002;
n=60; farx(n+1)=4.824249e+001; foe(n+1)=1.763974e+002; krok(n+1)=4.924301e-005; ng(n+1)=3.385429e+002;
n=61; farx(n+1)=4.734001e+001; foe(n+1)=1.755870e+002; krok(n+1)=8.368590e-005; ng(n+1)=2.468574e+002;
n=62; farx(n+1)=4.686746e+001; foe(n+1)=1.745191e+002; krok(n+1)=4.461789e-005; ng(n+1)=3.379765e+002;
n=63; farx(n+1)=4.601563e+001; foe(n+1)=1.737503e+002; krok(n+1)=8.147519e-005; ng(n+1)=2.425776e+002;
n=64; farx(n+1)=4.559469e+001; foe(n+1)=1.727374e+002; krok(n+1)=4.116566e-005; ng(n+1)=3.427658e+002;
n=65; farx(n+1)=4.485020e+001; foe(n+1)=1.720439e+002; krok(n+1)=7.324517e-005; ng(n+1)=2.428957e+002;
n=66; farx(n+1)=4.447529e+001; foe(n+1)=1.711269e+002; krok(n+1)=3.791849e-005; ng(n+1)=3.387426e+002;
n=67; farx(n+1)=4.379254e+001; foe(n+1)=1.704814e+002; krok(n+1)=6.918420e-005; ng(n+1)=2.410517e+002;
n=68; farx(n+1)=4.345426e+001; foe(n+1)=1.696227e+002; krok(n+1)=3.527374e-005; ng(n+1)=3.392704e+002;
n=69; farx(n+1)=4.283243e+001; foe(n+1)=1.690231e+002; krok(n+1)=6.485137e-005; ng(n+1)=2.397915e+002;
n=70; farx(n+1)=4.252563e+001; foe(n+1)=1.682224e+002; krok(n+1)=3.291660e-005; ng(n+1)=3.381758e+002;
n=71; farx(n+1)=4.195161e+001; foe(n+1)=1.676572e+002; krok(n+1)=6.158863e-005; ng(n+1)=2.376962e+002;
n=72; farx(n+1)=4.167211e+001; foe(n+1)=1.669068e+002; krok(n+1)=3.079431e-005; ng(n+1)=3.370781e+002;
n=73; farx(n+1)=4.112232e+001; foe(n+1)=1.663643e+002; krok(n+1)=6.078783e-005; ng(n+1)=2.345042e+002;
n=74; farx(n+1)=4.086388e+001; foe(n+1)=1.656374e+002; krok(n+1)=2.921834e-005; ng(n+1)=3.404508e+002;
n=75; farx(n+1)=4.036622e+001; foe(n+1)=1.651338e+002; krok(n+1)=5.658127e-005; ng(n+1)=2.341878e+002;
n=76; farx(n+1)=4.012892e+001; foe(n+1)=1.644570e+002; krok(n+1)=2.746916e-005; ng(n+1)=3.368206e+002;
n=77; farx(n+1)=3.963842e+001; foe(n+1)=1.639614e+002; krok(n+1)=5.745483e-005; ng(n+1)=2.297846e+002;
n=78; farx(n+1)=3.941736e+001; foe(n+1)=1.632937e+002; krok(n+1)=2.617217e-005; ng(n+1)=3.415212e+002;
n=79; farx(n+1)=3.894839e+001; foe(n+1)=1.628161e+002; krok(n+1)=5.658127e-005; ng(n+1)=2.271499e+002;
n=80; farx(n+1)=3.874147e+001; foe(n+1)=1.621690e+002; krok(n+1)=2.505415e-005; ng(n+1)=3.428680e+002;
n=81; farx(n+1)=3.830686e+001; foe(n+1)=1.617162e+002; krok(n+1)=5.394830e-005; ng(n+1)=2.258383e+002;
n=82; farx(n+1)=3.811166e+001; foe(n+1)=1.611045e+002; krok(n+1)=2.418894e-005; ng(n+1)=3.398659e+002;
n=83; farx(n+1)=3.773229e+001; foe(n+1)=1.606925e+002; krok(n+1)=4.837788e-005; ng(n+1)=2.277965e+002;
n=84; farx(n+1)=3.754880e+001; foe(n+1)=1.601355e+002; krok(n+1)=2.322634e-005; ng(n+1)=3.307571e+002;
n=85; farx(n+1)=3.720095e+001; foe(n+1)=1.597445e+002; krok(n+1)=4.553746e-005; ng(n+1)=2.271235e+002;
n=86; farx(n+1)=3.702800e+001; foe(n+1)=1.592234e+002; krok(n+1)=2.230894e-005; ng(n+1)=3.264892e+002;
n=87; farx(n+1)=3.668826e+001; foe(n+1)=1.588389e+002; krok(n+1)=4.572173e-005; ng(n+1)=2.239884e+002;
n=88; farx(n+1)=3.652500e+001; foe(n+1)=1.583269e+002; krok(n+1)=2.139878e-005; ng(n+1)=3.281613e+002;
n=89; farx(n+1)=3.616522e+001; foe(n+1)=1.579303e+002; krok(n+1)=4.989171e-005; ng(n+1)=2.178223e+002;
n=90; farx(n+1)=3.600730e+001; foe(n+1)=1.573930e+002; krok(n+1)=2.102242e-005; ng(n+1)=3.374689e+002;
n=91; farx(n+1)=3.567601e+001; foe(n+1)=1.570173e+002; krok(n+1)=4.724951e-005; ng(n+1)=2.173250e+002;
n=92; farx(n+1)=3.552564e+001; foe(n+1)=1.565104e+002; krok(n+1)=2.036880e-005; ng(n+1)=3.335373e+002;
n=93; farx(n+1)=3.520722e+001; foe(n+1)=1.561447e+002; krok(n+1)=4.673612e-005; ng(n+1)=2.152409e+002;
n=94; farx(n+1)=3.506337e+001; foe(n+1)=1.556500e+002; krok(n+1)=1.977673e-005; ng(n+1)=3.332930e+002;
n=95; farx(n+1)=3.474633e+001; foe(n+1)=1.552837e+002; krok(n+1)=4.789661e-005; ng(n+1)=2.114594e+002;
n=96; farx(n+1)=3.460693e+001; foe(n+1)=1.547901e+002; krok(n+1)=1.944453e-005; ng(n+1)=3.353112e+002;
n=97; farx(n+1)=3.431302e+001; foe(n+1)=1.544419e+002; krok(n+1)=4.572173e-005; ng(n+1)=2.108245e+002;
n=98; farx(n+1)=3.417883e+001; foe(n+1)=1.539725e+002; krok(n+1)=1.900886e-005; ng(n+1)=3.314894e+002;
n=99; farx(n+1)=3.390032e+001; foe(n+1)=1.536371e+002; krok(n+1)=4.462420e-005; ng(n+1)=2.094823e+002;
n=100; farx(n+1)=3.377048e+001; foe(n+1)=1.531822e+002; krok(n+1)=1.864430e-005; ng(n+1)=3.295487e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
