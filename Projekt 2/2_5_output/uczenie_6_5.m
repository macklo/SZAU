%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.283095e+003; foe(n+1)=4.382945e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
n=1; farx(n+1)=2.840131e+003; foe(n+1)=2.969108e+003; krok(n+1)=4.884852e-004; ng(n+1)=2.159885e+003;
n=2; farx(n+1)=5.050619e+002; foe(n+1)=5.170515e+002; krok(n+1)=5.433546e-003; ng(n+1)=6.907630e+002;
n=3; farx(n+1)=5.200278e+002; foe(n+1)=4.249566e+002; krok(n+1)=1.243743e-003; ng(n+1)=5.662658e+002;
n=4; farx(n+1)=4.351834e+002; foe(n+1)=4.003464e+002; krok(n+1)=2.427276e-004; ng(n+1)=7.114474e+002;
n=5; farx(n+1)=4.165824e+002; foe(n+1)=3.843884e+002; krok(n+1)=6.459813e-004; ng(n+1)=2.160698e+002;
n=6; farx(n+1)=3.754228e+002; foe(n+1)=3.714014e+002; krok(n+1)=2.619761e-004; ng(n+1)=4.377612e+002;
n=7; farx(n+1)=3.677343e+002; foe(n+1)=3.594074e+002; krok(n+1)=6.694872e-004; ng(n+1)=2.341696e+002;
n=8; farx(n+1)=3.347894e+002; foe(n+1)=3.484198e+002; krok(n+1)=2.574102e-004; ng(n+1)=4.037424e+002;
n=9; farx(n+1)=3.256816e+002; foe(n+1)=3.384987e+002; krok(n+1)=6.935224e-004; ng(n+1)=2.127268e+002;
n=10; farx(n+1)=3.003823e+002; foe(n+1)=3.293988e+002; krok(n+1)=2.559644e-004; ng(n+1)=3.610313e+002;
n=11; farx(n+1)=2.908165e+002; foe(n+1)=3.212311e+002; krok(n+1)=7.138862e-004; ng(n+1)=2.029856e+002;
n=12; farx(n+1)=2.716966e+002; foe(n+1)=3.140045e+002; krok(n+1)=2.609151e-004; ng(n+1)=3.142222e+002;
n=13; farx(n+1)=2.628305e+002; foe(n+1)=3.076047e+002; krok(n+1)=7.125145e-004; ng(n+1)=1.934349e+002;
n=14; farx(n+1)=2.482173e+002; foe(n+1)=3.018597e+002; krok(n+1)=2.661414e-004; ng(n+1)=2.805817e+002;
n=15; farx(n+1)=2.400345e+002; foe(n+1)=2.967843e+002; krok(n+1)=7.101242e-004; ng(n+1)=1.826469e+002;
n=16; farx(n+1)=2.286684e+002; foe(n+1)=2.921641e+002; krok(n+1)=2.702675e-004; ng(n+1)=2.528653e+002;
n=17; farx(n+1)=2.209830e+002; foe(n+1)=2.880110e+002; krok(n+1)=7.102557e-004; ng(n+1)=1.719617e+002;
n=18; farx(n+1)=2.118616e+002; foe(n+1)=2.841888e+002; krok(n+1)=2.752027e-004; ng(n+1)=2.311628e+002;
n=19; farx(n+1)=2.049543e+002; foe(n+1)=2.807747e+002; krok(n+1)=6.767654e-004; ng(n+1)=1.648351e+002;
n=20; farx(n+1)=1.975004e+002; foe(n+1)=2.775561e+002; krok(n+1)=2.790064e-004; ng(n+1)=2.192365e+002;
n=21; farx(n+1)=1.910739e+002; foe(n+1)=2.746037e+002; krok(n+1)=6.500859e-004; ng(n+1)=1.595354e+002;
n=22; farx(n+1)=1.847979e+002; foe(n+1)=2.717404e+002; krok(n+1)=2.778542e-004; ng(n+1)=2.145353e+002;
n=23; farx(n+1)=1.786734e+002; foe(n+1)=2.690748e+002; krok(n+1)=6.250570e-004; ng(n+1)=1.560202e+002;
n=24; farx(n+1)=1.732993e+002; foe(n+1)=2.664207e+002; krok(n+1)=2.714211e-004; ng(n+1)=2.145768e+002;
n=25; farx(n+1)=1.673885e+002; foe(n+1)=2.639264e+002; krok(n+1)=5.970040e-004; ng(n+1)=1.538862e+002;
n=26; farx(n+1)=1.627422e+002; foe(n+1)=2.614147e+002; krok(n+1)=2.618406e-004; ng(n+1)=2.167879e+002;
n=27; farx(n+1)=1.570769e+002; foe(n+1)=2.590556e+002; krok(n+1)=5.617933e-004; ng(n+1)=1.536188e+002;
n=28; farx(n+1)=1.530471e+002; foe(n+1)=2.566420e+002; krok(n+1)=2.478682e-004; ng(n+1)=2.212190e+002;
n=29; farx(n+1)=1.476167e+002; foe(n+1)=2.543836e+002; krok(n+1)=5.252174e-004; ng(n+1)=1.541305e+002;
n=30; farx(n+1)=1.441203e+002; foe(n+1)=2.520562e+002; krok(n+1)=2.312276e-004; ng(n+1)=2.267746e+002;
n=31; farx(n+1)=1.389863e+002; foe(n+1)=2.499065e+002; krok(n+1)=4.838390e-004; ng(n+1)=1.585928e+002;
n=32; farx(n+1)=1.359615e+002; foe(n+1)=2.476700e+002; krok(n+1)=2.125792e-004; ng(n+1)=2.327794e+002;
n=33; farx(n+1)=1.311759e+002; foe(n+1)=2.456365e+002; krok(n+1)=4.395065e-004; ng(n+1)=1.628122e+002;
n=34; farx(n+1)=1.285608e+002; foe(n+1)=2.435135e+002; krok(n+1)=1.942746e-004; ng(n+1)=2.379326e+002;
n=35; farx(n+1)=1.241965e+002; foe(n+1)=2.416368e+002; krok(n+1)=3.920873e-004; ng(n+1)=1.676871e+002;
n=36; farx(n+1)=1.219602e+002; foe(n+1)=2.396303e+002; krok(n+1)=1.736749e-004; ng(n+1)=2.433606e+002;
n=37; farx(n+1)=1.177521e+002; foe(n+1)=2.378353e+002; krok(n+1)=3.716215e-004; ng(n+1)=1.691771e+002;
n=38; farx(n+1)=1.158402e+002; foe(n+1)=2.358870e+002; krok(n+1)=1.540602e-004; ng(n+1)=2.546488e+002;
n=39; farx(n+1)=1.120106e+002; foe(n+1)=2.342320e+002; krok(n+1)=3.325065e-004; ng(n+1)=1.720847e+002;
n=40; farx(n+1)=1.103859e+002; foe(n+1)=2.324177e+002; krok(n+1)=1.367364e-004; ng(n+1)=2.598507e+002;
n=41; farx(n+1)=1.068684e+002; foe(n+1)=2.308933e+002; krok(n+1)=3.018569e-004; ng(n+1)=1.742089e+002;
n=42; farx(n+1)=1.054886e+002; foe(n+1)=2.291866e+002; krok(n+1)=1.209232e-004; ng(n+1)=2.671875e+002;
n=43; farx(n+1)=1.023008e+002; foe(n+1)=2.277924e+002; krok(n+1)=2.710039e-004; ng(n+1)=1.761171e+002;
n=44; farx(n+1)=1.011246e+002; foe(n+1)=2.262101e+002; krok(n+1)=1.075741e-004; ng(n+1)=2.725510e+002;
n=45; farx(n+1)=9.834256e+001; foe(n+1)=2.249711e+002; krok(n+1)=2.345612e-004; ng(n+1)=1.791342e+002;
n=46; farx(n+1)=9.733634e+001; foe(n+1)=2.235406e+002; krok(n+1)=9.580102e-005; ng(n+1)=2.737837e+002;
n=47; farx(n+1)=9.487501e+001; foe(n+1)=2.224196e+002; krok(n+1)=2.064885e-004; ng(n+1)=1.808503e+002;
n=48; farx(n+1)=9.399961e+001; foe(n+1)=2.211266e+002; krok(n+1)=8.669030e-005; ng(n+1)=2.747276e+002;
n=49; farx(n+1)=9.196149e+001; foe(n+1)=2.201681e+002; krok(n+1)=1.697983e-004; ng(n+1)=1.858805e+002;
n=50; farx(n+1)=9.119790e+001; foe(n+1)=2.190359e+002; krok(n+1)=7.773395e-005; ng(n+1)=2.695791e+002;
n=51; farx(n+1)=8.927941e+001; foe(n+1)=2.181331e+002; krok(n+1)=1.602589e-004; ng(n+1)=1.850740e+002;
n=52; farx(n+1)=8.861215e+001; foe(n+1)=2.170632e+002; krok(n+1)=7.022416e-005; ng(n+1)=2.758574e+002;
n=53; farx(n+1)=8.689299e+001; foe(n+1)=2.162422e+002; krok(n+1)=1.436854e-004; ng(n+1)=1.863904e+002;
n=54; farx(n+1)=8.630595e+001; foe(n+1)=2.152617e+002; krok(n+1)=6.352069e-005; ng(n+1)=2.768453e+002;
n=55; farx(n+1)=8.470826e+001; foe(n+1)=2.144951e+002; krok(n+1)=1.340573e-004; ng(n+1)=1.861937e+002;
n=56; farx(n+1)=8.419196e+001; foe(n+1)=2.135744e+002; krok(n+1)=5.755042e-005; ng(n+1)=2.810231e+002;
n=57; farx(n+1)=8.268644e+001; foe(n+1)=2.128493e+002; krok(n+1)=1.270414e-004; ng(n+1)=1.853358e+002;
n=58; farx(n+1)=8.222850e+001; foe(n+1)=2.119778e+002; krok(n+1)=5.270595e-005; ng(n+1)=2.859816e+002;
n=59; farx(n+1)=8.087307e+001; foe(n+1)=2.113154e+002; krok(n+1)=1.147197e-004; ng(n+1)=1.865125e+002;
n=60; farx(n+1)=8.046652e+001; foe(n+1)=2.105135e+002; krok(n+1)=4.789661e-005; ng(n+1)=2.860358e+002;
n=61; farx(n+1)=7.912037e+001; foe(n+1)=2.098621e+002; krok(n+1)=1.151008e-004; ng(n+1)=1.838207e+002;
n=62; farx(n+1)=7.875977e+001; foe(n+1)=2.090722e+002; krok(n+1)=4.417967e-005; ng(n+1)=2.957806e+002;
n=63; farx(n+1)=7.752958e+001; foe(n+1)=2.084708e+002; krok(n+1)=1.057004e-004; ng(n+1)=1.844915e+002;
n=64; farx(n+1)=7.720569e+001; foe(n+1)=2.077354e+002; krok(n+1)=4.073759e-005; ng(n+1)=2.967917e+002;
n=65; farx(n+1)=7.605144e+001; foe(n+1)=2.071675e+002; krok(n+1)=9.985203e-005; ng(n+1)=1.841681e+002;
n=66; farx(n+1)=7.575940e+001; foe(n+1)=2.064706e+002; krok(n+1)=3.782543e-005; ng(n+1)=2.997780e+002;
n=67; farx(n+1)=7.469172e+001; foe(n+1)=2.059397e+002; krok(n+1)=9.290537e-005; ng(n+1)=1.844719e+002;
n=68; farx(n+1)=7.442650e+001; foe(n+1)=2.052857e+002; krok(n+1)=3.523097e-005; ng(n+1)=3.007529e+002;
n=69; farx(n+1)=7.342602e+001; foe(n+1)=2.047853e+002; krok(n+1)=8.763006e-005; ng(n+1)=1.845091e+002;
n=70; farx(n+1)=7.318435e+001; foe(n+1)=2.041644e+002; krok(n+1)=3.298006e-005; ng(n+1)=3.028440e+002;
n=71; farx(n+1)=7.225955e+001; foe(n+1)=2.036928e+002; krok(n+1)=8.147699e-005; ng(n+1)=1.846993e+002;
n=72; farx(n+1)=7.203825e+001; foe(n+1)=2.031135e+002; krok(n+1)=3.079431e-005; ng(n+1)=3.018644e+002;
n=73; farx(n+1)=7.112981e+001; foe(n+1)=2.026518e+002; krok(n+1)=8.082043e-005; ng(n+1)=1.831466e+002;
n=74; farx(n+1)=7.092879e+001; foe(n+1)=2.020820e+002; krok(n+1)=2.897758e-005; ng(n+1)=3.081091e+002;
n=75; farx(n+1)=7.006201e+001; foe(n+1)=2.016373e+002; krok(n+1)=7.773395e-005; ng(n+1)=1.824680e+002;
n=76; farx(n+1)=6.987726e+001; foe(n+1)=2.010910e+002; krok(n+1)=2.737449e-005; ng(n+1)=3.103225e+002;
n=77; farx(n+1)=6.906163e+001; foe(n+1)=2.006692e+002; krok(n+1)=7.363092e-005; ng(n+1)=1.826236e+002;
n=78; farx(n+1)=6.889013e+001; foe(n+1)=2.001485e+002; krok(n+1)=2.605170e-005; ng(n+1)=3.112665e+002;
n=79; farx(n+1)=6.815512e+001; foe(n+1)=1.997588e+002; krok(n+1)=6.657663e-005; ng(n+1)=1.841809e+002;
n=80; farx(n+1)=6.799319e+001; foe(n+1)=1.992808e+002; krok(n+1)=2.480406e-005; ng(n+1)=3.061762e+002;
n=81; farx(n+1)=6.731624e+001; foe(n+1)=1.989141e+002; krok(n+1)=6.158863e-005; ng(n+1)=1.850767e+002;
n=82; farx(n+1)=6.716350e+001; foe(n+1)=1.984679e+002; krok(n+1)=2.364703e-005; ng(n+1)=3.028213e+002;
n=83; farx(n+1)=6.651101e+001; foe(n+1)=1.981143e+002; krok(n+1)=5.977404e-005; ng(n+1)=1.848645e+002;
n=84; farx(n+1)=6.636829e+001; foe(n+1)=1.976797e+002; krok(n+1)=2.266263e-005; ng(n+1)=3.052918e+002;
n=85; farx(n+1)=6.574749e+001; foe(n+1)=1.973410e+002; krok(n+1)=5.720244e-005; ng(n+1)=1.850947e+002;
n=86; farx(n+1)=6.561232e+001; foe(n+1)=1.969227e+002; krok(n+1)=2.190752e-005; ng(n+1)=3.057278e+002;
n=87; farx(n+1)=6.505700e+001; foe(n+1)=1.966094e+002; krok(n+1)=5.121654e-005; ng(n+1)=1.873813e+002;
n=88; farx(n+1)=6.492618e+001; foe(n+1)=1.962275e+002; krok(n+1)=2.109655e-005; ng(n+1)=2.977021e+002;
n=89; farx(n+1)=6.439751e+001; foe(n+1)=1.959257e+002; krok(n+1)=4.901091e-005; ng(n+1)=1.876348e+002;
n=90; farx(n+1)=6.427313e+001; foe(n+1)=1.955588e+002; krok(n+1)=2.034152e-005; ng(n+1)=2.969363e+002;
n=91; farx(n+1)=6.376223e+001; foe(n+1)=1.952649e+002; krok(n+1)=4.766085e-005; ng(n+1)=1.873451e+002;
n=92; farx(n+1)=6.364408e+001; foe(n+1)=1.949084e+002; krok(n+1)=1.967248e-005; ng(n+1)=2.976377e+002;
n=93; farx(n+1)=6.315362e+001; foe(n+1)=1.946239e+002; krok(n+1)=4.599898e-005; ng(n+1)=1.874692e+002;
n=94; farx(n+1)=6.304102e+001; foe(n+1)=1.942786e+002; krok(n+1)=1.904168e-005; ng(n+1)=2.975469e+002;
n=95; farx(n+1)=6.256339e+001; foe(n+1)=1.940018e+002; krok(n+1)=4.506169e-005; ng(n+1)=1.873679e+002;
n=96; farx(n+1)=6.245671e+001; foe(n+1)=1.936624e+002; krok(n+1)=1.843074e-005; ng(n+1)=2.995376e+002;
n=97; farx(n+1)=6.199366e+001; foe(n+1)=1.933909e+002; krok(n+1)=4.396667e-005; ng(n+1)=1.868869e+002;
n=98; farx(n+1)=6.189157e+001; foe(n+1)=1.930612e+002; krok(n+1)=1.791657e-005; ng(n+1)=2.995608e+002;
n=99; farx(n+1)=6.144655e+001; foe(n+1)=1.927988e+002; krok(n+1)=4.244958e-005; ng(n+1)=1.873061e+002;
n=100; farx(n+1)=6.134867e+001; foe(n+1)=1.924783e+002; krok(n+1)=1.743790e-005; ng(n+1)=2.994408e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
