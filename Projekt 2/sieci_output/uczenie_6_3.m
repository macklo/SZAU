%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.845045e+003; foe(n+1)=4.753288e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.793509e+003; foe(n+1)=3.634721e+003; krok(n+1)=3.831728e-004; ng(n+1)=4.310667e+003;
n=2; farx(n+1)=1.216943e+003; foe(n+1)=9.953344e+002; krok(n+1)=1.358387e-003; ng(n+1)=3.867285e+003;
n=3; farx(n+1)=1.437290e+003; foe(n+1)=7.910893e+002; krok(n+1)=2.659553e-004; ng(n+1)=4.800297e+003;
n=4; farx(n+1)=5.749311e+002; foe(n+1)=5.471117e+002; krok(n+1)=5.485375e-003; ng(n+1)=7.798172e+002;
n=5; farx(n+1)=3.568533e+002; foe(n+1)=4.721239e+002; krok(n+1)=5.440190e-004; ng(n+1)=2.058268e+003;
n=6; farx(n+1)=2.919614e+002; foe(n+1)=4.015172e+002; krok(n+1)=1.303603e-003; ng(n+1)=3.282005e+003;
n=7; farx(n+1)=2.506364e+002; foe(n+1)=3.656730e+002; krok(n+1)=1.310706e-003; ng(n+1)=2.111720e+003;
n=8; farx(n+1)=2.217406e+002; foe(n+1)=3.512301e+002; krok(n+1)=9.407062e-004; ng(n+1)=1.415805e+003;
n=9; farx(n+1)=2.150731e+002; foe(n+1)=3.444655e+002; krok(n+1)=9.917166e-004; ng(n+1)=3.824093e+002;
n=10; farx(n+1)=1.799018e+002; foe(n+1)=3.178145e+002; krok(n+1)=3.127439e-003; ng(n+1)=5.946958e+002;
n=11; farx(n+1)=1.593575e+002; foe(n+1)=3.099911e+002; krok(n+1)=3.233551e-004; ng(n+1)=2.555668e+003;
n=12; farx(n+1)=8.381950e+001; foe(n+1)=2.426528e+002; krok(n+1)=5.915528e-003; ng(n+1)=2.050005e+003;
n=13; farx(n+1)=6.770007e+001; foe(n+1)=2.314041e+002; krok(n+1)=7.239334e-005; ng(n+1)=2.037747e+003;
n=14; farx(n+1)=6.764170e+001; foe(n+1)=2.283621e+002; krok(n+1)=4.187548e-004; ng(n+1)=3.883263e+003;
n=15; farx(n+1)=6.439684e+001; foe(n+1)=2.265858e+002; krok(n+1)=1.874926e-004; ng(n+1)=4.253479e+003;
n=16; farx(n+1)=6.350281e+001; foe(n+1)=2.212924e+002; krok(n+1)=1.812955e-003; ng(n+1)=3.934138e+003;
n=17; farx(n+1)=6.634167e+001; foe(n+1)=1.908251e+002; krok(n+1)=4.974525e-003; ng(n+1)=4.181838e+003;
n=18; farx(n+1)=6.453671e+001; foe(n+1)=1.837501e+002; krok(n+1)=2.638590e-004; ng(n+1)=4.894461e+003;
n=19; farx(n+1)=6.346801e+001; foe(n+1)=1.736096e+002; krok(n+1)=4.532388e-004; ng(n+1)=5.545636e+003;
n=20; farx(n+1)=5.775160e+001; foe(n+1)=1.677068e+002; krok(n+1)=7.499704e-004; ng(n+1)=7.069307e+003;
n=21; farx(n+1)=6.028607e+001; foe(n+1)=1.457093e+002; krok(n+1)=3.409550e-003; ng(n+1)=6.054794e+003;
n=22; farx(n+1)=6.161968e+001; foe(n+1)=1.450742e+002; krok(n+1)=1.845902e-004; ng(n+1)=3.927152e+003;
n=23; farx(n+1)=6.463343e+001; foe(n+1)=1.437512e+002; krok(n+1)=7.648589e-004; ng(n+1)=2.987705e+003;
n=24; farx(n+1)=6.831031e+001; foe(n+1)=1.407973e+002; krok(n+1)=2.129131e-003; ng(n+1)=2.545704e+003;
n=25; farx(n+1)=6.648550e+001; foe(n+1)=1.398258e+002; krok(n+1)=2.450413e-003; ng(n+1)=2.556205e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=6.645029e+001; foe(n+1)=1.386596e+002; krok(n+1)=1.521735e-006; ng(n+1)=2.667773e+003;
n=27; farx(n+1)=6.328796e+001; foe(n+1)=1.326059e+002; krok(n+1)=9.456357e-006; ng(n+1)=3.065889e+003;
n=28; farx(n+1)=6.108783e+001; foe(n+1)=1.275193e+002; krok(n+1)=7.624450e-005; ng(n+1)=1.225366e+003;
n=29; farx(n+1)=6.090928e+001; foe(n+1)=1.004674e+002; krok(n+1)=1.775310e-004; ng(n+1)=1.708743e+003;
n=30; farx(n+1)=5.915251e+001; foe(n+1)=9.294809e+001; krok(n+1)=1.297027e-004; ng(n+1)=5.581835e+003;
n=31; farx(n+1)=4.721516e+001; foe(n+1)=7.957459e+001; krok(n+1)=6.305647e-004; ng(n+1)=6.500968e+003;
n=32; farx(n+1)=2.840323e+001; foe(n+1)=5.611062e+001; krok(n+1)=1.146884e-003; ng(n+1)=3.889502e+003;
n=33; farx(n+1)=2.447112e+001; foe(n+1)=4.818291e+001; krok(n+1)=6.298318e-004; ng(n+1)=8.126825e+003;
n=34; farx(n+1)=2.182349e+001; foe(n+1)=3.734190e+001; krok(n+1)=5.433546e-003; ng(n+1)=1.858699e+003;
n=35; farx(n+1)=1.031054e+001; foe(n+1)=2.584588e+001; krok(n+1)=2.123614e-003; ng(n+1)=3.380096e+003;
n=36; farx(n+1)=7.877094e+000; foe(n+1)=1.769839e+001; krok(n+1)=1.387045e-003; ng(n+1)=2.905701e+003;
n=37; farx(n+1)=5.706806e+000; foe(n+1)=1.367256e+001; krok(n+1)=1.220985e-003; ng(n+1)=1.523452e+003;
n=38; farx(n+1)=5.211917e+000; foe(n+1)=1.265262e+001; krok(n+1)=2.409214e-003; ng(n+1)=7.388518e+002;
n=39; farx(n+1)=4.455465e+000; foe(n+1)=1.131633e+001; krok(n+1)=8.704305e-003; ng(n+1)=1.388618e+003;
n=40; farx(n+1)=4.236498e+000; foe(n+1)=1.065904e+001; krok(n+1)=8.085176e-003; ng(n+1)=5.740726e+002;
n=41; farx(n+1)=4.168188e+000; foe(n+1)=1.035002e+001; krok(n+1)=8.300976e-003; ng(n+1)=3.363316e+002;
n=42; farx(n+1)=4.227851e+000; foe(n+1)=1.020690e+001; krok(n+1)=1.297031e-002; ng(n+1)=6.020979e+002;
n=43; farx(n+1)=4.154628e+000; foe(n+1)=1.012119e+001; krok(n+1)=9.949946e-003; ng(n+1)=5.584028e+002;
n=44; farx(n+1)=3.940471e+000; foe(n+1)=9.905700e+000; krok(n+1)=4.486659e-002; ng(n+1)=4.769972e+002;
n=45; farx(n+1)=3.904963e+000; foe(n+1)=9.752926e+000; krok(n+1)=1.464382e-002; ng(n+1)=8.529490e+002;
n=46; farx(n+1)=3.620540e+000; foe(n+1)=9.497469e+000; krok(n+1)=5.529705e-002; ng(n+1)=6.871998e+002;
n=47; farx(n+1)=3.173475e+000; foe(n+1)=8.779264e+000; krok(n+1)=4.821088e-002; ng(n+1)=7.953849e+002;
n=48; farx(n+1)=3.248725e+000; foe(n+1)=8.418800e+000; krok(n+1)=6.266827e-002; ng(n+1)=7.496523e+002;
n=49; farx(n+1)=3.062512e+000; foe(n+1)=7.940567e+000; krok(n+1)=1.943277e-001; ng(n+1)=4.259080e+002;
n=50; farx(n+1)=3.041292e+000; foe(n+1)=7.417965e+000; krok(n+1)=4.568872e-002; ng(n+1)=1.345744e+003;
%odnowa zmiennej metryki
n=51; farx(n+1)=3.032676e+000; foe(n+1)=7.326408e+000; krok(n+1)=1.226134e-006; ng(n+1)=8.311536e+002;
n=52; farx(n+1)=3.022018e+000; foe(n+1)=7.274242e+000; krok(n+1)=1.636503e-005; ng(n+1)=2.302368e+002;
n=53; farx(n+1)=2.984369e+000; foe(n+1)=7.230584e+000; krok(n+1)=2.885108e-005; ng(n+1)=1.682727e+002;
n=54; farx(n+1)=2.945589e+000; foe(n+1)=7.106330e+000; krok(n+1)=3.226735e-004; ng(n+1)=8.542424e+001;
n=55; farx(n+1)=2.906017e+000; foe(n+1)=6.993431e+000; krok(n+1)=2.852289e-004; ng(n+1)=9.783413e+001;
n=56; farx(n+1)=2.830496e+000; foe(n+1)=6.657031e+000; krok(n+1)=1.308315e-003; ng(n+1)=8.700036e+001;
n=57; farx(n+1)=2.748436e+000; foe(n+1)=6.452300e+000; krok(n+1)=4.752741e-003; ng(n+1)=1.579626e+002;
n=58; farx(n+1)=2.611037e+000; foe(n+1)=6.118598e+000; krok(n+1)=3.576264e-003; ng(n+1)=1.362391e+002;
n=59; farx(n+1)=2.613319e+000; foe(n+1)=5.897395e+000; krok(n+1)=6.727259e-003; ng(n+1)=8.258256e+002;
n=60; farx(n+1)=2.505993e+000; foe(n+1)=5.712455e+000; krok(n+1)=5.767896e-003; ng(n+1)=6.372786e+002;
n=61; farx(n+1)=2.395732e+000; foe(n+1)=5.460653e+000; krok(n+1)=2.113019e-002; ng(n+1)=9.267838e+002;
n=62; farx(n+1)=2.366014e+000; foe(n+1)=5.320878e+000; krok(n+1)=2.074719e-002; ng(n+1)=1.104777e+003;
n=63; farx(n+1)=2.346545e+000; foe(n+1)=5.270351e+000; krok(n+1)=1.405854e-002; ng(n+1)=6.817231e+002;
n=64; farx(n+1)=2.316303e+000; foe(n+1)=5.247094e+000; krok(n+1)=8.799720e-003; ng(n+1)=3.966776e+002;
n=65; farx(n+1)=2.122923e+000; foe(n+1)=5.154701e+000; krok(n+1)=1.785858e-002; ng(n+1)=3.882326e+002;
n=66; farx(n+1)=1.878871e+000; foe(n+1)=4.943517e+000; krok(n+1)=2.254004e-002; ng(n+1)=3.223555e+002;
n=67; farx(n+1)=1.680275e+000; foe(n+1)=4.579771e+000; krok(n+1)=5.558094e-002; ng(n+1)=6.863926e+002;
n=68; farx(n+1)=1.355834e+000; foe(n+1)=4.339283e+000; krok(n+1)=2.954351e-002; ng(n+1)=1.089663e+003;
n=69; farx(n+1)=1.219492e+000; foe(n+1)=4.203670e+000; krok(n+1)=1.694661e-002; ng(n+1)=1.206578e+003;
n=70; farx(n+1)=1.179413e+000; foe(n+1)=3.983526e+000; krok(n+1)=3.571716e-002; ng(n+1)=1.135021e+003;
n=71; farx(n+1)=1.150360e+000; foe(n+1)=3.801637e+000; krok(n+1)=2.243330e-002; ng(n+1)=5.986976e+002;
n=72; farx(n+1)=1.132922e+000; foe(n+1)=3.597090e+000; krok(n+1)=5.120727e-002; ng(n+1)=9.834454e+002;
n=73; farx(n+1)=1.152923e+000; foe(n+1)=3.243867e+000; krok(n+1)=3.618094e-002; ng(n+1)=7.484281e+002;
n=74; farx(n+1)=1.141642e+000; foe(n+1)=2.951370e+000; krok(n+1)=7.727538e-002; ng(n+1)=7.901427e+002;
n=75; farx(n+1)=1.031536e+000; foe(n+1)=2.791988e+000; krok(n+1)=6.078596e-002; ng(n+1)=5.097624e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.029408e+000; foe(n+1)=2.766451e+000; krok(n+1)=3.792618e-006; ng(n+1)=2.794391e+002;
n=77; farx(n+1)=1.030267e+000; foe(n+1)=2.706734e+000; krok(n+1)=1.615744e-006; ng(n+1)=5.620063e+002;
n=78; farx(n+1)=1.030532e+000; foe(n+1)=2.701790e+000; krok(n+1)=3.917546e-006; ng(n+1)=1.369040e+002;
n=79; farx(n+1)=1.028058e+000; foe(n+1)=2.687022e+000; krok(n+1)=2.681467e-004; ng(n+1)=2.559690e+001;
n=80; farx(n+1)=1.026045e+000; foe(n+1)=2.681216e+000; krok(n+1)=1.205665e-004; ng(n+1)=3.105770e+001;
n=81; farx(n+1)=1.012324e+000; foe(n+1)=2.664607e+000; krok(n+1)=8.336543e-004; ng(n+1)=2.342160e+001;
n=82; farx(n+1)=9.692212e-001; foe(n+1)=2.638338e+000; krok(n+1)=2.827577e-003; ng(n+1)=1.519506e+001;
n=83; farx(n+1)=9.192754e-001; foe(n+1)=2.593209e+000; krok(n+1)=3.858127e-003; ng(n+1)=1.665609e+001;
n=84; farx(n+1)=9.094664e-001; foe(n+1)=2.575770e+000; krok(n+1)=4.290347e-003; ng(n+1)=3.904293e+001;
n=85; farx(n+1)=8.895597e-001; foe(n+1)=2.539450e+000; krok(n+1)=1.111520e-002; ng(n+1)=7.765762e+001;
n=86; farx(n+1)=8.986437e-001; foe(n+1)=2.522628e+000; krok(n+1)=4.375565e-003; ng(n+1)=2.453306e+002;
n=87; farx(n+1)=8.953127e-001; foe(n+1)=2.467770e+000; krok(n+1)=1.384117e-002; ng(n+1)=3.286322e+002;
n=88; farx(n+1)=8.991068e-001; foe(n+1)=2.434426e+000; krok(n+1)=2.550930e-002; ng(n+1)=4.265726e+002;
n=89; farx(n+1)=8.822641e-001; foe(n+1)=2.411948e+000; krok(n+1)=1.111520e-002; ng(n+1)=4.479571e+002;
n=90; farx(n+1)=8.645310e-001; foe(n+1)=2.387383e+000; krok(n+1)=2.223039e-002; ng(n+1)=2.417609e+002;
n=91; farx(n+1)=8.296322e-001; foe(n+1)=2.347136e+000; krok(n+1)=4.000365e-002; ng(n+1)=5.021560e+002;
n=92; farx(n+1)=7.302977e-001; foe(n+1)=2.253697e+000; krok(n+1)=3.221469e-002; ng(n+1)=5.103072e+002;
n=93; farx(n+1)=6.876264e-001; foe(n+1)=2.201938e+000; krok(n+1)=5.142093e-002; ng(n+1)=4.966886e+002;
n=94; farx(n+1)=6.225043e-001; foe(n+1)=2.156782e+000; krok(n+1)=3.221469e-002; ng(n+1)=3.190672e+002;
n=95; farx(n+1)=5.888432e-001; foe(n+1)=2.120415e+000; krok(n+1)=1.975828e-002; ng(n+1)=2.546008e+002;
n=96; farx(n+1)=5.896479e-001; foe(n+1)=2.067854e+000; krok(n+1)=3.482091e-002; ng(n+1)=6.794680e+002;
n=97; farx(n+1)=5.990091e-001; foe(n+1)=2.004429e+000; krok(n+1)=1.489153e-001; ng(n+1)=2.904880e+002;
n=98; farx(n+1)=5.677408e-001; foe(n+1)=1.956705e+000; krok(n+1)=4.446078e-002; ng(n+1)=7.027487e+002;
n=99; farx(n+1)=5.984243e-001; foe(n+1)=1.787501e+000; krok(n+1)=3.238984e-001; ng(n+1)=2.085780e+002;
n=100; farx(n+1)=5.966585e-001; foe(n+1)=1.734923e+000; krok(n+1)=8.662598e-002; ng(n+1)=6.008563e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
