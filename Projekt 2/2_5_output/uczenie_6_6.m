%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.526765e+003; foe(n+1)=4.446337e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
n=1; farx(n+1)=3.036785e+003; foe(n+1)=3.046966e+003; krok(n+1)=4.957365e-004; ng(n+1)=2.329321e+003;
n=2; farx(n+1)=3.076822e+003; foe(n+1)=3.040326e+003; krok(n+1)=8.881378e-004; ng(n+1)=1.386040e+002;
n=3; farx(n+1)=6.661307e+002; foe(n+1)=6.216137e+002; krok(n+1)=3.300845e-002; ng(n+1)=1.805127e+002;
n=4; farx(n+1)=5.126071e+002; foe(n+1)=4.069593e+002; krok(n+1)=1.389400e-003; ng(n+1)=8.808301e+002;
n=5; farx(n+1)=4.842324e+002; foe(n+1)=3.717958e+002; krok(n+1)=6.767654e-004; ng(n+1)=3.049970e+002;
n=6; farx(n+1)=4.614885e+002; foe(n+1)=3.575071e+002; krok(n+1)=3.367365e-003; ng(n+1)=9.363508e+001;
n=7; farx(n+1)=4.337148e+002; foe(n+1)=3.500357e+002; krok(n+1)=7.312079e-004; ng(n+1)=2.000220e+002;
n=8; farx(n+1)=4.074390e+002; foe(n+1)=3.449566e+002; krok(n+1)=2.394037e-003; ng(n+1)=9.716583e+001;
n=9; farx(n+1)=3.967536e+002; foe(n+1)=3.408083e+002; krok(n+1)=7.102557e-004; ng(n+1)=1.585370e+002;
n=10; farx(n+1)=3.824038e+002; foe(n+1)=3.373379e+002; krok(n+1)=2.282102e-003; ng(n+1)=8.333820e+001;
n=11; farx(n+1)=3.755560e+002; foe(n+1)=3.344771e+002; krok(n+1)=6.935224e-004; ng(n+1)=1.348836e+002;
n=12; farx(n+1)=3.671613e+002; foe(n+1)=3.321690e+002; krok(n+1)=2.063028e-003; ng(n+1)=7.399136e+001;
n=13; farx(n+1)=3.625549e+002; foe(n+1)=3.302455e+002; krok(n+1)=6.935224e-004; ng(n+1)=1.139169e+002;
n=14; farx(n+1)=3.573661e+002; foe(n+1)=3.286800e+002; krok(n+1)=1.854565e-003; ng(n+1)=6.589503e+001;
n=15; farx(n+1)=3.540662e+002; foe(n+1)=3.273403e+002; krok(n+1)=6.847610e-004; ng(n+1)=9.666496e+001;
n=16; farx(n+1)=3.502503e+002; foe(n+1)=3.261432e+002; krok(n+1)=1.895132e-003; ng(n+1)=5.778431e+001;
n=17; farx(n+1)=3.476909e+002; foe(n+1)=3.251089e+002; krok(n+1)=6.696507e-004; ng(n+1)=8.568442e+001;
n=18; farx(n+1)=3.449199e+002; foe(n+1)=3.242021e+002; krok(n+1)=1.860699e-003; ng(n+1)=5.109682e+001;
n=19; farx(n+1)=3.429030e+002; foe(n+1)=3.234053e+002; krok(n+1)=6.621248e-004; ng(n+1)=7.538766e+001;
n=20; farx(n+1)=3.408209e+002; foe(n+1)=3.226957e+002; krok(n+1)=1.835515e-003; ng(n+1)=4.546677e+001;
n=21; farx(n+1)=3.391798e+002; foe(n+1)=3.220609e+002; krok(n+1)=6.586506e-004; ng(n+1)=6.702131e+001;
n=22; farx(n+1)=3.376065e+002; foe(n+1)=3.214943e+002; krok(n+1)=1.776276e-003; ng(n+1)=4.094566e+001;
n=23; farx(n+1)=3.362448e+002; foe(n+1)=3.209783e+002; krok(n+1)=6.586506e-004; ng(n+1)=5.982794e+001;
n=24; farx(n+1)=3.350169e+002; foe(n+1)=3.205080e+002; krok(n+1)=1.736904e-003; ng(n+1)=3.721166e+001;
n=25; farx(n+1)=3.338539e+002; foe(n+1)=3.200713e+002; krok(n+1)=6.586506e-004; ng(n+1)=5.423637e+001;
n=26; farx(n+1)=3.328844e+002; foe(n+1)=3.196693e+002; krok(n+1)=1.683682e-003; ng(n+1)=3.415974e+001;
n=27; farx(n+1)=3.318766e+002; foe(n+1)=3.192927e+002; krok(n+1)=6.603184e-004; ng(n+1)=4.940530e+001;
n=28; farx(n+1)=3.310895e+002; foe(n+1)=3.189368e+002; krok(n+1)=1.651908e-003; ng(n+1)=3.176707e+001;
n=29; farx(n+1)=3.301970e+002; foe(n+1)=3.186011e+002; krok(n+1)=6.603743e-004; ng(n+1)=4.596479e+001;
n=30; farx(n+1)=3.295425e+002; foe(n+1)=3.182777e+002; krok(n+1)=1.651654e-003; ng(n+1)=3.049278e+001;
n=31; farx(n+1)=3.287295e+002; foe(n+1)=3.179664e+002; krok(n+1)=6.586506e-004; ng(n+1)=4.448227e+001;
n=32; farx(n+1)=3.281937e+002; foe(n+1)=3.176696e+002; krok(n+1)=1.606327e-003; ng(n+1)=2.957065e+001;
n=33; farx(n+1)=3.274552e+002; foe(n+1)=3.173830e+002; krok(n+1)=6.603743e-004; ng(n+1)=4.248093e+001;
n=34; farx(n+1)=3.270024e+002; foe(n+1)=3.171031e+002; krok(n+1)=1.602367e-003; ng(n+1)=2.872031e+001;
n=35; farx(n+1)=3.263169e+002; foe(n+1)=3.168303e+002; krok(n+1)=6.553529e-004; ng(n+1)=4.144089e+001;
n=36; farx(n+1)=3.259237e+002; foe(n+1)=3.165598e+002; krok(n+1)=1.629095e-003; ng(n+1)=2.788747e+001;
n=37; farx(n+1)=3.252751e+002; foe(n+1)=3.162943e+002; krok(n+1)=6.502194e-004; ng(n+1)=4.095882e+001;
n=38; farx(n+1)=3.249388e+002; foe(n+1)=3.160326e+002; krok(n+1)=1.629095e-003; ng(n+1)=2.729555e+001;
n=39; farx(n+1)=3.243239e+002; foe(n+1)=3.157746e+002; krok(n+1)=6.467102e-004; ng(n+1)=4.028074e+001;
n=40; farx(n+1)=3.240328e+002; foe(n+1)=3.155181e+002; krok(n+1)=1.635770e-003; ng(n+1)=2.676724e+001;
n=41; farx(n+1)=3.234429e+002; foe(n+1)=3.152655e+002; krok(n+1)=6.465635e-004; ng(n+1)=4.022165e+001;
n=42; farx(n+1)=3.231979e+002; foe(n+1)=3.150161e+002; krok(n+1)=1.620122e-003; ng(n+1)=2.651190e+001;
n=43; farx(n+1)=3.226282e+002; foe(n+1)=3.147673e+002; krok(n+1)=6.459813e-004; ng(n+1)=4.024860e+001;
n=44; farx(n+1)=3.224197e+002; foe(n+1)=3.145225e+002; krok(n+1)=1.597633e-003; ng(n+1)=2.629739e+001;
n=45; farx(n+1)=3.218698e+002; foe(n+1)=3.142784e+002; krok(n+1)=6.465635e-004; ng(n+1)=4.005912e+001;
n=46; farx(n+1)=3.216901e+002; foe(n+1)=3.140357e+002; krok(n+1)=1.593673e-003; ng(n+1)=2.624706e+001;
n=47; farx(n+1)=3.211518e+002; foe(n+1)=3.137929e+002; krok(n+1)=6.459813e-004; ng(n+1)=4.018122e+001;
n=48; farx(n+1)=3.209972e+002; foe(n+1)=3.135522e+002; krok(n+1)=1.572201e-003; ng(n+1)=2.642439e+001;
n=49; farx(n+1)=3.204710e+002; foe(n+1)=3.133121e+002; krok(n+1)=6.483498e-004; ng(n+1)=4.005797e+001;
n=50; farx(n+1)=3.203390e+002; foe(n+1)=3.130732e+002; krok(n+1)=1.554197e-003; ng(n+1)=2.668783e+001;
n=51; farx(n+1)=3.198209e+002; foe(n+1)=3.128337e+002; krok(n+1)=6.467102e-004; ng(n+1)=4.012408e+001;
n=52; farx(n+1)=3.197040e+002; foe(n+1)=3.125926e+002; krok(n+1)=1.567808e-003; ng(n+1)=2.683891e+001;
n=53; farx(n+1)=3.191854e+002; foe(n+1)=3.123503e+002; krok(n+1)=6.459813e-004; ng(n+1)=4.054756e+001;
n=54; farx(n+1)=3.190860e+002; foe(n+1)=3.121100e+002; krok(n+1)=1.548092e-003; ng(n+1)=2.722108e+001;
n=55; farx(n+1)=3.185699e+002; foe(n+1)=3.118667e+002; krok(n+1)=6.459813e-004; ng(n+1)=4.078833e+001;
n=56; farx(n+1)=3.184814e+002; foe(n+1)=3.116247e+002; krok(n+1)=1.528320e-003; ng(n+1)=2.755454e+001;
n=57; farx(n+1)=3.179693e+002; foe(n+1)=3.113815e+002; krok(n+1)=6.467102e-004; ng(n+1)=4.098837e+001;
n=58; farx(n+1)=3.178916e+002; foe(n+1)=3.111360e+002; krok(n+1)=1.548092e-003; ng(n+1)=2.787046e+001;
n=59; farx(n+1)=3.173717e+002; foe(n+1)=3.108860e+002; krok(n+1)=6.384605e-004; ng(n+1)=4.195695e+001;
n=60; farx(n+1)=3.172964e+002; foe(n+1)=3.106327e+002; krok(n+1)=1.573435e-003; ng(n+1)=2.805687e+001;
n=61; farx(n+1)=3.167697e+002; foe(n+1)=3.103775e+002; krok(n+1)=6.359070e-004; ng(n+1)=4.269438e+001;
n=62; farx(n+1)=3.167008e+002; foe(n+1)=3.101193e+002; krok(n+1)=1.573435e-003; ng(n+1)=2.848204e+001;
n=63; farx(n+1)=3.161677e+002; foe(n+1)=3.098592e+002; krok(n+1)=6.328553e-004; ng(n+1)=4.336561e+001;
n=64; farx(n+1)=3.161033e+002; foe(n+1)=3.095933e+002; krok(n+1)=1.597633e-003; ng(n+1)=2.887471e+001;
n=65; farx(n+1)=3.155563e+002; foe(n+1)=3.093242e+002; krok(n+1)=6.275936e-004; ng(n+1)=4.443950e+001;
n=66; farx(n+1)=3.154959e+002; foe(n+1)=3.090500e+002; krok(n+1)=1.620122e-003; ng(n+1)=2.933570e+001;
n=67; farx(n+1)=3.149317e+002; foe(n+1)=3.087706e+002; krok(n+1)=6.250570e-004; ng(n+1)=4.558023e+001;
n=68; farx(n+1)=3.148754e+002; foe(n+1)=3.084904e+002; krok(n+1)=1.597633e-003; ng(n+1)=3.001338e+001;
n=69; farx(n+1)=3.143021e+002; foe(n+1)=3.082056e+002; krok(n+1)=6.277808e-004; ng(n+1)=4.617294e+001;
n=70; farx(n+1)=3.142503e+002; foe(n+1)=3.079198e+002; krok(n+1)=1.563831e-003; ng(n+1)=3.081438e+001;
n=71; farx(n+1)=3.136688e+002; foe(n+1)=3.076297e+002; krok(n+1)=6.295194e-004; ng(n+1)=4.671234e+001;
n=72; farx(n+1)=3.136209e+002; foe(n+1)=3.073346e+002; krok(n+1)=1.564230e-003; ng(n+1)=3.156130e+001;
n=73; farx(n+1)=3.130213e+002; foe(n+1)=3.070330e+002; krok(n+1)=6.254538e-004; ng(n+1)=4.785254e+001;
n=74; farx(n+1)=3.129730e+002; foe(n+1)=3.067242e+002; krok(n+1)=1.581052e-003; ng(n+1)=3.224034e+001;
n=75; farx(n+1)=3.123497e+002; foe(n+1)=3.064087e+002; krok(n+1)=6.250570e-004; ng(n+1)=4.917007e+001;
n=76; farx(n+1)=3.123026e+002; foe(n+1)=3.060908e+002; krok(n+1)=1.532691e-003; ng(n+1)=3.327038e+001;
n=77; farx(n+1)=3.116692e+002; foe(n+1)=3.057693e+002; krok(n+1)=6.326141e-004; ng(n+1)=4.962908e+001;
n=78; farx(n+1)=3.116305e+002; foe(n+1)=3.054434e+002; krok(n+1)=1.499367e-003; ng(n+1)=3.448168e+001;
n=79; farx(n+1)=3.109767e+002; foe(n+1)=3.051089e+002; krok(n+1)=6.326141e-004; ng(n+1)=5.074128e+001;
n=80; farx(n+1)=3.109375e+002; foe(n+1)=3.047684e+002; krok(n+1)=1.478815e-003; ng(n+1)=3.552711e+001;
n=81; farx(n+1)=3.102658e+002; foe(n+1)=3.044215e+002; krok(n+1)=6.328553e-004; ng(n+1)=5.171263e+001;
n=82; farx(n+1)=3.102285e+002; foe(n+1)=3.040633e+002; krok(n+1)=1.486749e-003; ng(n+1)=3.655883e+001;
n=83; farx(n+1)=3.095241e+002; foe(n+1)=3.036969e+002; krok(n+1)=6.326141e-004; ng(n+1)=5.335040e+001;
n=84; farx(n+1)=3.094931e+002; foe(n+1)=3.033259e+002; krok(n+1)=1.457389e-003; ng(n+1)=3.790610e+001;
n=85; farx(n+1)=3.087642e+002; foe(n+1)=3.029442e+002; krok(n+1)=6.326141e-004; ng(n+1)=5.453769e+001;
n=86; farx(n+1)=3.087346e+002; foe(n+1)=3.025587e+002; krok(n+1)=1.441974e-003; ng(n+1)=3.905075e+001;
n=87; farx(n+1)=3.079841e+002; foe(n+1)=3.021644e+002; krok(n+1)=6.277588e-004; ng(n+1)=5.556170e+001;
n=88; farx(n+1)=3.079469e+002; foe(n+1)=3.017697e+002; krok(n+1)=1.441974e-003; ng(n+1)=3.963553e+001;
n=89; farx(n+1)=3.071900e+002; foe(n+1)=3.013791e+002; krok(n+1)=6.174461e-004; ng(n+1)=5.560542e+001;
n=90; farx(n+1)=3.071224e+002; foe(n+1)=3.010209e+002; krok(n+1)=1.358387e-003; ng(n+1)=3.865648e+001;
n=91; farx(n+1)=3.064701e+002; foe(n+1)=3.007135e+002; krok(n+1)=6.106204e-004; ng(n+1)=4.895457e+001;
n=92; farx(n+1)=3.064090e+002; foe(n+1)=3.004592e+002; krok(n+1)=1.293420e-003; ng(n+1)=3.238180e+001;
n=93; farx(n+1)=3.058946e+002; foe(n+1)=3.002343e+002; krok(n+1)=6.453470e-004; ng(n+1)=3.908971e+001;
n=94; farx(n+1)=3.057728e+002; foe(n+1)=3.000124e+002; krok(n+1)=1.297613e-003; ng(n+1)=2.821270e+001;
n=95; farx(n+1)=3.052117e+002; foe(n+1)=2.997929e+002; krok(n+1)=6.459813e-004; ng(n+1)=3.980262e+001;
n=96; farx(n+1)=3.049798e+002; foe(n+1)=2.995755e+002; krok(n+1)=1.280854e-003; ng(n+1)=2.757176e+001;
n=97; farx(n+1)=3.044036e+002; foe(n+1)=2.993583e+002; krok(n+1)=6.423816e-004; ng(n+1)=4.026636e+001;
n=98; farx(n+1)=3.040859e+002; foe(n+1)=2.991389e+002; krok(n+1)=1.300172e-003; ng(n+1)=2.720028e+001;
n=99; farx(n+1)=3.035089e+002; foe(n+1)=2.989194e+002; krok(n+1)=6.359070e-004; ng(n+1)=4.103843e+001;
n=100; farx(n+1)=3.031209e+002; foe(n+1)=2.986982e+002; krok(n+1)=1.303632e-003; ng(n+1)=2.743383e+001;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)