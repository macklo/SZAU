%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.685641e+003; foe(n+1)=4.927584e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.835331e+003; foe(n+1)=3.823212e+003; krok(n+1)=3.325065e-004; ng(n+1)=6.119334e+003;
n=2; farx(n+1)=8.909541e+002; foe(n+1)=8.042168e+002; krok(n+1)=2.586841e-003; ng(n+1)=2.833518e+003;
n=3; farx(n+1)=2.745166e+002; foe(n+1)=5.388356e+002; krok(n+1)=6.480489e-003; ng(n+1)=1.058407e+003;
n=4; farx(n+1)=2.050712e+002; foe(n+1)=4.396762e+002; krok(n+1)=2.952483e-005; ng(n+1)=1.225524e+004;
n=5; farx(n+1)=1.727821e+002; foe(n+1)=3.381995e+002; krok(n+1)=2.215284e-004; ng(n+1)=8.393296e+003;
n=6; farx(n+1)=1.665356e+002; foe(n+1)=3.310467e+002; krok(n+1)=7.770983e-004; ng(n+1)=1.240331e+003;
n=7; farx(n+1)=1.824025e+002; foe(n+1)=3.038796e+002; krok(n+1)=1.675019e-003; ng(n+1)=1.763188e+003;
n=8; farx(n+1)=1.543334e+002; foe(n+1)=2.902117e+002; krok(n+1)=7.989279e-004; ng(n+1)=1.482448e+003;
n=9; farx(n+1)=1.258443e+002; foe(n+1)=2.668238e+002; krok(n+1)=8.259539e-004; ng(n+1)=2.156281e+003;
n=10; farx(n+1)=1.176897e+002; foe(n+1)=2.610262e+002; krok(n+1)=6.710424e-004; ng(n+1)=3.613476e+003;
n=11; farx(n+1)=7.745722e+001; foe(n+1)=2.288585e+002; krok(n+1)=9.053002e-004; ng(n+1)=3.701413e+003;
n=12; farx(n+1)=6.684954e+001; foe(n+1)=2.211055e+002; krok(n+1)=9.050205e-005; ng(n+1)=4.509076e+003;
n=13; farx(n+1)=4.691905e+001; foe(n+1)=2.107071e+002; krok(n+1)=4.689376e-004; ng(n+1)=5.599931e+003;
n=14; farx(n+1)=3.925914e+001; foe(n+1)=2.063966e+002; krok(n+1)=1.874208e-004; ng(n+1)=4.669476e+003;
n=15; farx(n+1)=2.779712e+001; foe(n+1)=1.937819e+002; krok(n+1)=5.593002e-004; ng(n+1)=5.844414e+003;
n=16; farx(n+1)=2.566241e+001; foe(n+1)=1.902261e+002; krok(n+1)=5.534736e-004; ng(n+1)=3.539231e+003;
n=17; farx(n+1)=1.469171e+001; foe(n+1)=1.609479e+002; krok(n+1)=1.970836e-003; ng(n+1)=4.111446e+003;
n=18; farx(n+1)=1.178586e+001; foe(n+1)=1.480031e+002; krok(n+1)=5.081655e-004; ng(n+1)=4.094923e+003;
n=19; farx(n+1)=1.708182e+001; foe(n+1)=1.170914e+002; krok(n+1)=5.333277e-004; ng(n+1)=8.418412e+003;
n=20; farx(n+1)=2.066076e+001; foe(n+1)=1.123403e+002; krok(n+1)=3.933588e-004; ng(n+1)=3.254652e+003;
n=21; farx(n+1)=2.392570e+001; foe(n+1)=1.095356e+002; krok(n+1)=2.732662e-004; ng(n+1)=4.702482e+003;
n=22; farx(n+1)=2.506922e+001; foe(n+1)=1.063179e+002; krok(n+1)=3.334617e-003; ng(n+1)=2.472704e+003;
n=23; farx(n+1)=3.198013e+001; foe(n+1)=9.462832e+001; krok(n+1)=5.547107e-004; ng(n+1)=4.027858e+003;
n=24; farx(n+1)=3.228723e+001; foe(n+1)=8.875850e+001; krok(n+1)=1.310706e-003; ng(n+1)=4.672142e+003;
n=25; farx(n+1)=3.495135e+001; foe(n+1)=8.562466e+001; krok(n+1)=4.261436e-004; ng(n+1)=7.491995e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=3.154310e+001; foe(n+1)=7.045082e+001; krok(n+1)=1.095389e-005; ng(n+1)=7.319922e+003;
n=27; farx(n+1)=3.145640e+001; foe(n+1)=6.775449e+001; krok(n+1)=2.870650e-006; ng(n+1)=6.117917e+003;
n=28; farx(n+1)=2.421611e+001; foe(n+1)=5.704291e+001; krok(n+1)=2.025153e-004; ng(n+1)=1.500256e+003;
n=29; farx(n+1)=2.407671e+001; foe(n+1)=5.247364e+001; krok(n+1)=7.318567e-005; ng(n+1)=1.462015e+003;
n=30; farx(n+1)=2.198943e+001; foe(n+1)=4.696015e+001; krok(n+1)=3.870231e-004; ng(n+1)=1.427368e+003;
n=31; farx(n+1)=1.892804e+001; foe(n+1)=4.188020e+001; krok(n+1)=5.295816e-004; ng(n+1)=1.940996e+003;
n=32; farx(n+1)=1.576301e+001; foe(n+1)=3.470892e+001; krok(n+1)=2.141344e-003; ng(n+1)=4.520985e+003;
n=33; farx(n+1)=1.460687e+001; foe(n+1)=3.228393e+001; krok(n+1)=1.037622e-003; ng(n+1)=1.630875e+003;
n=34; farx(n+1)=1.164457e+001; foe(n+1)=2.742919e+001; krok(n+1)=5.557598e-003; ng(n+1)=1.723559e+003;
n=35; farx(n+1)=1.060401e+001; foe(n+1)=2.560429e+001; krok(n+1)=1.751967e-003; ng(n+1)=1.359526e+003;
n=36; farx(n+1)=9.567558e+000; foe(n+1)=2.335240e+001; krok(n+1)=1.040544e-003; ng(n+1)=2.436427e+003;
n=37; farx(n+1)=7.583366e+000; foe(n+1)=1.920964e+001; krok(n+1)=7.517005e-003; ng(n+1)=8.745615e+002;
n=38; farx(n+1)=6.361883e+000; foe(n+1)=1.689039e+001; krok(n+1)=1.532691e-003; ng(n+1)=1.073897e+003;
n=39; farx(n+1)=5.164325e+000; foe(n+1)=1.366529e+001; krok(n+1)=2.156451e-003; ng(n+1)=1.820167e+003;
n=40; farx(n+1)=4.979219e+000; foe(n+1)=1.330428e+001; krok(n+1)=8.751134e-004; ng(n+1)=7.261275e+002;
n=41; farx(n+1)=4.411379e+000; foe(n+1)=1.273886e+001; krok(n+1)=1.310706e-003; ng(n+1)=9.927432e+002;
n=42; farx(n+1)=3.825548e+000; foe(n+1)=1.191148e+001; krok(n+1)=6.544675e-003; ng(n+1)=4.712012e+002;
n=43; farx(n+1)=3.270915e+000; foe(n+1)=1.064065e+001; krok(n+1)=7.739087e-003; ng(n+1)=4.500775e+002;
n=44; farx(n+1)=3.188227e+000; foe(n+1)=1.013742e+001; krok(n+1)=1.041670e-002; ng(n+1)=7.002492e+002;
n=45; farx(n+1)=3.101408e+000; foe(n+1)=9.781847e+000; krok(n+1)=3.253129e-003; ng(n+1)=5.778237e+002;
n=46; farx(n+1)=3.046467e+000; foe(n+1)=9.392279e+000; krok(n+1)=3.683227e-003; ng(n+1)=9.000375e+002;
n=47; farx(n+1)=2.939058e+000; foe(n+1)=8.944704e+000; krok(n+1)=6.475792e-003; ng(n+1)=6.864067e+002;
n=48; farx(n+1)=2.885766e+000; foe(n+1)=8.655458e+000; krok(n+1)=8.252113e-003; ng(n+1)=6.300743e+002;
n=49; farx(n+1)=2.899345e+000; foe(n+1)=8.105219e+000; krok(n+1)=8.806488e-003; ng(n+1)=7.226793e+002;
n=50; farx(n+1)=2.896778e+000; foe(n+1)=7.834260e+000; krok(n+1)=1.039581e-002; ng(n+1)=3.937182e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=2.898667e+000; foe(n+1)=7.761419e+000; krok(n+1)=8.746968e-006; ng(n+1)=4.635147e+002;
n=52; farx(n+1)=2.901707e+000; foe(n+1)=7.650043e+000; krok(n+1)=3.418410e-005; ng(n+1)=3.030412e+002;
n=53; farx(n+1)=2.865736e+000; foe(n+1)=7.384996e+000; krok(n+1)=8.147699e-005; ng(n+1)=2.905879e+002;
n=54; farx(n+1)=2.754844e+000; foe(n+1)=7.059784e+000; krok(n+1)=1.120093e-004; ng(n+1)=3.019430e+002;
n=55; farx(n+1)=2.626174e+000; foe(n+1)=6.654049e+000; krok(n+1)=2.218714e-004; ng(n+1)=2.524780e+002;
n=56; farx(n+1)=2.574807e+000; foe(n+1)=6.389437e+000; krok(n+1)=1.025308e-003; ng(n+1)=1.756016e+002;
n=57; farx(n+1)=2.476552e+000; foe(n+1)=5.892346e+000; krok(n+1)=2.241713e-003; ng(n+1)=2.634714e+002;
n=58; farx(n+1)=2.094777e+000; foe(n+1)=5.170663e+000; krok(n+1)=1.417782e-002; ng(n+1)=6.099994e+002;
n=59; farx(n+1)=1.934538e+000; foe(n+1)=4.691293e+000; krok(n+1)=7.354218e-003; ng(n+1)=3.943823e+002;
n=60; farx(n+1)=1.915008e+000; foe(n+1)=4.543296e+000; krok(n+1)=1.772228e-003; ng(n+1)=5.866592e+002;
n=61; farx(n+1)=1.825273e+000; foe(n+1)=4.390977e+000; krok(n+1)=4.698731e-003; ng(n+1)=4.044327e+002;
n=62; farx(n+1)=1.593271e+000; foe(n+1)=3.901188e+000; krok(n+1)=4.893842e-003; ng(n+1)=6.833260e+002;
n=63; farx(n+1)=1.407099e+000; foe(n+1)=3.584160e+000; krok(n+1)=4.334641e-003; ng(n+1)=3.961786e+002;
n=64; farx(n+1)=1.315953e+000; foe(n+1)=3.375059e+000; krok(n+1)=5.557598e-003; ng(n+1)=5.593710e+002;
n=65; farx(n+1)=1.282980e+000; foe(n+1)=3.316518e+000; krok(n+1)=1.425029e-003; ng(n+1)=3.331783e+002;
n=66; farx(n+1)=1.221579e+000; foe(n+1)=3.186814e+000; krok(n+1)=6.746361e-003; ng(n+1)=5.068306e+002;
n=67; farx(n+1)=1.190470e+000; foe(n+1)=3.017576e+000; krok(n+1)=1.204622e-002; ng(n+1)=1.943853e+002;
n=68; farx(n+1)=1.114057e+000; foe(n+1)=2.811367e+000; krok(n+1)=1.555551e-002; ng(n+1)=3.425595e+002;
n=69; farx(n+1)=1.053784e+000; foe(n+1)=2.616327e+000; krok(n+1)=2.147336e-002; ng(n+1)=3.682978e+002;
n=70; farx(n+1)=1.025698e+000; foe(n+1)=2.428597e+000; krok(n+1)=3.715831e-003; ng(n+1)=3.697577e+002;
n=71; farx(n+1)=1.011719e+000; foe(n+1)=2.337836e+000; krok(n+1)=3.026819e-003; ng(n+1)=4.083971e+002;
n=72; farx(n+1)=9.856517e-001; foe(n+1)=2.183221e+000; krok(n+1)=6.543079e-003; ng(n+1)=3.056978e+002;
n=73; farx(n+1)=9.307663e-001; foe(n+1)=1.920576e+000; krok(n+1)=2.083339e-002; ng(n+1)=4.651693e+002;
n=74; farx(n+1)=9.046208e-001; foe(n+1)=1.868985e+000; krok(n+1)=9.871862e-003; ng(n+1)=2.187859e+002;
n=75; farx(n+1)=8.597626e-001; foe(n+1)=1.785945e+000; krok(n+1)=8.580693e-003; ng(n+1)=3.265280e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=8.557318e-001; foe(n+1)=1.765947e+000; krok(n+1)=6.639014e-006; ng(n+1)=2.899545e+002;
n=77; farx(n+1)=8.488309e-001; foe(n+1)=1.711715e+000; krok(n+1)=3.014162e-005; ng(n+1)=2.552236e+002;
n=78; farx(n+1)=8.454953e-001; foe(n+1)=1.687639e+000; krok(n+1)=9.600090e-006; ng(n+1)=2.662662e+002;
n=79; farx(n+1)=8.470017e-001; foe(n+1)=1.681048e+000; krok(n+1)=3.616555e-005; ng(n+1)=7.983926e+001;
n=80; farx(n+1)=8.508682e-001; foe(n+1)=1.644154e+000; krok(n+1)=6.935224e-004; ng(n+1)=4.161010e+001;
n=81; farx(n+1)=8.545871e-001; foe(n+1)=1.624241e+000; krok(n+1)=4.437427e-004; ng(n+1)=4.452369e+001;
n=82; farx(n+1)=8.618646e-001; foe(n+1)=1.576425e+000; krok(n+1)=1.427772e-003; ng(n+1)=5.478047e+001;
n=83; farx(n+1)=8.674508e-001; foe(n+1)=1.536740e+000; krok(n+1)=2.620821e-003; ng(n+1)=1.576705e+002;
n=84; farx(n+1)=8.651609e-001; foe(n+1)=1.509763e+000; krok(n+1)=2.885366e-003; ng(n+1)=3.550270e+002;
n=85; farx(n+1)=8.304888e-001; foe(n+1)=1.422636e+000; krok(n+1)=1.813820e-002; ng(n+1)=5.261929e+002;
n=86; farx(n+1)=8.087664e-001; foe(n+1)=1.338920e+000; krok(n+1)=8.510569e-003; ng(n+1)=4.223207e+002;
n=87; farx(n+1)=7.687812e-001; foe(n+1)=1.283149e+000; krok(n+1)=9.299198e-003; ng(n+1)=8.911737e+001;
n=88; farx(n+1)=7.547752e-001; foe(n+1)=1.263140e+000; krok(n+1)=3.065383e-003; ng(n+1)=3.047306e+002;
n=89; farx(n+1)=7.303723e-001; foe(n+1)=1.237693e+000; krok(n+1)=8.443488e-003; ng(n+1)=2.135261e+002;
n=90; farx(n+1)=6.989901e-001; foe(n+1)=1.170770e+000; krok(n+1)=1.170476e-002; ng(n+1)=3.623526e+002;
n=91; farx(n+1)=6.697257e-001; foe(n+1)=1.131965e+000; krok(n+1)=1.748904e-002; ng(n+1)=2.778477e+002;
n=92; farx(n+1)=6.522488e-001; foe(n+1)=1.084807e+000; krok(n+1)=1.109636e-002; ng(n+1)=2.599634e+002;
n=93; farx(n+1)=6.380812e-001; foe(n+1)=1.029158e+000; krok(n+1)=2.447549e-002; ng(n+1)=3.331794e+002;
n=94; farx(n+1)=6.301442e-001; foe(n+1)=9.970788e-001; krok(n+1)=6.026360e-003; ng(n+1)=2.783297e+002;
n=95; farx(n+1)=6.248376e-001; foe(n+1)=9.673315e-001; krok(n+1)=1.214548e-002; ng(n+1)=1.848015e+002;
n=96; farx(n+1)=6.058153e-001; foe(n+1)=9.345533e-001; krok(n+1)=1.450364e-002; ng(n+1)=3.238465e+002;
n=97; farx(n+1)=5.981577e-001; foe(n+1)=9.216151e-001; krok(n+1)=4.683883e-003; ng(n+1)=2.218853e+002;
n=98; farx(n+1)=6.073519e-001; foe(n+1)=8.986842e-001; krok(n+1)=3.413297e-002; ng(n+1)=1.773662e+002;
n=99; farx(n+1)=6.258791e-001; foe(n+1)=8.687567e-001; krok(n+1)=1.355378e-002; ng(n+1)=1.551359e+002;
n=100; farx(n+1)=5.964920e-001; foe(n+1)=8.320985e-001; krok(n+1)=3.847317e-002; ng(n+1)=3.600625e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
