%uczenie predyktora arx
clear all;
n=0; farx(n+1)=4.934405e+003; foe(n+1)=4.765867e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.585896e+003; foe(n+1)=3.890223e+003; krok(n+1)=2.884695e-004; ng(n+1)=8.706314e+003;
n=2; farx(n+1)=1.902549e+003; foe(n+1)=8.004362e+003; krok(n+1)=2.989017e-004; ng(n+1)=9.261672e+003;
n=3; farx(n+1)=1.174372e+003; foe(n+1)=1.728259e+004; krok(n+1)=1.691914e-004; ng(n+1)=1.113463e+004;
n=4; farx(n+1)=8.484644e+002; foe(n+1)=1.552219e+004; krok(n+1)=1.011508e-003; ng(n+1)=7.792905e+003;
n=5; farx(n+1)=3.628236e+002; foe(n+1)=9.541040e+003; krok(n+1)=1.409080e-003; ng(n+1)=7.255536e+003;
n=6; farx(n+1)=1.347643e+002; foe(n+1)=9.925065e+003; krok(n+1)=2.501815e-003; ng(n+1)=2.578706e+003;
n=7; farx(n+1)=1.137083e+002; foe(n+1)=9.103114e+003; krok(n+1)=1.389271e-004; ng(n+1)=2.208620e+003;
n=8; farx(n+1)=7.915717e+001; foe(n+1)=8.984229e+003; krok(n+1)=1.045638e-003; ng(n+1)=1.206313e+003;
n=9; farx(n+1)=1.910963e+001; foe(n+1)=1.274518e+003; krok(n+1)=1.174683e-003; ng(n+1)=1.445926e+003;
n=10; farx(n+1)=9.675489e+000; foe(n+1)=1.207410e+003; krok(n+1)=3.065383e-003; ng(n+1)=1.128620e+003;
n=11; farx(n+1)=5.524469e+000; foe(n+1)=1.378657e+002; krok(n+1)=6.748373e-003; ng(n+1)=3.562518e+002;
n=12; farx(n+1)=3.219055e+000; foe(n+1)=8.166376e+001; krok(n+1)=4.449709e-003; ng(n+1)=4.070434e+002;
n=13; farx(n+1)=2.594393e+000; foe(n+1)=7.716824e+001; krok(n+1)=5.445782e-003; ng(n+1)=2.121679e+002;
n=14; farx(n+1)=2.094123e+000; foe(n+1)=1.148824e+002; krok(n+1)=6.900644e-002; ng(n+1)=1.432116e+002;
n=15; farx(n+1)=1.864495e+000; foe(n+1)=1.093496e+002; krok(n+1)=8.877087e-002; ng(n+1)=9.578631e+001;
n=16; farx(n+1)=1.719244e+000; foe(n+1)=1.055399e+002; krok(n+1)=3.836134e-002; ng(n+1)=2.660356e+001;
n=17; farx(n+1)=1.638811e+000; foe(n+1)=1.057202e+002; krok(n+1)=6.524238e-002; ng(n+1)=3.317907e+001;
n=18; farx(n+1)=1.492104e+000; foe(n+1)=9.723342e+001; krok(n+1)=1.105271e-001; ng(n+1)=8.964560e+001;
n=19; farx(n+1)=1.203207e+000; foe(n+1)=1.042087e+002; krok(n+1)=1.006999e-001; ng(n+1)=1.555406e+002;
n=20; farx(n+1)=1.044949e+000; foe(n+1)=1.032908e+002; krok(n+1)=1.115271e-001; ng(n+1)=4.612191e+001;
n=21; farx(n+1)=9.703054e-001; foe(n+1)=1.101627e+002; krok(n+1)=6.640781e-002; ng(n+1)=9.500239e+001;
n=22; farx(n+1)=9.296843e-001; foe(n+1)=9.309323e+001; krok(n+1)=1.253365e-001; ng(n+1)=8.186766e+001;
n=23; farx(n+1)=8.496474e-001; foe(n+1)=8.385253e+001; krok(n+1)=1.638325e-001; ng(n+1)=8.747227e+001;
n=24; farx(n+1)=7.509966e-001; foe(n+1)=7.508167e+001; krok(n+1)=3.477470e-001; ng(n+1)=1.151050e+002;
n=25; farx(n+1)=7.247591e-001; foe(n+1)=7.318302e+001; krok(n+1)=1.865388e-001; ng(n+1)=6.098322e+001;
%odnowa zmiennej metryki
n=26; farx(n+1)=7.179998e-001; foe(n+1)=7.062539e+001; krok(n+1)=5.434585e-005; ng(n+1)=4.260901e+001;
n=27; farx(n+1)=7.170089e-001; foe(n+1)=7.113670e+001; krok(n+1)=1.338340e-004; ng(n+1)=1.223338e+001;
n=28; farx(n+1)=7.034522e-001; foe(n+1)=6.452561e+001; krok(n+1)=7.088910e-003; ng(n+1)=6.206871e+000;
n=29; farx(n+1)=6.739218e-001; foe(n+1)=8.593348e+001; krok(n+1)=4.829711e-003; ng(n+1)=1.257140e+001;
n=30; farx(n+1)=6.628460e-001; foe(n+1)=8.709461e+001; krok(n+1)=8.360316e-003; ng(n+1)=1.940085e+001;
n=31; farx(n+1)=6.507023e-001; foe(n+1)=7.465320e+001; krok(n+1)=1.121678e-002; ng(n+1)=2.826898e+001;
n=32; farx(n+1)=6.400601e-001; foe(n+1)=7.516493e+001; krok(n+1)=2.835564e-002; ng(n+1)=4.515317e+001;
n=33; farx(n+1)=6.257565e-001; foe(n+1)=7.141697e+001; krok(n+1)=9.909413e-002; ng(n+1)=6.569551e+001;
n=34; farx(n+1)=5.991334e-001; foe(n+1)=6.285419e+001; krok(n+1)=3.311156e-001; ng(n+1)=5.250596e+001;
n=35; farx(n+1)=5.924703e-001; foe(n+1)=6.283182e+001; krok(n+1)=5.927401e-002; ng(n+1)=3.875286e+001;
n=36; farx(n+1)=5.809265e-001; foe(n+1)=5.445851e+001; krok(n+1)=2.577175e-001; ng(n+1)=1.461830e+001;
n=37; farx(n+1)=5.644249e-001; foe(n+1)=5.405076e+001; krok(n+1)=6.754317e-002; ng(n+1)=2.005514e+001;
n=38; farx(n+1)=5.484231e-001; foe(n+1)=5.802690e+001; krok(n+1)=8.268560e-002; ng(n+1)=3.945820e+001;
n=39; farx(n+1)=5.268622e-001; foe(n+1)=5.604958e+001; krok(n+1)=3.923690e-001; ng(n+1)=4.317686e+001;
n=40; farx(n+1)=5.124104e-001; foe(n+1)=6.625903e+001; krok(n+1)=1.261335e-001; ng(n+1)=4.051452e+001;
n=41; farx(n+1)=4.998464e-001; foe(n+1)=6.127825e+001; krok(n+1)=3.057702e-001; ng(n+1)=1.930981e+001;
n=42; farx(n+1)=4.810455e-001; foe(n+1)=5.905426e+001; krok(n+1)=2.745494e-001; ng(n+1)=1.990212e+001;
n=43; farx(n+1)=4.763466e-001; foe(n+1)=6.332110e+001; krok(n+1)=7.266281e-002; ng(n+1)=3.207815e+001;
n=44; farx(n+1)=4.666628e-001; foe(n+1)=5.505704e+001; krok(n+1)=1.516292e-001; ng(n+1)=2.328867e+001;
n=45; farx(n+1)=4.519210e-001; foe(n+1)=4.885522e+001; krok(n+1)=5.131430e-001; ng(n+1)=1.189983e+001;
n=46; farx(n+1)=4.452870e-001; foe(n+1)=4.412904e+001; krok(n+1)=3.064368e-001; ng(n+1)=1.559590e+001;
n=47; farx(n+1)=4.405753e-001; foe(n+1)=5.916083e+001; krok(n+1)=2.760258e-001; ng(n+1)=4.891946e+000;
n=48; farx(n+1)=4.356457e-001; foe(n+1)=8.846496e+001; krok(n+1)=3.444314e-001; ng(n+1)=1.207215e+001;
n=49; farx(n+1)=4.287792e-001; foe(n+1)=1.046176e+002; krok(n+1)=7.021280e-001; ng(n+1)=1.187961e+001;
n=50; farx(n+1)=4.242171e-001; foe(n+1)=1.427217e+002; krok(n+1)=4.063898e-001; ng(n+1)=6.077636e+000;
%odnowa zmiennej metryki
n=51; farx(n+1)=4.229193e-001; foe(n+1)=9.911895e+001; krok(n+1)=6.946355e-005; ng(n+1)=2.026561e+001;
n=52; farx(n+1)=4.224668e-001; foe(n+1)=9.642609e+001; krok(n+1)=7.910691e-005; ng(n+1)=9.692260e+000;
n=53; farx(n+1)=4.220013e-001; foe(n+1)=9.581412e+001; krok(n+1)=2.707062e-003; ng(n+1)=2.120030e+000;
n=54; farx(n+1)=4.214446e-001; foe(n+1)=7.789369e+001; krok(n+1)=3.111744e-003; ng(n+1)=2.014385e+000;
n=55; farx(n+1)=4.209566e-001; foe(n+1)=4.726484e+001; krok(n+1)=4.687691e-003; ng(n+1)=1.416104e+000;
n=56; farx(n+1)=4.194908e-001; foe(n+1)=4.798834e+001; krok(n+1)=9.091272e-002; ng(n+1)=6.662994e-001;
n=57; farx(n+1)=4.164903e-001; foe(n+1)=4.188739e+001; krok(n+1)=3.571281e-002; ng(n+1)=2.038116e+000;
n=58; farx(n+1)=4.140630e-001; foe(n+1)=3.430393e+001; krok(n+1)=4.215364e-002; ng(n+1)=1.062876e+001;
n=59; farx(n+1)=4.114832e-001; foe(n+1)=2.909072e+001; krok(n+1)=2.670983e-001; ng(n+1)=2.182625e+001;
n=60; farx(n+1)=4.052347e-001; foe(n+1)=3.281948e+001; krok(n+1)=3.444314e-001; ng(n+1)=2.140075e+001;
n=61; farx(n+1)=3.970709e-001; foe(n+1)=2.164041e+001; krok(n+1)=3.103642e-001; ng(n+1)=1.170377e+001;
n=62; farx(n+1)=3.902469e-001; foe(n+1)=1.557296e+001; krok(n+1)=1.439513e-001; ng(n+1)=1.576135e+001;
n=63; farx(n+1)=3.874462e-001; foe(n+1)=1.455449e+001; krok(n+1)=1.024273e-001; ng(n+1)=1.974335e+001;
n=64; farx(n+1)=3.797634e-001; foe(n+1)=1.175643e+001; krok(n+1)=2.962682e-001; ng(n+1)=5.687589e+000;
n=65; farx(n+1)=3.752907e-001; foe(n+1)=1.363406e+001; krok(n+1)=1.738735e-001; ng(n+1)=6.462923e+000;
n=66; farx(n+1)=3.658796e-001; foe(n+1)=1.558729e+001; krok(n+1)=2.656028e-001; ng(n+1)=1.264908e+001;
n=67; farx(n+1)=3.611428e-001; foe(n+1)=1.624663e+001; krok(n+1)=1.943277e-001; ng(n+1)=1.444048e+001;
n=68; farx(n+1)=3.557414e-001; foe(n+1)=1.439266e+001; krok(n+1)=1.775417e-001; ng(n+1)=1.041833e+001;
n=69; farx(n+1)=3.519865e-001; foe(n+1)=1.406988e+001; krok(n+1)=8.509079e-002; ng(n+1)=1.285664e+001;
n=70; farx(n+1)=3.484690e-001; foe(n+1)=1.284086e+001; krok(n+1)=1.660944e-001; ng(n+1)=4.266491e+000;
n=71; farx(n+1)=3.395826e-001; foe(n+1)=9.978221e+000; krok(n+1)=2.656312e-001; ng(n+1)=2.597651e+001;
n=72; farx(n+1)=3.314379e-001; foe(n+1)=1.007627e+001; krok(n+1)=4.386596e-001; ng(n+1)=2.098240e+001;
n=73; farx(n+1)=3.263020e-001; foe(n+1)=8.805722e+000; krok(n+1)=4.681360e-001; ng(n+1)=7.349754e+000;
n=74; farx(n+1)=3.224141e-001; foe(n+1)=8.086913e+000; krok(n+1)=3.072232e-001; ng(n+1)=1.128437e+001;
n=75; farx(n+1)=3.169901e-001; foe(n+1)=6.760866e+000; krok(n+1)=1.663330e-001; ng(n+1)=1.068200e+001;
%odnowa zmiennej metryki
n=76; farx(n+1)=3.166130e-001; foe(n+1)=6.786828e+000; krok(n+1)=4.041939e-005; ng(n+1)=1.133763e+001;
n=77; farx(n+1)=3.160990e-001; foe(n+1)=6.793173e+000; krok(n+1)=1.084015e-003; ng(n+1)=3.116550e+000;
n=78; farx(n+1)=3.158997e-001; foe(n+1)=6.872443e+000; krok(n+1)=3.550621e-004; ng(n+1)=3.490483e+000;
n=79; farx(n+1)=3.134465e-001; foe(n+1)=7.694602e+000; krok(n+1)=4.672280e-003; ng(n+1)=3.111958e+000;
n=80; farx(n+1)=3.129476e-001; foe(n+1)=7.994814e+000; krok(n+1)=1.848181e-002; ng(n+1)=8.230463e-001;
n=81; farx(n+1)=3.127421e-001; foe(n+1)=7.875328e+000; krok(n+1)=4.541426e-003; ng(n+1)=9.737899e-001;
n=82; farx(n+1)=3.104116e-001; foe(n+1)=7.875688e+000; krok(n+1)=1.036878e-001; ng(n+1)=7.847442e-001;
n=83; farx(n+1)=3.086354e-001; foe(n+1)=8.375227e+000; krok(n+1)=1.399123e-001; ng(n+1)=5.552045e+000;
n=84; farx(n+1)=3.064778e-001; foe(n+1)=7.893637e+000; krok(n+1)=3.047856e-001; ng(n+1)=1.537725e+001;
n=85; farx(n+1)=3.050587e-001; foe(n+1)=7.376666e+000; krok(n+1)=9.048247e-002; ng(n+1)=9.891487e+000;
n=86; farx(n+1)=2.975251e-001; foe(n+1)=4.920103e+000; krok(n+1)=6.101385e-001; ng(n+1)=1.129582e+001;
n=87; farx(n+1)=2.926390e-001; foe(n+1)=4.244576e+000; krok(n+1)=1.305790e-001; ng(n+1)=2.144599e+001;
n=88; farx(n+1)=2.881063e-001; foe(n+1)=3.959683e+000; krok(n+1)=1.305790e-001; ng(n+1)=2.692047e+001;
n=89; farx(n+1)=2.850608e-001; foe(n+1)=3.736333e+000; krok(n+1)=1.139041e-001; ng(n+1)=1.486654e+001;
n=90; farx(n+1)=2.818376e-001; foe(n+1)=3.613125e+000; krok(n+1)=1.084302e-001; ng(n+1)=8.905595e+000;
n=91; farx(n+1)=2.761138e-001; foe(n+1)=3.125709e+000; krok(n+1)=2.561085e-001; ng(n+1)=9.847079e+000;
n=92; farx(n+1)=2.682923e-001; foe(n+1)=3.271031e+000; krok(n+1)=5.945315e-001; ng(n+1)=1.437759e+001;
n=93; farx(n+1)=2.613163e-001; foe(n+1)=2.673907e+000; krok(n+1)=3.268918e-001; ng(n+1)=1.654855e+001;
n=94; farx(n+1)=2.571336e-001; foe(n+1)=2.356767e+000; krok(n+1)=1.259999e-001; ng(n+1)=1.191133e+001;
n=95; farx(n+1)=2.551325e-001; foe(n+1)=2.213424e+000; krok(n+1)=1.591991e-001; ng(n+1)=7.103610e+000;
n=96; farx(n+1)=2.515337e-001; foe(n+1)=1.878108e+000; krok(n+1)=1.196783e-001; ng(n+1)=1.381005e+001;
n=97; farx(n+1)=2.493522e-001; foe(n+1)=1.608573e+000; krok(n+1)=3.655097e-001; ng(n+1)=1.095391e+001;
n=98; farx(n+1)=2.466487e-001; foe(n+1)=1.646919e+000; krok(n+1)=2.489869e-001; ng(n+1)=7.295676e+000;
n=99; farx(n+1)=2.433248e-001; foe(n+1)=1.397941e+000; krok(n+1)=5.131430e-001; ng(n+1)=7.974791e+000;
n=100; farx(n+1)=2.395257e-001; foe(n+1)=1.340396e+000; krok(n+1)=4.164005e-001; ng(n+1)=3.157457e+000;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora ARX');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)