%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.454231e+003; foe(n+1)=4.513363e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.621735e+003; foe(n+1)=2.715807e+003; krok(n+1)=4.957630e-004; ng(n+1)=4.615843e+003;
n=2; farx(n+1)=1.047947e+003; foe(n+1)=8.942483e+002; krok(n+1)=1.797919e-003; ng(n+1)=3.100537e+003;
n=3; farx(n+1)=1.284970e+003; foe(n+1)=5.307514e+002; krok(n+1)=2.563269e-004; ng(n+1)=6.381282e+003;
n=4; farx(n+1)=1.262750e+003; foe(n+1)=4.739297e+002; krok(n+1)=1.111520e-002; ng(n+1)=1.185795e+003;
n=5; farx(n+1)=1.235750e+003; foe(n+1)=4.708455e+002; krok(n+1)=6.738910e-005; ng(n+1)=9.884441e+002;
n=6; farx(n+1)=7.594532e+002; foe(n+1)=4.206757e+002; krok(n+1)=5.945944e-003; ng(n+1)=9.800965e+002;
n=7; farx(n+1)=3.872384e+002; foe(n+1)=3.471760e+002; krok(n+1)=1.043660e-003; ng(n+1)=1.437052e+003;
n=8; farx(n+1)=3.342041e+002; foe(n+1)=3.320969e+002; krok(n+1)=3.691804e-004; ng(n+1)=1.157297e+003;
n=9; farx(n+1)=3.287235e+002; foe(n+1)=3.250755e+002; krok(n+1)=3.187345e-003; ng(n+1)=1.028877e+003;
n=10; farx(n+1)=3.007993e+002; foe(n+1)=3.185035e+002; krok(n+1)=1.121665e-002; ng(n+1)=5.019608e+002;
n=11; farx(n+1)=2.812743e+002; foe(n+1)=3.135368e+002; krok(n+1)=4.838912e-003; ng(n+1)=7.048065e+002;
n=12; farx(n+1)=2.542695e+002; foe(n+1)=3.099363e+002; krok(n+1)=2.684169e-003; ng(n+1)=5.767457e+002;
n=13; farx(n+1)=2.111137e+002; foe(n+1)=2.967482e+002; krok(n+1)=2.531421e-003; ng(n+1)=1.174266e+003;
n=14; farx(n+1)=2.025570e+002; foe(n+1)=2.958206e+002; krok(n+1)=3.164277e-004; ng(n+1)=6.585008e+002;
n=15; farx(n+1)=1.554188e+002; foe(n+1)=2.883989e+002; krok(n+1)=9.328627e-003; ng(n+1)=1.007020e+003;
n=16; farx(n+1)=1.371245e+002; foe(n+1)=2.842171e+002; krok(n+1)=1.203868e-003; ng(n+1)=7.941268e+002;
n=17; farx(n+1)=1.074085e+002; foe(n+1)=2.752158e+002; krok(n+1)=1.548092e-003; ng(n+1)=1.077101e+003;
n=18; farx(n+1)=8.956238e+001; foe(n+1)=2.678658e+002; krok(n+1)=2.240817e-003; ng(n+1)=1.827096e+003;
n=19; farx(n+1)=7.967568e+001; foe(n+1)=2.654892e+002; krok(n+1)=1.330026e-003; ng(n+1)=2.209016e+003;
n=20; farx(n+1)=5.941568e+001; foe(n+1)=2.538989e+002; krok(n+1)=1.041670e-002; ng(n+1)=2.853465e+003;
n=21; farx(n+1)=5.020882e+001; foe(n+1)=2.485173e+002; krok(n+1)=1.944482e-003; ng(n+1)=2.823884e+003;
n=22; farx(n+1)=4.485780e+001; foe(n+1)=2.419991e+002; krok(n+1)=5.945944e-003; ng(n+1)=2.594087e+003;
n=23; farx(n+1)=4.413043e+001; foe(n+1)=2.391985e+002; krok(n+1)=1.573435e-003; ng(n+1)=2.964862e+003;
n=24; farx(n+1)=4.140014e+001; foe(n+1)=2.307429e+002; krok(n+1)=3.966104e-003; ng(n+1)=2.763270e+003;
n=25; farx(n+1)=3.426463e+001; foe(n+1)=2.126980e+002; krok(n+1)=2.107351e-002; ng(n+1)=2.853464e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=3.070304e+001; foe(n+1)=1.700157e+002; krok(n+1)=8.632369e-006; ng(n+1)=7.638605e+003;
n=27; farx(n+1)=3.059454e+001; foe(n+1)=9.647254e+001; krok(n+1)=7.095977e-005; ng(n+1)=3.023087e+003;
n=28; farx(n+1)=2.698514e+001; foe(n+1)=8.851147e+001; krok(n+1)=5.714551e-005; ng(n+1)=6.663074e+003;
n=29; farx(n+1)=2.704263e+001; foe(n+1)=8.342152e+001; krok(n+1)=1.042068e-004; ng(n+1)=9.286755e+003;
n=30; farx(n+1)=2.361780e+001; foe(n+1)=7.334559e+001; krok(n+1)=4.988482e-004; ng(n+1)=6.813427e+003;
n=31; farx(n+1)=2.138111e+001; foe(n+1)=6.463308e+001; krok(n+1)=2.171131e-004; ng(n+1)=6.623872e+003;
n=32; farx(n+1)=2.090847e+001; foe(n+1)=5.755850e+001; krok(n+1)=7.383609e-004; ng(n+1)=2.967382e+003;
n=33; farx(n+1)=1.824090e+001; foe(n+1)=4.395018e+001; krok(n+1)=2.021294e-003; ng(n+1)=3.346278e+003;
n=34; farx(n+1)=1.344657e+001; foe(n+1)=3.757932e+001; krok(n+1)=5.894618e-003; ng(n+1)=1.457961e+003;
n=35; farx(n+1)=1.106535e+001; foe(n+1)=3.273192e+001; krok(n+1)=7.589753e-004; ng(n+1)=2.912392e+003;
n=36; farx(n+1)=9.338842e+000; foe(n+1)=2.896027e+001; krok(n+1)=9.732819e-003; ng(n+1)=1.455735e+003;
n=37; farx(n+1)=9.109776e+000; foe(n+1)=2.802137e+001; krok(n+1)=8.106441e-004; ng(n+1)=1.520220e+003;
n=38; farx(n+1)=8.873988e+000; foe(n+1)=2.559172e+001; krok(n+1)=2.275225e-003; ng(n+1)=6.106937e+002;
n=39; farx(n+1)=8.166961e+000; foe(n+1)=2.328939e+001; krok(n+1)=1.835014e-002; ng(n+1)=1.039282e+003;
n=40; farx(n+1)=8.264813e+000; foe(n+1)=2.226219e+001; krok(n+1)=1.165911e-002; ng(n+1)=1.404243e+003;
n=41; farx(n+1)=7.755596e+000; foe(n+1)=2.070593e+001; krok(n+1)=2.992730e-002; ng(n+1)=1.060580e+003;
n=42; farx(n+1)=7.658762e+000; foe(n+1)=2.029614e+001; krok(n+1)=6.752063e-003; ng(n+1)=5.205809e+002;
n=43; farx(n+1)=7.675787e+000; foe(n+1)=1.951799e+001; krok(n+1)=5.851700e-002; ng(n+1)=4.693539e+002;
n=44; farx(n+1)=7.594413e+000; foe(n+1)=1.895983e+001; krok(n+1)=1.086709e-002; ng(n+1)=2.850200e+002;
n=45; farx(n+1)=6.566936e+000; foe(n+1)=1.605251e+001; krok(n+1)=2.485530e-001; ng(n+1)=5.602271e+002;
n=46; farx(n+1)=6.264413e+000; foe(n+1)=1.545137e+001; krok(n+1)=4.254539e-002; ng(n+1)=2.677438e+002;
n=47; farx(n+1)=5.788973e+000; foe(n+1)=1.445608e+001; krok(n+1)=1.114169e-001; ng(n+1)=3.140996e+002;
n=48; farx(n+1)=5.675723e+000; foe(n+1)=1.383763e+001; krok(n+1)=1.471339e-001; ng(n+1)=2.855169e+002;
n=49; farx(n+1)=5.390828e+000; foe(n+1)=1.333344e+001; krok(n+1)=1.044603e-001; ng(n+1)=4.244450e+002;
n=50; farx(n+1)=5.083540e+000; foe(n+1)=1.244427e+001; krok(n+1)=3.870680e-001; ng(n+1)=6.087031e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=5.112397e+000; foe(n+1)=1.224609e+001; krok(n+1)=1.322977e-005; ng(n+1)=4.647657e+002;
n=52; farx(n+1)=5.094179e+000; foe(n+1)=1.215962e+001; krok(n+1)=8.923578e-005; ng(n+1)=1.066428e+002;
n=53; farx(n+1)=5.124932e+000; foe(n+1)=1.208445e+001; krok(n+1)=4.430569e-004; ng(n+1)=4.866423e+001;
n=54; farx(n+1)=5.075612e+000; foe(n+1)=1.187290e+001; krok(n+1)=7.209870e-004; ng(n+1)=6.110422e+001;
n=55; farx(n+1)=4.974900e+000; foe(n+1)=1.157833e+001; krok(n+1)=1.704775e-003; ng(n+1)=6.836787e+001;
n=56; farx(n+1)=4.858124e+000; foe(n+1)=1.129924e+001; krok(n+1)=4.002465e-003; ng(n+1)=3.525641e+001;
n=57; farx(n+1)=4.806723e+000; foe(n+1)=1.099418e+001; krok(n+1)=3.816075e-003; ng(n+1)=1.309418e+002;
n=58; farx(n+1)=4.749196e+000; foe(n+1)=1.019651e+001; krok(n+1)=2.812841e-002; ng(n+1)=2.680189e+002;
n=59; farx(n+1)=4.407300e+000; foe(n+1)=9.799732e+000; krok(n+1)=1.611092e-002; ng(n+1)=1.924177e+002;
n=60; farx(n+1)=3.876106e+000; foe(n+1)=8.794332e+000; krok(n+1)=3.394563e-002; ng(n+1)=3.283367e+002;
n=61; farx(n+1)=3.833146e+000; foe(n+1)=8.697018e+000; krok(n+1)=1.353531e-003; ng(n+1)=1.824370e+002;
n=62; farx(n+1)=3.750359e+000; foe(n+1)=8.475743e+000; krok(n+1)=4.625835e-003; ng(n+1)=2.856991e+002;
n=63; farx(n+1)=3.556267e+000; foe(n+1)=8.302254e+000; krok(n+1)=1.980819e-002; ng(n+1)=1.001133e+002;
n=64; farx(n+1)=3.391752e+000; foe(n+1)=8.104393e+000; krok(n+1)=2.284436e-002; ng(n+1)=1.328416e+002;
n=65; farx(n+1)=3.231624e+000; foe(n+1)=7.973943e+000; krok(n+1)=4.954706e-002; ng(n+1)=2.273731e+002;
n=66; farx(n+1)=3.040405e+000; foe(n+1)=7.679887e+000; krok(n+1)=3.320390e-002; ng(n+1)=2.549571e+002;
n=67; farx(n+1)=2.760553e+000; foe(n+1)=7.108688e+000; krok(n+1)=1.308949e-001; ng(n+1)=2.335368e+002;
n=68; farx(n+1)=2.631537e+000; foe(n+1)=6.838265e+000; krok(n+1)=2.750088e-002; ng(n+1)=1.558303e+002;
n=69; farx(n+1)=2.552273e+000; foe(n+1)=6.581184e+000; krok(n+1)=2.243330e-002; ng(n+1)=3.221229e+002;
n=70; farx(n+1)=2.474336e+000; foe(n+1)=6.394869e+000; krok(n+1)=3.001959e-002; ng(n+1)=9.355918e+001;
n=71; farx(n+1)=2.252338e+000; foe(n+1)=5.925935e+000; krok(n+1)=1.028419e-001; ng(n+1)=1.978757e+002;
n=72; farx(n+1)=2.185250e+000; foe(n+1)=5.653694e+000; krok(n+1)=3.269900e-002; ng(n+1)=1.066803e+002;
n=73; farx(n+1)=2.175338e+000; foe(n+1)=5.518020e+000; krok(n+1)=7.274432e-002; ng(n+1)=2.388093e+002;
n=74; farx(n+1)=2.149701e+000; foe(n+1)=5.347347e+000; krok(n+1)=7.727538e-002; ng(n+1)=1.147257e+002;
n=75; farx(n+1)=2.012608e+000; foe(n+1)=4.894452e+000; krok(n+1)=1.101515e-001; ng(n+1)=5.748072e+001;
%odnowa zmiennej metryki
n=76; farx(n+1)=2.001356e+000; foe(n+1)=4.860266e+000; krok(n+1)=1.751338e-005; ng(n+1)=1.766639e+002;
n=77; farx(n+1)=1.997769e+000; foe(n+1)=4.826206e+000; krok(n+1)=1.976919e-005; ng(n+1)=1.697053e+002;
n=78; farx(n+1)=1.987661e+000; foe(n+1)=4.806262e+000; krok(n+1)=6.648882e-005; ng(n+1)=6.284712e+001;
n=79; farx(n+1)=1.976569e+000; foe(n+1)=4.735667e+000; krok(n+1)=4.588786e-004; ng(n+1)=5.427177e+001;
n=80; farx(n+1)=1.961810e+000; foe(n+1)=4.696522e+000; krok(n+1)=6.750897e-004; ng(n+1)=3.454364e+001;
n=81; farx(n+1)=1.970113e+000; foe(n+1)=4.572265e+000; krok(n+1)=4.736926e-003; ng(n+1)=2.470347e+001;
n=82; farx(n+1)=1.931596e+000; foe(n+1)=4.415023e+000; krok(n+1)=6.802533e-003; ng(n+1)=4.601251e+001;
n=83; farx(n+1)=1.924171e+000; foe(n+1)=4.373375e+000; krok(n+1)=4.664313e-003; ng(n+1)=1.769807e+002;
n=84; farx(n+1)=1.913010e+000; foe(n+1)=4.335595e+000; krok(n+1)=1.477175e-002; ng(n+1)=2.505230e+002;
n=85; farx(n+1)=1.889636e+000; foe(n+1)=4.214815e+000; krok(n+1)=2.959126e-002; ng(n+1)=3.183439e+002;
n=86; farx(n+1)=1.907961e+000; foe(n+1)=4.111424e+000; krok(n+1)=6.970444e-003; ng(n+1)=4.009187e+002;
n=87; farx(n+1)=1.887067e+000; foe(n+1)=4.056749e+000; krok(n+1)=6.947618e-003; ng(n+1)=1.674041e+002;
n=88; farx(n+1)=1.851375e+000; foe(n+1)=3.968662e+000; krok(n+1)=6.669234e-003; ng(n+1)=6.533574e+001;
n=89; farx(n+1)=1.819442e+000; foe(n+1)=3.882152e+000; krok(n+1)=2.419175e-002; ng(n+1)=1.201623e+002;
n=90; farx(n+1)=1.813763e+000; foe(n+1)=3.831137e+000; krok(n+1)=2.219272e-002; ng(n+1)=1.113898e+002;
n=91; farx(n+1)=1.793334e+000; foe(n+1)=3.797925e+000; krok(n+1)=1.895365e-002; ng(n+1)=1.132473e+002;
n=92; farx(n+1)=1.771974e+000; foe(n+1)=3.742518e+000; krok(n+1)=5.908701e-002; ng(n+1)=1.193542e+002;
n=93; farx(n+1)=1.746638e+000; foe(n+1)=3.643592e+000; krok(n+1)=3.661089e-002; ng(n+1)=1.067295e+002;
n=94; farx(n+1)=1.718778e+000; foe(n+1)=3.444794e+000; krok(n+1)=9.790078e-002; ng(n+1)=1.766478e+002;
n=95; farx(n+1)=1.703698e+000; foe(n+1)=3.374775e+000; krok(n+1)=9.228634e-002; ng(n+1)=1.345002e+002;
n=96; farx(n+1)=1.697763e+000; foe(n+1)=3.278887e+000; krok(n+1)=8.472165e-002; ng(n+1)=2.931806e+002;
n=97; farx(n+1)=1.730190e+000; foe(n+1)=3.212108e+000; krok(n+1)=1.541897e-001; ng(n+1)=1.517556e+002;
n=98; farx(n+1)=1.703926e+000; foe(n+1)=3.150313e+000; krok(n+1)=1.347359e-001; ng(n+1)=1.188826e+002;
n=99; farx(n+1)=1.685818e+000; foe(n+1)=3.037719e+000; krok(n+1)=1.668649e-001; ng(n+1)=1.693328e+002;
n=100; farx(n+1)=1.590329e+000; foe(n+1)=2.953382e+000; krok(n+1)=4.065930e-001; ng(n+1)=1.209444e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
