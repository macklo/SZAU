%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.048345e+003; foe(n+1)=4.046763e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.959553e+003; foe(n+1)=3.023792e+003; krok(n+1)=5.005982e-004; ng(n+1)=1.556019e+003;
n=2; farx(n+1)=6.084544e+002; foe(n+1)=8.148840e+002; krok(n+1)=9.280118e-003; ng(n+1)=3.408420e+002;
n=3; farx(n+1)=6.186350e+002; foe(n+1)=7.811160e+002; krok(n+1)=1.246120e-005; ng(n+1)=3.750844e+003;
n=4; farx(n+1)=1.403810e+003; foe(n+1)=7.506085e+002; krok(n+1)=5.357205e-003; ng(n+1)=3.893639e+003;
n=5; farx(n+1)=1.024216e+003; foe(n+1)=7.212202e+002; krok(n+1)=4.069805e-002; ng(n+1)=2.135529e+003;
n=6; farx(n+1)=1.034752e+003; foe(n+1)=7.142189e+002; krok(n+1)=2.452306e-002; ng(n+1)=1.916117e+003;
n=7; farx(n+1)=8.181211e+002; foe(n+1)=6.106536e+002; krok(n+1)=8.877087e-002; ng(n+1)=1.994485e+003;
n=8; farx(n+1)=4.552539e+002; foe(n+1)=5.410120e+002; krok(n+1)=3.607652e-002; ng(n+1)=1.508037e+003;
n=9; farx(n+1)=2.615426e+002; foe(n+1)=4.883443e+002; krok(n+1)=1.536047e-002; ng(n+1)=2.164931e+003;
n=10; farx(n+1)=2.055071e+002; foe(n+1)=4.158062e+002; krok(n+1)=3.607908e-001; ng(n+1)=4.313446e+003;
n=11; farx(n+1)=9.297490e+001; foe(n+1)=2.374880e+002; krok(n+1)=1.454604e+000; ng(n+1)=2.744330e+003;
n=12; farx(n+1)=6.710498e+001; foe(n+1)=2.098382e+002; krok(n+1)=2.373843e-001; ng(n+1)=2.572540e+003;
n=13; farx(n+1)=4.339161e+001; foe(n+1)=1.638715e+002; krok(n+1)=1.728699e-001; ng(n+1)=3.660962e+003;
n=14; farx(n+1)=3.588343e+001; foe(n+1)=1.395205e+002; krok(n+1)=1.350958e-001; ng(n+1)=2.717249e+003;
n=15; farx(n+1)=2.815650e+001; foe(n+1)=1.189248e+002; krok(n+1)=8.457767e-001; ng(n+1)=7.957349e+002;
n=16; farx(n+1)=1.743049e+001; foe(n+1)=1.096999e+002; krok(n+1)=5.636152e-001; ng(n+1)=5.098442e+002;
n=17; farx(n+1)=1.601099e+001; foe(n+1)=1.055121e+002; krok(n+1)=1.948045e+000; ng(n+1)=4.044446e+002;
n=18; farx(n+1)=1.382459e+001; foe(n+1)=1.009772e+002; krok(n+1)=1.905690e+000; ng(n+1)=1.473875e+002;
n=19; farx(n+1)=1.376550e+001; foe(n+1)=9.809957e+001; krok(n+1)=3.456460e-001; ng(n+1)=2.640488e+002;
n=20; farx(n+1)=1.105591e+001; foe(n+1)=9.011877e+001; krok(n+1)=1.422745e+000; ng(n+1)=8.911279e+002;
n=21; farx(n+1)=1.062727e+001; foe(n+1)=8.881401e+001; krok(n+1)=3.261226e-001; ng(n+1)=9.329438e+002;
n=22; farx(n+1)=9.914363e+000; foe(n+1)=8.792707e+001; krok(n+1)=2.429623e-001; ng(n+1)=7.168773e+002;
n=23; farx(n+1)=9.735171e+000; foe(n+1)=8.704802e+001; krok(n+1)=3.785765e-001; ng(n+1)=2.359646e+002;
n=24; farx(n+1)=9.605652e+000; foe(n+1)=8.639393e+001; krok(n+1)=4.787134e-001; ng(n+1)=2.815958e+002;
n=25; farx(n+1)=9.874630e+000; foe(n+1)=8.494487e+001; krok(n+1)=1.104871e+000; ng(n+1)=5.725735e+002;
%odnowa zmiennej metryki
n=26; farx(n+1)=9.725118e+000; foe(n+1)=8.466210e+001; krok(n+1)=4.269784e-006; ng(n+1)=5.852877e+002;
n=27; farx(n+1)=9.666694e+000; foe(n+1)=8.459780e+001; krok(n+1)=2.288097e-004; ng(n+1)=3.982067e+001;
n=28; farx(n+1)=9.741469e+000; foe(n+1)=8.438478e+001; krok(n+1)=1.116407e-004; ng(n+1)=7.819368e+001;
n=29; farx(n+1)=1.003392e+001; foe(n+1)=8.398431e+001; krok(n+1)=6.847610e-004; ng(n+1)=8.414674e+001;
n=30; farx(n+1)=1.000374e+001; foe(n+1)=8.394546e+001; krok(n+1)=2.189959e-004; ng(n+1)=3.993261e+001;
n=31; farx(n+1)=9.975542e+000; foe(n+1)=8.377376e+001; krok(n+1)=7.368272e-004; ng(n+1)=8.649652e+001;
n=32; farx(n+1)=9.984994e+000; foe(n+1)=8.340769e+001; krok(n+1)=5.654998e-003; ng(n+1)=1.090158e+002;
n=33; farx(n+1)=9.954444e+000; foe(n+1)=8.319725e+001; krok(n+1)=9.282331e-001; ng(n+1)=1.949060e+002;
n=34; farx(n+1)=1.002668e+001; foe(n+1)=8.306823e+001; krok(n+1)=9.961132e-001; ng(n+1)=1.146926e+002;
n=35; farx(n+1)=1.001301e+001; foe(n+1)=8.305056e+001; krok(n+1)=7.215816e-001; ng(n+1)=6.196124e+001;
n=36; farx(n+1)=1.000711e+001; foe(n+1)=8.304681e+001; krok(n+1)=1.183754e+000; ng(n+1)=2.718272e+001;
n=37; farx(n+1)=1.001141e+001; foe(n+1)=8.304654e+001; krok(n+1)=7.310195e-001; ng(n+1)=1.201567e+001;
n=38; farx(n+1)=1.000770e+001; foe(n+1)=8.304649e+001; krok(n+1)=9.857126e-001; ng(n+1)=4.854708e+000;
n=39; farx(n+1)=1.000488e+001; foe(n+1)=8.304647e+001; krok(n+1)=2.684325e+000; ng(n+1)=1.657060e+000;
n=40; farx(n+1)=1.000498e+001; foe(n+1)=8.304647e+001; krok(n+1)=1.585247e+000; ng(n+1)=2.878110e-001;
n=41; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.492366e+000; ng(n+1)=3.026095e-001;
n=42; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.420334e+000; ng(n+1)=9.282719e-002;
n=43; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=8.579941e-006; ng(n+1)=1.140957e-002;
n=44; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=9.888442e-008; ng(n+1)=1.140947e-002;
n=45; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.802903e-009; ng(n+1)=1.140947e-002;
n=46; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=4.868077e-007; ng(n+1)=1.140947e-002;
n=47; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=3.629487e-006; ng(n+1)=1.140947e-002;
n=48; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=7.243029e-006; ng(n+1)=1.140942e-002;
n=49; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.483310e-005; ng(n+1)=1.140934e-002;
n=50; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.111900e-008; ng(n+1)=1.140917e-002;
%odnowa zmiennej metryki
n=51; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.141062e-005; ng(n+1)=1.140917e-002;
n=52; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=6.273389e-006; ng(n+1)=1.085037e-002;
n=53; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.253881e-004; ng(n+1)=1.611778e-003;
n=54; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=3.870231e-004; ng(n+1)=3.745053e-004;
n=55; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.504753e-005; ng(n+1)=3.534966e-004;
n=56; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.209447e-005; ng(n+1)=3.528893e-004;
n=57; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=3.616520e-005; ng(n+1)=3.527956e-004;
n=58; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=6.010240e-006; ng(n+1)=3.531007e-004;
n=59; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=3.320428e-006; ng(n+1)=3.530986e-004;
n=60; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.523663e-005; ng(n+1)=3.530974e-004;
n=61; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=8.600086e-006; ng(n+1)=3.530920e-004;
n=62; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=7.353122e-011; ng(n+1)=3.530890e-004;
n=63; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=6.333551e-009; ng(n+1)=3.530889e-004;
n=64; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=4.468940e-007; ng(n+1)=3.530889e-004;
n=65; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=5.860955e-007; ng(n+1)=3.530944e-004;
n=66; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=5.586175e-008; ng(n+1)=3.530936e-004;
 % z�y kierunek w metodzie zm - odnowa 
n=67; farx(n+1)=1.000504e+001; foe(n+1)=8.304646e+001; krok(n+1)=5.671074e-008; ng(n+1)=3.530936e-004;
n=68; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=2.603533e-004; ng(n+1)=3.501485e-004;
n=69; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=5.720244e-005; ng(n+1)=5.209971e-004;
n=70; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.215721e-003; ng(n+1)=5.002681e-004;
n=71; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=2.279497e-004; ng(n+1)=5.681372e-004;
n=72; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=8.637483e-006; ng(n+1)=4.371217e-004;
n=73; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=3.022327e-005; ng(n+1)=4.370047e-004;
n=74; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.069939e-005; ng(n+1)=4.366391e-004;
n=75; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=3.213587e-005; ng(n+1)=4.366344e-004;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=5.500912e-006; ng(n+1)=4.366203e-004;
n=77; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.281243e-007; ng(n+1)=3.424106e-004;
n=78; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=8.545717e-004; ng(n+1)=3.419313e-004;
n=79; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=2.054632e-004; ng(n+1)=5.118981e-004;
n=80; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.353531e-003; ng(n+1)=4.002887e-004;
n=81; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=2.617323e-003; ng(n+1)=4.772986e-004;
n=82; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.592817e-005; ng(n+1)=4.081689e-004;
n=83; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=1.149157e-005; ng(n+1)=4.081403e-004;
n=84; farx(n+1)=1.000505e+001; foe(n+1)=8.304646e+001; krok(n+1)=2.565543e-008; ng(n+1)=4.081357e-004;
 % z�y kierunek w metodzie zm - odnowa 
 % z�y kierunek w metodzie zm po wykonaniu odnowy

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
