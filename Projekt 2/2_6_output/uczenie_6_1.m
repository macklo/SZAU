%uczenie predyktora arx
clear all;
n=0; farx(n+1)=4.265830e+003; foe(n+1)=4.225145e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.928137e+003; foe(n+1)=3.018736e+003; krok(n+1)=4.927090e-004; ng(n+1)=2.792734e+003;
n=2; farx(n+1)=2.160149e+003; foe(n+1)=4.868063e+003; krok(n+1)=7.988163e-004; ng(n+1)=1.202098e+003;
n=3; farx(n+1)=1.396011e+003; foe(n+1)=1.073884e+004; krok(n+1)=3.824295e-004; ng(n+1)=3.772536e+003;
n=4; farx(n+1)=4.878038e+002; foe(n+1)=2.101018e+004; krok(n+1)=8.487097e-003; ng(n+1)=2.585309e+003;
n=5; farx(n+1)=4.864913e+002; foe(n+1)=2.100950e+004; krok(n+1)=4.924301e-005; ng(n+1)=1.735954e+003;
n=6; farx(n+1)=2.212116e+002; foe(n+1)=2.415532e+004; krok(n+1)=5.636152e-001; ng(n+1)=1.641495e+003;
n=7; farx(n+1)=1.569762e+002; foe(n+1)=2.376684e+004; krok(n+1)=4.108715e-001; ng(n+1)=7.819961e+002;
n=8; farx(n+1)=8.176890e+001; foe(n+1)=6.203480e+002; krok(n+1)=8.892951e-001; ng(n+1)=2.551277e+002;
n=9; farx(n+1)=6.389924e+001; foe(n+1)=1.199739e+004; krok(n+1)=3.854790e-001; ng(n+1)=3.354088e+002;
n=10; farx(n+1)=5.277485e+001; foe(n+1)=1.142670e+004; krok(n+1)=7.763546e-001; ng(n+1)=2.043684e+002;
n=11; farx(n+1)=3.822535e+001; foe(n+1)=3.472880e+003; krok(n+1)=7.741360e-001; ng(n+1)=1.591678e+002;
n=12; farx(n+1)=2.254798e+001; foe(n+1)=2.896088e+002; krok(n+1)=1.458994e+000; ng(n+1)=1.545049e+002;
n=13; farx(n+1)=9.132679e+000; foe(n+1)=1.969901e+002; krok(n+1)=1.307567e+000; ng(n+1)=4.153969e+002;
n=14; farx(n+1)=5.588812e+000; foe(n+1)=1.293199e+002; krok(n+1)=7.807902e-001; ng(n+1)=7.528236e+001;
n=15; farx(n+1)=4.731026e+000; foe(n+1)=4.047621e+002; krok(n+1)=1.601772e+000; ng(n+1)=7.848802e+001;
n=16; farx(n+1)=4.178615e+000; foe(n+1)=1.750178e+002; krok(n+1)=8.131859e-001; ng(n+1)=4.015726e+001;
n=17; farx(n+1)=4.120110e+000; foe(n+1)=1.684719e+002; krok(n+1)=7.021280e-001; ng(n+1)=4.130454e+001;
n=18; farx(n+1)=4.086083e+000; foe(n+1)=1.991844e+002; krok(n+1)=1.180120e+000; ng(n+1)=2.251540e+001;
n=19; farx(n+1)=4.080560e+000; foe(n+1)=1.968448e+002; krok(n+1)=1.009068e+000; ng(n+1)=1.021826e+001;
n=20; farx(n+1)=4.058881e+000; foe(n+1)=1.874757e+002; krok(n+1)=6.124441e+000; ng(n+1)=9.052898e+000;
n=21; farx(n+1)=3.886563e+000; foe(n+1)=1.614416e+002; krok(n+1)=7.344400e+000; ng(n+1)=1.997980e+001;
n=22; farx(n+1)=3.644101e+000; foe(n+1)=1.599996e+002; krok(n+1)=2.485026e+000; ng(n+1)=9.574595e+001;
n=23; farx(n+1)=3.445773e+000; foe(n+1)=1.584828e+002; krok(n+1)=9.270274e-001; ng(n+1)=1.130247e+002;
n=24; farx(n+1)=3.269351e+000; foe(n+1)=1.319172e+002; krok(n+1)=1.049915e+000; ng(n+1)=6.936120e+000;
n=25; farx(n+1)=3.171258e+000; foe(n+1)=1.297180e+002; krok(n+1)=1.625559e+000; ng(n+1)=4.789540e+001;
%odnowa zmiennej metryki
n=26; farx(n+1)=3.111498e+000; foe(n+1)=1.373129e+002; krok(n+1)=3.502676e-005; ng(n+1)=8.585786e+001;
n=27; farx(n+1)=3.093174e+000; foe(n+1)=1.383492e+002; krok(n+1)=1.733806e-004; ng(n+1)=2.090081e+001;
n=28; farx(n+1)=3.088006e+000; foe(n+1)=1.324169e+002; krok(n+1)=2.914397e-003; ng(n+1)=3.225816e+000;
n=29; farx(n+1)=3.006512e+000; foe(n+1)=1.497452e+002; krok(n+1)=3.481722e-002; ng(n+1)=3.539081e+000;
n=30; farx(n+1)=2.964577e+000; foe(n+1)=1.400567e+002; krok(n+1)=5.753341e-002; ng(n+1)=1.423571e+001;
n=31; farx(n+1)=2.940235e+000; foe(n+1)=1.262982e+002; krok(n+1)=3.768236e-001; ng(n+1)=4.545147e+001;
n=32; farx(n+1)=2.914805e+000; foe(n+1)=1.251795e+002; krok(n+1)=1.677326e-001; ng(n+1)=5.721946e+001;
n=33; farx(n+1)=2.888425e+000; foe(n+1)=1.290822e+002; krok(n+1)=1.011873e-001; ng(n+1)=4.882867e+001;
n=34; farx(n+1)=2.809346e+000; foe(n+1)=1.405006e+002; krok(n+1)=1.420334e+000; ng(n+1)=5.183359e+001;
n=35; farx(n+1)=2.789213e+000; foe(n+1)=1.297629e+002; krok(n+1)=5.596492e-001; ng(n+1)=2.690082e+001;
n=36; farx(n+1)=2.759513e+000; foe(n+1)=1.186719e+002; krok(n+1)=5.919790e-001; ng(n+1)=3.731668e+001;
n=37; farx(n+1)=2.713539e+000; foe(n+1)=1.304295e+002; krok(n+1)=2.150960e+000; ng(n+1)=4.814244e+001;
n=38; farx(n+1)=2.697024e+000; foe(n+1)=1.245231e+002; krok(n+1)=9.282331e-001; ng(n+1)=1.729022e+001;
n=39; farx(n+1)=2.686917e+000; foe(n+1)=1.251212e+002; krok(n+1)=2.845227e-001; ng(n+1)=3.368939e+001;
n=40; farx(n+1)=2.673056e+000; foe(n+1)=1.162949e+002; krok(n+1)=1.050152e+000; ng(n+1)=2.763620e+001;
n=41; farx(n+1)=2.662769e+000; foe(n+1)=1.215213e+002; krok(n+1)=1.164915e+000; ng(n+1)=2.999155e+001;
n=42; farx(n+1)=2.655454e+000; foe(n+1)=1.195344e+002; krok(n+1)=1.044632e+000; ng(n+1)=1.063507e+001;
n=43; farx(n+1)=2.650624e+000; foe(n+1)=1.203771e+002; krok(n+1)=5.014424e-001; ng(n+1)=2.449338e+001;
n=44; farx(n+1)=2.645495e+000; foe(n+1)=1.169798e+002; krok(n+1)=5.714050e-001; ng(n+1)=1.571346e+001;
n=45; farx(n+1)=2.641208e+000; foe(n+1)=1.149973e+002; krok(n+1)=1.356032e+000; ng(n+1)=2.416973e+001;
n=46; farx(n+1)=2.638604e+000; foe(n+1)=1.179583e+002; krok(n+1)=9.155235e-001; ng(n+1)=1.577107e+001;
n=47; farx(n+1)=2.635316e+000; foe(n+1)=1.160570e+002; krok(n+1)=7.679697e-001; ng(n+1)=9.092569e+000;
n=48; farx(n+1)=2.633715e+000; foe(n+1)=1.139559e+002; krok(n+1)=1.095445e+000; ng(n+1)=1.532791e+001;
n=49; farx(n+1)=2.632300e+000; foe(n+1)=1.138256e+002; krok(n+1)=8.252141e-001; ng(n+1)=1.862420e+001;
n=50; farx(n+1)=2.630707e+000; foe(n+1)=1.127279e+002; krok(n+1)=1.285695e+000; ng(n+1)=2.734670e+000;
%odnowa zmiennej metryki
n=51; farx(n+1)=2.630461e+000; foe(n+1)=1.132507e+002; krok(n+1)=1.272771e-005; ng(n+1)=8.627786e+000;
n=52; farx(n+1)=2.630389e+000; foe(n+1)=1.134258e+002; krok(n+1)=7.560801e-005; ng(n+1)=1.555786e+000;
n=53; farx(n+1)=2.630373e+000; foe(n+1)=1.139823e+002; krok(n+1)=1.564230e-003; ng(n+1)=2.305365e-001;
n=54; farx(n+1)=2.630027e+000; foe(n+1)=1.135997e+002; krok(n+1)=1.843937e-002; ng(n+1)=2.929786e-001;
n=55; farx(n+1)=2.629820e+000; foe(n+1)=1.134347e+002; krok(n+1)=3.951655e-002; ng(n+1)=1.611140e-001;
n=56; farx(n+1)=2.629661e+000; foe(n+1)=1.129322e+002; krok(n+1)=8.277891e-002; ng(n+1)=1.556766e-001;
n=57; farx(n+1)=2.629115e+000; foe(n+1)=1.119174e+002; krok(n+1)=3.414648e+000; ng(n+1)=3.657234e-001;
n=58; farx(n+1)=2.628873e+000; foe(n+1)=1.119964e+002; krok(n+1)=1.872544e+000; ng(n+1)=6.901673e+000;
n=59; farx(n+1)=2.628795e+000; foe(n+1)=1.117931e+002; krok(n+1)=1.230345e+000; ng(n+1)=1.170257e+000;
n=60; farx(n+1)=2.628786e+000; foe(n+1)=1.116750e+002; krok(n+1)=1.126364e+000; ng(n+1)=8.771931e-001;
n=61; farx(n+1)=2.628781e+000; foe(n+1)=1.117495e+002; krok(n+1)=4.641165e-001; ng(n+1)=9.108170e-001;
n=62; farx(n+1)=2.628781e+000; foe(n+1)=1.117309e+002; krok(n+1)=1.150553e+000; ng(n+1)=1.851859e-001;
n=63; farx(n+1)=2.628781e+000; foe(n+1)=1.117337e+002; krok(n+1)=1.462039e+000; ng(n+1)=9.678237e-003;
n=64; farx(n+1)=2.628781e+000; foe(n+1)=1.117340e+002; krok(n+1)=1.106940e+000; ng(n+1)=8.672903e-003;
n=65; farx(n+1)=2.628781e+000; foe(n+1)=1.117340e+002; krok(n+1)=1.646043e-005; ng(n+1)=5.710480e-004;
n=66; farx(n+1)=2.628781e+000; foe(n+1)=1.117340e+002; krok(n+1)=3.688974e-006; ng(n+1)=5.710382e-004;
n=67; farx(n+1)=2.628781e+000; foe(n+1)=1.117340e+002; krok(n+1)=6.099114e-009; ng(n+1)=5.710361e-004;
n=68; farx(n+1)=2.628781e+000; foe(n+1)=1.117340e+002; krok(n+1)=2.452157e-009; ng(n+1)=5.710361e-004;
 % z造 kierunek w metodzie zm - odnowa 
n=69; farx(n+1)=2.628781e+000; foe(n+1)=1.117340e+002; krok(n+1)=1.356750e-005; ng(n+1)=5.710361e-004;
n=70; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=5.611139e-005; ng(n+1)=1.969363e-004;
n=71; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=2.099864e-011; ng(n+1)=2.511613e-005;
 % z造 kierunek w metodzie zm - odnowa 
n=72; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=7.611819e-009; ng(n+1)=2.511613e-005;
n=73; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=3.324083e-005; ng(n+1)=2.510365e-005;
n=74; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=5.385865e-006; ng(n+1)=3.547944e-005;
 % z造 kierunek w metodzie zm - odnowa 
n=75; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.383480e-005; ng(n+1)=3.546190e-005;
%odnowa zmiennej metryki
n=76; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=6.844574e-007; ng(n+1)=9.985851e-006;
n=77; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.557473e-005; ng(n+1)=9.625852e-006;
n=78; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=2.547731e-005; ng(n+1)=9.022929e-006;
n=79; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.430061e-005; ng(n+1)=9.186581e-006;
 % z造 kierunek w metodzie zm - odnowa 
n=80; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.290553e-005; ng(n+1)=9.186437e-006;
n=81; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.611240e-005; ng(n+1)=4.828075e-006;
n=82; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=9.917887e-006; ng(n+1)=4.247190e-006;
n=83; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=7.165071e-007; ng(n+1)=4.205271e-006;
n=84; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=3.338291e-008; ng(n+1)=4.205271e-006;
n=85; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.042950e-009; ng(n+1)=4.205271e-006;
n=86; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=5.349696e-006; ng(n+1)=4.205271e-006;
n=87; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=1.065195e-005; ng(n+1)=4.205207e-006;
n=88; farx(n+1)=2.628781e+000; foe(n+1)=1.117339e+002; krok(n+1)=5.214084e-009; ng(n+1)=4.205125e-006;
 % z造 kierunek w metodzie zm - odnowa 
 % z造 kierunek w metodzie zm po wykonaniu odnowy

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora ARX');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
