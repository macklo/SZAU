%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.877836e+003; foe(n+1)=4.915473e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.943751e+003; foe(n+1)=3.928509e+003; krok(n+1)=3.293253e-004; ng(n+1)=2.311741e+003;
n=2; farx(n+1)=3.841597e+003; foe(n+1)=3.832903e+003; krok(n+1)=1.153579e-002; ng(n+1)=1.237668e+002;
n=3; farx(n+1)=3.660237e+003; foe(n+1)=3.650818e+003; krok(n+1)=4.219311e-005; ng(n+1)=2.431384e+003;
n=4; farx(n+1)=3.656820e+003; foe(n+1)=3.595113e+003; krok(n+1)=6.746361e-003; ng(n+1)=3.821814e+003;
n=5; farx(n+1)=2.215424e+003; foe(n+1)=2.556802e+003; krok(n+1)=1.074785e+000; ng(n+1)=4.009052e+003;
n=6; farx(n+1)=2.167214e+003; foe(n+1)=2.527381e+003; krok(n+1)=3.104538e-002; ng(n+1)=1.772829e+003;
n=7; farx(n+1)=1.095549e+003; foe(n+1)=1.976862e+003; krok(n+1)=3.334142e-002; ng(n+1)=1.080409e+004;
n=8; farx(n+1)=8.231678e+002; foe(n+1)=9.328491e+002; krok(n+1)=2.537142e+001; ng(n+1)=1.346026e+003;
n=9; farx(n+1)=8.189369e+002; foe(n+1)=9.264039e+002; krok(n+1)=7.609550e-003; ng(n+1)=8.911250e+002;
n=10; farx(n+1)=7.148799e+002; foe(n+1)=8.385624e+002; krok(n+1)=5.962951e-001; ng(n+1)=8.961713e+002;
n=11; farx(n+1)=6.993246e+002; foe(n+1)=8.227872e+002; krok(n+1)=4.438543e-002; ng(n+1)=8.599348e+003;
n=12; farx(n+1)=6.185702e+002; foe(n+1)=7.815983e+002; krok(n+1)=6.088646e+000; ng(n+1)=3.323525e+002;
n=13; farx(n+1)=6.117291e+002; foe(n+1)=7.749698e+002; krok(n+1)=1.355274e+000; ng(n+1)=3.862404e+003;
n=14; farx(n+1)=5.905742e+002; foe(n+1)=7.673726e+002; krok(n+1)=1.202875e+000; ng(n+1)=3.575267e+002;
n=15; farx(n+1)=5.882855e+002; foe(n+1)=7.645136e+002; krok(n+1)=1.814761e+000; ng(n+1)=1.948548e+003;
n=16; farx(n+1)=5.886623e+002; foe(n+1)=7.638631e+002; krok(n+1)=9.907791e-002; ng(n+1)=1.411162e+003;
n=17; farx(n+1)=5.837112e+002; foe(n+1)=7.632278e+002; krok(n+1)=6.710813e-001; ng(n+1)=3.772960e+002;
n=18; farx(n+1)=5.828747e+002; foe(n+1)=7.630443e+002; krok(n+1)=4.063898e-001; ng(n+1)=6.766571e+001;
n=19; farx(n+1)=5.833477e+002; foe(n+1)=7.629964e+002; krok(n+1)=1.420334e+000; ng(n+1)=4.459352e+002;
n=20; farx(n+1)=5.856764e+002; foe(n+1)=7.628429e+002; krok(n+1)=1.507294e+000; ng(n+1)=7.730907e+002;
n=21; farx(n+1)=5.857902e+002; foe(n+1)=7.628232e+002; krok(n+1)=1.015268e+000; ng(n+1)=7.754621e+002;
n=22; farx(n+1)=5.861499e+002; foe(n+1)=7.628184e+002; krok(n+1)=1.348916e+000; ng(n+1)=3.872914e+001;
n=23; farx(n+1)=5.861761e+002; foe(n+1)=7.628093e+002; krok(n+1)=1.071771e+001; ng(n+1)=1.861425e+001;
n=24; farx(n+1)=5.852815e+002; foe(n+1)=7.626794e+002; krok(n+1)=1.079133e+001; ng(n+1)=1.874124e+002;
n=25; farx(n+1)=5.800590e+002; foe(n+1)=7.598085e+002; krok(n+1)=2.373504e+000; ng(n+1)=8.569386e+002;
%odnowa zmiennej metryki
n=26; farx(n+1)=5.798548e+002; foe(n+1)=7.594513e+002; krok(n+1)=2.846463e-006; ng(n+1)=1.260290e+003;
n=27; farx(n+1)=5.784368e+002; foe(n+1)=7.594114e+002; krok(n+1)=7.553219e-004; ng(n+1)=1.659980e+001;
n=28; farx(n+1)=5.798713e+002; foe(n+1)=7.588300e+002; krok(n+1)=8.648560e-003; ng(n+1)=2.727672e+001;
n=29; farx(n+1)=5.814481e+002; foe(n+1)=7.587278e+002; krok(n+1)=2.563431e-002; ng(n+1)=2.179090e+002;
n=30; farx(n+1)=5.818065e+002; foe(n+1)=7.586562e+002; krok(n+1)=2.367770e-001; ng(n+1)=1.166520e+002;
n=31; farx(n+1)=5.818498e+002; foe(n+1)=7.585869e+002; krok(n+1)=1.098329e+000; ng(n+1)=1.030907e+002;
n=32; farx(n+1)=5.818990e+002; foe(n+1)=7.583917e+002; krok(n+1)=1.402391e+000; ng(n+1)=2.796494e+002;
n=33; farx(n+1)=5.760458e+002; foe(n+1)=7.560790e+002; krok(n+1)=4.585930e+000; ng(n+1)=2.036425e+002;
n=34; farx(n+1)=3.383384e+002; foe(n+1)=6.656391e+002; krok(n+1)=1.799746e+000; ng(n+1)=3.083397e+002;
n=35; farx(n+1)=3.382812e+002; foe(n+1)=6.656388e+002; krok(n+1)=5.724201e-007; ng(n+1)=3.300315e+002;
n=36; farx(n+1)=3.382812e+002; foe(n+1)=6.656388e+002; krok(n+1)=4.877727e-016; ng(n+1)=3.347396e+002;
n=37; farx(n+1)=3.485788e+002; foe(n+1)=5.710209e+002; krok(n+1)=9.891278e-001; ng(n+1)=3.347396e+002;
n=38; farx(n+1)=3.486699e+002; foe(n+1)=5.710063e+002; krok(n+1)=2.686576e-006; ng(n+1)=6.517598e+003;
n=39; farx(n+1)=3.577935e+002; foe(n+1)=5.660986e+002; krok(n+1)=1.638325e-001; ng(n+1)=6.559443e+003;
n=40; farx(n+1)=2.628137e+002; foe(n+1)=5.306940e+002; krok(n+1)=6.595347e-001; ng(n+1)=4.144799e+003;
n=41; farx(n+1)=2.580682e+002; foe(n+1)=5.299318e+002; krok(n+1)=2.977641e-002; ng(n+1)=1.034853e+003;
n=42; farx(n+1)=2.455805e+002; foe(n+1)=5.263750e+002; krok(n+1)=6.024636e-002; ng(n+1)=9.780871e+002;
 % z�y kierunek w metodzie zm - odnowa 
n=43; farx(n+1)=2.350357e+002; foe(n+1)=5.159510e+002; krok(n+1)=1.167393e-005; ng(n+1)=1.415918e+003;
n=44; farx(n+1)=1.850511e+002; foe(n+1)=4.542345e+002; krok(n+1)=1.239407e-004; ng(n+1)=1.504064e+003;
n=45; farx(n+1)=7.333488e+001; foe(n+1)=2.802400e+002; krok(n+1)=3.259080e-004; ng(n+1)=1.338704e+003;
n=46; farx(n+1)=6.445542e+001; foe(n+1)=2.623967e+002; krok(n+1)=4.927090e-004; ng(n+1)=2.467855e+003;
n=47; farx(n+1)=6.166703e+001; foe(n+1)=2.566473e+002; krok(n+1)=8.923578e-005; ng(n+1)=3.430606e+003;
n=48; farx(n+1)=4.723470e+001; foe(n+1)=2.026931e+002; krok(n+1)=2.487263e-003; ng(n+1)=2.765565e+003;
n=49; farx(n+1)=2.746908e+001; foe(n+1)=1.415356e+002; krok(n+1)=7.667796e-003; ng(n+1)=5.840874e+002;
n=50; farx(n+1)=2.978864e+001; foe(n+1)=1.337126e+002; krok(n+1)=7.571876e-001; ng(n+1)=2.772201e+003;
%odnowa zmiennej metryki
n=51; farx(n+1)=2.949445e+001; foe(n+1)=1.317608e+002; krok(n+1)=7.653364e-007; ng(n+1)=3.470198e+003;
n=52; farx(n+1)=2.985647e+001; foe(n+1)=1.286172e+002; krok(n+1)=8.082043e-005; ng(n+1)=4.706652e+002;
n=53; farx(n+1)=2.376643e+001; foe(n+1)=1.223031e+002; krok(n+1)=5.840350e-004; ng(n+1)=4.140538e+002;
n=54; farx(n+1)=2.040780e+001; foe(n+1)=1.169017e+002; krok(n+1)=1.275185e-004; ng(n+1)=9.410612e+002;
n=55; farx(n+1)=1.721464e+001; foe(n+1)=1.039600e+002; krok(n+1)=9.289576e-004; ng(n+1)=2.220116e+003;
n=56; farx(n+1)=1.730162e+001; foe(n+1)=1.039553e+002; krok(n+1)=7.821151e-004; ng(n+1)=1.863166e+003;
n=57; farx(n+1)=1.655807e+001; foe(n+1)=1.018244e+002; krok(n+1)=3.202323e-001; ng(n+1)=1.889727e+003;
n=58; farx(n+1)=1.829688e+001; foe(n+1)=9.711727e+001; krok(n+1)=1.181740e-001; ng(n+1)=2.340707e+003;
n=59; farx(n+1)=1.687247e+001; foe(n+1)=9.364242e+001; krok(n+1)=1.062525e+000; ng(n+1)=2.097575e+003;
n=60; farx(n+1)=1.456994e+001; foe(n+1)=8.785831e+001; krok(n+1)=1.707324e+000; ng(n+1)=2.922996e+003;
n=61; farx(n+1)=1.285980e+001; foe(n+1)=8.398195e+001; krok(n+1)=1.216368e+000; ng(n+1)=1.018224e+003;
n=62; farx(n+1)=1.252565e+001; foe(n+1)=8.274390e+001; krok(n+1)=5.861088e-001; ng(n+1)=1.561225e+003;
n=63; farx(n+1)=1.189166e+001; foe(n+1)=8.196422e+001; krok(n+1)=4.164005e-001; ng(n+1)=9.079179e+002;
n=64; farx(n+1)=1.131037e+001; foe(n+1)=8.171881e+001; krok(n+1)=4.997929e-001; ng(n+1)=1.279534e+002;
n=65; farx(n+1)=1.122522e+001; foe(n+1)=8.161868e+001; krok(n+1)=4.756263e-001; ng(n+1)=2.853614e+002;
n=66; farx(n+1)=1.122700e+001; foe(n+1)=8.159335e+001; krok(n+1)=1.080691e+000; ng(n+1)=1.874165e+001;
n=67; farx(n+1)=1.121385e+001; foe(n+1)=8.158633e+001; krok(n+1)=1.691553e+000; ng(n+1)=4.886243e+001;
n=68; farx(n+1)=1.123043e+001; foe(n+1)=8.158393e+001; krok(n+1)=1.447680e+000; ng(n+1)=4.756057e+001;
n=69; farx(n+1)=1.124349e+001; foe(n+1)=8.158374e+001; krok(n+1)=1.348916e+000; ng(n+1)=3.695797e+000;
n=70; farx(n+1)=1.126261e+001; foe(n+1)=8.158316e+001; krok(n+1)=9.108337e+000; ng(n+1)=5.183221e+000;
n=71; farx(n+1)=1.129411e+001; foe(n+1)=8.157580e+001; krok(n+1)=1.048528e+001; ng(n+1)=9.061796e+000;
n=72; farx(n+1)=1.130489e+001; foe(n+1)=8.157480e+001; krok(n+1)=1.095445e+000; ng(n+1)=3.140992e+001;
n=73; farx(n+1)=1.130407e+001; foe(n+1)=8.157480e+001; krok(n+1)=7.807902e-001; ng(n+1)=6.975398e-001;
n=74; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=8.752990e-001; ng(n+1)=6.716849e-002;
n=75; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=1.284145e-005; ng(n+1)=3.884989e-003;
%odnowa zmiennej metryki
n=76; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=3.052279e-005; ng(n+1)=3.884939e-003;
n=77; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=1.698321e-006; ng(n+1)=3.378940e-003;
n=78; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=5.162212e-005; ng(n+1)=2.090932e-003;
n=79; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=9.989885e-004; ng(n+1)=4.956843e-004;
n=80; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=2.244025e-005; ng(n+1)=3.278265e-004;
n=81; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=3.052036e-007; ng(n+1)=3.275385e-004;
n=82; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=2.734729e-005; ng(n+1)=3.275377e-004;
n=83; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=2.371798e-005; ng(n+1)=3.275150e-004;
n=84; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=7.150304e-006; ng(n+1)=3.275059e-004;
n=85; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=3.262270e-006; ng(n+1)=3.275036e-004;
 % z�y kierunek w metodzie zm - odnowa 
n=86; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=1.870933e-011; ng(n+1)=3.275025e-004;
n=87; farx(n+1)=1.130425e+001; foe(n+1)=8.157480e+001; krok(n+1)=2.407322e-005; ng(n+1)=3.275012e-004;
n=88; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=1.100811e-003; ng(n+1)=6.063877e-004;
n=89; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=2.288097e-004; ng(n+1)=7.411635e-004;
n=90; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=2.274910e-003; ng(n+1)=4.626558e-004;
n=91; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=6.569632e-005; ng(n+1)=2.588902e-004;
n=92; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=1.171586e-005; ng(n+1)=2.572850e-004;
n=93; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=4.004915e-006; ng(n+1)=2.569660e-004;
n=94; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=7.579439e-006; ng(n+1)=2.569650e-004;
n=95; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=2.497580e-007; ng(n+1)=2.569631e-004;
n=96; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=5.237933e-006; ng(n+1)=2.569631e-004;
n=97; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=1.430061e-005; ng(n+1)=2.569617e-004;
n=98; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=2.838071e-007; ng(n+1)=2.569580e-004;
n=99; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=2.493323e-009; ng(n+1)=2.569580e-004;
n=100; farx(n+1)=1.130424e+001; foe(n+1)=8.157480e+001; krok(n+1)=4.749515e-011; ng(n+1)=2.569579e-004;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
