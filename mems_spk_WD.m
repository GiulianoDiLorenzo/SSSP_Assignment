clear; close all; clc

%% Simscape File for Ground Truth

load('ssc_output.mat')

%% Sampling Frequency
fs = 192e3;

%% Sampling Period
Ts = 1/fs;

%% Simulation Duration
stop_time = 2;  % [seconds]

%% Input Signal

% Time Axis
t = 0:Ts:stop_time;
t = t(1:end-1);

% Signal Amplitude
A = 10;
vin = A * sweeptone(1.99, 0.01, fs, 'SweepFrequencyRange', [500 20000]);

%% Circuit Parameters
% Piezoelectric MEMS loudspeaker in Free-field conditions

% Transformer Ratio
alpha = 3.7e-4;
Seff = 2.0e-5;

% Electrical Domain
Re = 4;
Cp = 2.4e-8;

% Mechanical Domain
Rm = 9.7e-3;
Mm = 1e-6;
Cm = 2.2e-3;

% Acoustic Domain
Cbc = 3.6e-13; 
Ltube1 = 1e2;
Ltube2 = 1e2;
Ctube = 6.5e-13;
Rac = 5e6;

%% Removing Ideal Transformers (Mechanical Domain only)

gamma1 = alpha^-1;
gamma2 = Seff;

% Resistive Elements
R1 = Re/gamma1^2;
R2 = Rm;
R3 = Rac*gamma2^2;


% Dynamic Elements
% Capacitors
C1 = Cp*gamma1^2;
C2 = Cm;
C3 = Cbc/gamma2^2;
C4 = Ctube/gamma2^2;


% Inductors
L1 = Mm;
L2 = Ltube1*gamma2^2;
L3 = Ltube2*gamma2^2;

% Source
Fin = alpha * vin;

%% Setting of Free Parameters (Adaptation Conditions)
%Single elements impedance
Z_R1 = R1;
Z_R2 = R2;
Z_R3 = R3;
Z_C1 = Ts/(2*C1);
Z_C2 = Ts/(2*C2);
Z_C3 = Ts/(2*C3);
Z_C4 = Ts/(2*C4);
Z_L1 = (2*L1)/Ts;
Z_L2 = (2*L2)/Ts;
Z_L3 = (2*L3)/Ts;

%Adaption Conditions
Z2 = Z_R1;
Z5 = Z_C1;
Z8 = Z_R2;
Z11 = Z_L1;
Z14 = Z_C2;
Z17 = Z_C3;
Z20 = Z_L2;
Z23 = Z_C4;
Z26 = Z_L3;
Z27 = Z_R3;

Z25 = Z26 + Z27;
Z24 = Z25;
Z22 = (Z23*Z24)/(Z23 + Z24);
Z21 = Z22;
Z19 = Z20 + Z21;
Z18 = Z19;
Z16 = (Z17*Z18)/(Z17 + Z18);
Z15 = Z16;
Z13 = Z14 + Z15;
Z12 = Z13;
Z10 = Z11 + Z12;
Z9 = Z10;
Z7 = Z8 + Z9;
Z6 = Z7;
Z4 = (Z5*Z6)/(Z5 + Z6);
Z3 = Z4;
Z1 = Z2 + Z3;


%% Computing Scattering Matrices
Bser = [1,1,1];
Bsermin = [-1,1,1];
Qpar = [1,1,1];

Zpar1 = diag([Z4,Z5,Z6]);
Zpar2 = diag([Z16,Z17,Z18]);
Zpar3 = diag([Z22,Z23,Z24]);


Zser1 = diag([Z1,Z2,Z3]);
Zser2 = diag([Z7,Z8,Z9]);
Zser3 = diag([Z10,Z11,Z12]);
Zser4 = diag([Z13,Z14,Z15]);
Zser5 = diag([Z19,Z20,Z21]);
Zser6 = diag([Z25,Z26,Z27]);

Sser1= eye(3) - 2*Zser1*Bser'*inv(Bser*Zser1*Bser')*Bser;
Sser2= eye(3) - 2*Zser2*Bser'*inv(Bser*Zser2*Bser')*Bser;
Sser3= eye(3) - 2*Zser3*Bser'*inv(Bser*Zser3*Bser')*Bser;
%Sser4= eye(3) - 2*Zser4*Bsermin'*inv(Bsermin*Zser4*Bsermin')*Bsermin;
Sser4= eye(3) - 2*Zser4*Bser'*inv(Bser*Zser4*Bser')*Bser;
Sser5= eye(3) - 2*Zser5*Bsermin'*inv(Bsermin*Zser5*Bsermin')*Bsermin;
Sser6= eye(3) - 2*Zser6*Bsermin'*inv(Bsermin*Zser6*Bsermin')*Bsermin;

Spar1= 2*Qpar'*inv(Qpar*inv(Zpar1)*Qpar')*Qpar*inv(Zpar1) - eye(3);
Spar2= 2*Qpar'*inv(Qpar*inv(Zpar2)*Qpar')*Qpar*inv(Zpar2) - eye(3);
Spar3= 2*Qpar'*inv(Qpar*inv(Zpar3)*Qpar')*Qpar*inv(Zpar3) - eye(3);

%% Initialization of Waves
a2 = 0;
a5 = 0; 
a8 = 0; 
a27 = 0;
a11 = 0;
a14 = 0;
a17 = 0;
a20 = 0;
a23 = 0;
a26 = 0;
b5 = 0; 
b11 = 0;
b14 = 0;
b17 = 0; 
b20 = 0;
b23 = 0; 
b26 = 0;

%% Initialization of Output Signals

Fout = zeros(1, length(t));

%% Simulation Algorithm

for n = 1 : length(Fin)
    a5 = b5;
    a11 = -b11;
    a14 = b14;
    a17 = b17;
    a20 = -b20;
    a23 = b23;
    a26 = -b26;

    % Forward Scan
    b25 = Sser6(1,:)*[0;a26;a27];
    a24 = b25;

    b22 = Spar3(1,:)*[0;a23;a24];
    a21 = b22;

    b19 = Sser5(1,:)*[0;a20;a21];
    a18 = b19;

    b16 = Spar2(1,:)*[0;a17;a18];
    a15 = b16;

    b13 = Sser4(1,:)*[0;a14;a15];
    a12 = b13;

    b10 = Sser3(1,:)*[0;a11;a12];
    a9 = b10;
    
    b7 = Sser2(1,:)*[0;a8;a9];
    a6 = b7;

    b4 = Spar1(1,:)*[0;a5;a6];
    a3 = b4;

    b1 = Sser1(1,:)*[0;a2;a3];

    % Local Root Scattering

     a_root = b1;
    b_root = 2*Fin(n) - a_root;

    % Backward Scan
    a1 = b_root;
    b2 = Sser1(2,:)*[a1;a2;a3];
    b3 = Sser1(3,:)*[a1;a2;a3];
    a4 = b3;

    b5 = Spar1(2,:)*[a4;a5;a6];
    b6 = Spar1(3,:)*[a4;a5;a6];
    a7 = b6;

    b8 = Sser2(2,:)*[a7;a8;a9];
    b9 = Sser2(3,:)*[a7;a8;a9];
    a10 = b9;

    b11 = Sser3(2,:)*[a10;a11;a12];
    b12 = Sser3(3,:)*[a10;a11;a12];
    a13 = b12;

    b14 = Sser4(2,:)*[a13;a14;a15];
    b15 = Sser4(3,:)*[a13;a14;a15];
    a16 = b15;

    b17 = Spar2(2,:)*[a16;a17;a18];
    b18 = Spar2(3,:)*[a16;a17;a18];
    a19 = b18;

    b20 = Sser5(2,:)*[a19;a20;a21];
    b21 = Sser5(3,:)*[a19;a20;a21];
    a22 = b21;

    b23 = Spar3(2,:)*[a22;a23;a24];
    b24 = Spar3(3,:)*[a22;a23;a24];
    a25 = b24;

    b26 = Sser6(2,:)*[a25;a26;a27];
    b27 = Sser6(3,:)*[a25;a26;a27];

    % Read Output
    Fout(n) = (a27+b27)/2;

end

%% Output Plots

% Computing acoustic pressure
pout = Fout ./ Seff;

% Time Domain Plots
figure
set(gcf, 'Color', 'w');
plot(t, pout, 'b', 'Linewidth', 2);
hold on
plot(gt(1, :), gt(2, :), 'r--', 'Linewidth', 2);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$p_{\mathrm{out}}$ [Pa]','Fontsize',16,'interpreter','latex');
legend('WDF','SSC','Fontsize',16,'interpreter','latex');
title('Output Pressure - Time Domain','Fontsize',18,'interpreter','latex');

% Frequency domain Plots
nfft = 2^20;
res = fs/nfft;
f = (0:nfft/2-1) * res;

ir = impzest(vin/A, pout');
ir_gt = impzest(vin/A, gt(2, 1:end-1)');

tf = fft(ir, nfft);
tf_gt = fft(ir_gt, nfft);

abs_tf = abs(tf(1:nfft/2));
abs_tf_gt = abs(tf_gt(1:nfft/2));

figure
set(gcf, 'Color', 'w');
semilogx(f, 20*log10(abs_tf/2e-5), 'b', 'Linewidth', 2);
hold on
semilogx(f, 20*log10(abs_tf_gt/2e-5), 'r--', 'Linewidth', 2);
grid on;
xlim([500, 20000])
xlabel('Frequency [Hz]','Fontsize',16,'interpreter','latex');
ylabel('$\mathrm{SPL}\,[\mathrm{dB}_\mathrm{SPL}]$','Fontsize',16,'interpreter','latex');
legend('WDF','SSC','Fontsize',16,'interpreter','latex');
title('Output Sound Pressure Level - Frequency Domain','Fontsize',16,'interpreter','latex');

%% Error Plots

figure
set(gcf, 'Color', 'w');
hold on;
plot(t, pout - gt(2, 1:end-1), 'k', 'Linewidth', 2);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$\mathcal{E}_{\mathrm{out}}$ [Pa]','Fontsize',16,'interpreter','latex');
title('Error Signal','Fontsize',16,'interpreter','latex');

%% Compute Mean Squared Error (MSE)

mse = mean((pout - gt(2, 1:end-1)).^2);
disp('MSE = ')
disp(mse)
