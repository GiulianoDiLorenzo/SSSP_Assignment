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
% Series adaptor (I)


% Parallel adaptor (H)


% Series adaptor (G)


% Parallel adaptor (F)


% Series adaptor (E)

% Series adaptor (D)


% Series adaptor (C)


% Parallel adaptor (B)


% Series adaptor (A)

%% Computing Scattering Matrices
% Series adaptor (I)
Z27 = R3;
Z26 = 2*L3/Ts;
Z25 = Z27 + Z26;
Ri = [Z25; Z26; Z27];
Si = eye(3) - 2*Ri*ones(1,3)/sum(Ri);

% Parallel adaptor (H)
Z24 = Z25;
Z23 = Ts/(2*C4);
Z22 = Z23*Z24/(Z23+Z24);
Gh = [1/Z22; 1/Z23; 1/Z24];
Sh = 2*Gh*ones(1,3)/sum(Gh) - eye(3);

% Series adaptor (G)
Z21 = Z22;
Z20 = 2*L2/Ts;
Z19 = Z21 + Z20;
Rg = [Z19; Z20; Z21];
Sg = eye(3) - 2*Rg*ones(1,3)/sum(Rg);

% Parallel adaptor (F)
Z18 = Z19;
Z17 = Ts/(2*C3);
Z16 = Z17*Z18/(Z17+Z18);
Gf = [1/Z16; 1/Z17; 1/Z18];
Sf = 2*Gf*ones(1,3)/sum(Gf) - eye(3);

% Series adaptor (E)
Z15 = Z16;
Z14 = Ts/(2*C2);
Z13 = Z15 + Z14;
Re = [Z13; Z14; Z15];
Se = eye(3) - 2*Re*ones(1,3)/sum(Re);

% Series adaptor (D)
Z12 = Z13;
Z11 = 2*L1/Ts;
Z10 = Z12 + Z11;
Rd = [Z10; Z11; Z12];
Sd = eye(3) - 2*Rd*ones(1,3)/sum(Rd);

% Series adaptor (C)
Z9 = Z10;
Z8 = R2;
Z7 = Z9 + Z8;
Rc = [Z7; Z8; Z9];
Sc = eye(3) - 2*Rc*ones(1,3)/sum(Rc);

% Parallel adaptor (B)
Z6 = Z7;
Z5 = Ts/(2*C1);
Z4 = Z5*Z6/(Z5+Z6);
Gb = [1/Z4; 1/Z5; 1/Z6];
Sb = 2*Gb*ones(1,3)/sum(Gb) - eye(3);

% Series adaptor (A)
Z3 = Z4;
Z2 = R1;
Z1 = Z2 + Z3;
Ra = [Z1; Z2; Z3];
Sa = eye(3) - 2*Ra*ones(1,3)/sum(Ra);

%% Initialization of Waves
a = zeros(27, 1);
b = zeros(27, 1);
s = zeros(27, 1);

% b(5) = 0;
% b(11) = 0;
% b(14) = 0;
% b(17) = 0;
% b(20) = 0;
% b(23) = 0;
% b(26) = 0;

%% Initialization of Output Signals

Fout = zeros(1, length(t));

%% Simulation Algorithm

for n = 1 : length(Fin)

    % Leaves update
    % a(2) = 0;
    % a(5) = b(5);
    % a(8) = 0;
    % a(11) = -b(11);
    % a(14) = b(14);
    % a(17) = b(17);
    % a(20) = -b(20);
    % a(23) = b(23);
    % a(26) = -b(26);
    % a(27) = 0;
    a = s;

    % Forward Scan
    b(25) = Si(1,:)*[a(25);a(26);a(27)];
    a(24) = b(25);
    b(22) = Sh(1,:)*[a(22);a(23);a(24)];
    a(21) = b(22);
    b(19) = Sg(1,:)*[a(19);a(20);a(21)];
    a(18) = b(19);
    b(16) = Sf(1,:)*[a(16);a(17);a(18)];
    a(15) = b(16);
    b(13) = Se(1,:)*[a(13);a(14);a(15)];
    a(12) = b(13);
    b(10) = Sd(1,:)*[a(10);a(11);a(12)];
    a(9) = b(10);
    b(7) = Sc(1,:)*[a(7);a(8);a(9)];
    a(6) = b(7);
    b(4) = Sb(1,:)*[a(4);a(5);a(6)];
    a(3) = b(4);
    b(1) = Sa(1,:)*[a(1);a(2);a(3)];

    % Local Root Scattering
    a_root = b(1);
    b_root = 2*Fin(n) - a_root;
    a(1) = b_root;

    % Backward Scan
    b(2) = Sa(2,:)*[a(1);a(2);a(3)];
    b(3) = Sa(3,:)*[a(1);a(2);a(3)];
    a(4) = b(3);
    b(5) = Sb(2,:)*[a(4);a(5);a(6)];
    b(6) = Sb(3,:)*[a(4);a(5);a(6)];
    a(7) = b(6);
    b(8) = Sc(2,:)*[a(7);a(8);a(9)];
    b(9) = Sc(3,:)*[a(7);a(8);a(9)];
    a(10) = b(9);
    b(11) = Sd(2,:)*[a(10);a(11);a(12)];
    b(12) = Sd(3,:)*[a(10);a(11);a(12)];
    a(13) = b(12);
    b(14) = Se(2,:)*[a(13);a(14);a(15)];
    b(15) = Se(3,:)*[a(13);a(14);a(15)];
    a(16) = b(15);
    b(17) = Sf(2,:)*[a(16);a(17);a(18)];
    b(18) = Sf(3,:)*[a(16);a(17);a(18)];
    a(19) = b(18);
    b(20) = Sg(2,:)*[a(19);a(20);a(21)];
    b(21) = Sg(3,:)*[a(19);a(20);a(21)];
    a(22) = b(21);
    b(23) = Sh(2,:)*[a(22);a(23);a(24)];
    b(24) = Sh(3,:)*[a(22);a(23);a(24)];
    a(25) = b(24);
    b(26) = Si(2,:)*[a(25);a(26);a(27)];
    b(27) = Si(3,:)*[a(25);a(26);a(27)];

    % Read Output
    Fout(n) = (a(27) + b(27))/2;

    % Previous state save
    s = b .* [0;0;0;0;1;0;0;0;0;0;-1;0;0;1;0;0;1;0;0;-1;0;0;1;0;0;-1;0];
    % s = b .* [1;0;1;1;1;1;1;0;1;1;-1;1;1;1;1;1;1;1;1;-1;1;1;1;1;1;-1;0];
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
