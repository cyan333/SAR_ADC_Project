%% Initialization for setting
clear
clc
close all
format long

%% DAC DNL & INL 1000 runs
for i=1:1000
   [INL(i,:), DNL(i,:)] = getDNLINL();
end

sigma_DNL = std(DNL);
sigma_INL = std(INL);

sigma_3_DNL = 3.*max(sigma_DNL);
sigma_3_INL = 3.*max(sigma_INL);

figure(1)
subplot(2,1,1);
plot(sigma_DNL,'DisplayName','DNL','LineWidth',2);
ylabel('\sigma_{DNL} [LSB]','FontSize',12,'FontWeight','bold');
xlabel('Output Code','FontSize',12,'FontWeight','bold');

grid on
legend('show');
xlim([0,2^11]);

subplot(2,1,2);
plot(sigma_INL,'DisplayName','INL','LineWidth',2);
ylabel('\sigma_{INL} [LSB]','FontSize',12,'FontWeight','bold');
xlabel('Output Code','FontSize',12,'FontWeight','bold');

grid on
legend('show');
xlim([0,2^11]);

%%%% Simulink Model %%%%
%% Initialization for parameters

%process parameter
Avth_pmos=2.5e-3;   % V.um
Abita_pmos=0.5;     % %um
Avth_nmos=3.3e-3;   % v.um
Abita_nmos=0.6;     % %.um

td = 0.001;         %sampling time
fclk =  2*pi*6e6;   %clock frequency
Tclk = 1/fclk;      %clock period
bw = 20e3*1;        %Bandwidth
OSR = fclk/(2*bw);

Fs= 2*pi*6e6/13;    %sampling frequency
Ts=1/Fs;            %sampling period

fin = 2*pi*20000;
% fin=Fs./2;

sigma_thermal = 0.1;

%Comparator Parameters
%This is also defined in DAC
W = 2;                      % transistor width  um
L  =1;                      % transistor length um
Cp = W*L*6e-15;             %input parasitic capacitor
C_input_cancel = 1000*Cp;   %input cancelation capacitor
Vos = Avth_pmos/sqrt(W*L)   % offset
input_cancel_gain = C_input_cancel./(C_input_cancel + Cp);

%DAC Parameters
Ru = 200;       % Unit Resistor
Cu = 5.75e-14;  % Unit Capacitor
sigma_C = 0.001322917989585;
sigma_R = 0.001841423909340;

% R_mismatch = normrnd(1,sigma_R,1,64);
% C_mismatch = normrnd(1,sigma_C,1,32);
R_mismatch = normrnd(Ru,Ru*sigma_R,1,64);
C_mismatch = normrnd(Cu,Cu*sigma_C,1,32);
% R_mismatch = normrnd(200,0,1,64);
% C_mismatch = normrnd(5.75e-14,0,1,32);

sim('SAR_ADC')

%% SNR calculation & plotting
z=sar_adc_out;

SNR = snr(z,Fs)

% win_z=blackmanharris(N).*z;
% SNR = snr(win_z,Fs)

figure(2)
snr(z,Fs);

%% SNDR 
% R_mismatch = normrnd(Ru,Ru*sigma_R,1,64);
% C_mismatch = normrnd(5.75e-14,Cu*sigma_C,1,32);
% sim('SAR_ADC')
% 
% %% SNR calculation & plotting
% z=sar_adc_out;

THD = thd(z,Fs,6)
SNDR = sinad(z,Fs)
DR = sfdr(z,Fs)
ENOB=(SNDR-1.761)/6.02

figure(3)
thd(z,Fs,6,'aliased');

figure(4)
sinad(z,Fs);

figure(5)
sfdr(z,Fs);




