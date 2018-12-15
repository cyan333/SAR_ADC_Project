%% Initialization for setting
clear
clc
close all
format long

%% Initialization for parameters
Avth_pmos=2.5e-3; %V.um
Abita_pmos=0.5;  % %um
Avth_nmos=3.3e-3 % v.um
Abita_nmos=0.6; % %.um

td = 2^16-1;
%clock frequency
fclk =  2*pi*6e6;
Tclk = 1/fclk;
bw = 20e3*1;    % Bandwidth
OSR = fclk/(2*bw);

%sampling frequency
Fs= 2*pi*6e6/13;
Ts=1/Fs;

fin = 2*pi*200000;
% fin=Fs./2;

%also defined in DAC
W = 0.13; %comp input
L=0.065;
Cp = W*L*6e-15;
C_input_cancel = 1000*Cp;
input_cancel_gain = C_input_cancel./(C_input_cancel + Cp);
Vos = Avth_pmos*sqrt(W*L);

Ru = 200;
Cu = 5.75e-14;
sigma_C = 0.001322917989585;
sigma_R = 0.001841423909340;
% R_mismatch = normrnd(Ru,Ru*sigma_R,1,64);
% C_mismatch = normrnd(5.75e-14,Cu*sigma_C,1,32);
R_mismatch = normrnd(200,0,1,64);
C_mismatch = normrnd(5.75e-14,0,1,32);

sim('SAR_ADC_ideal')

%% FFT plotting and SNR calculation
z=sar_adc_out;
L = length(z);             % Length of signal
z=blackmanharris(L).*z;
figure(1)
plot(z);

Y=fft(z);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*2*P1(2:end-1);
figure(2)
f = Fs*(0:(L/2))/L;
plot(f,P1) 

figure(3)
fft_out = 20.*log10(P1);
SNR = snr(z,Fs)
 plot(f./Fs,fft_out); 
% t=1:1:N/2;
% 
% win_z=blackmanharris(N).*z;
% 
% Yfft=fft(win_z);
% % y_abs = 2.*abs(Yfft(1:N/2));
% 
% P2 = abs(Yfft/N);
% P1 = P2(1:N/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% % P1 = y_abs(1:L/2+1);
% % P1(2:end-1) = 2*P1(2:end-1);
% % f = Fs./N.*(t);
% f = (0:(N/2))/N;
% 
% fft_out = 20.*log10(P1);
% 
% %% Generate SNR
% SNR = snr(win_z,Fs)
% 
% figure(1)
% plot(f,fft_out); 
% title('Single-Sided Amplitude Spectrum of SAR Output')
% xlabel('Fin/Fs')
% ylabel('Normalized Amplitute (dBFS)')
% 
% figure(2)
% snr(z,Fs);
% 
% R_mismatch = normrnd(Ru,Ru*sigma_R,1,64);
% C_mismatch = normrnd(5.75e-14,Cu*sigma_C,1,32);
% 
% sim('SAR_ADC_ideal')
% z=sar_adc_out;
% THD = thd(z,Fs,6)
% SNDR = sinad(z,Fs)
% DR = sfdr(z,Fs)
% 
% figure(3)
% thd(z,Fs,6,'aliased');
% 
% figure(4)
% sinad(z,Fs);
% 
% figure(5)
% sfdr(z,Fs);




