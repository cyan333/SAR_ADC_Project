%% Initialization for setting
clear
clc
close all
format long

%% Initialization for parameters
td = 2^16-1;
fclk = 6e6;
Tclk = 1/fclk;
bw = 20e3*1;    % Bandwidth
OSR = fclk/(2*bw);

Fs=6e6/13;
Ts=1/Fs;

fin = 2*pi*2000;

sim('SAR_ADC_ideal')

%% FFT plotting and SNR calculation
z=sar_adc_out;
L = length(z);             % Length of signal

% figure(3)
% plot(z)
Yfft=fft(z);
P2 = abs(Yfft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

SNR = snr(z,Fs)
THD = thd(z,Fs,6)
SNDR = sinad(z,Fs)
DR = sfdr(z,Fs)
figure(1)
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of SAR Output')
xlabel('freq (kHz)')
ylabel('|P1(f)|')

figure(2)
snr(z,Fs);

figure(3)
thd(z,Fs,6,'aliased');

figure(4)
sinad(z,Fs);

figure(5)
sfdr(z,Fs);




