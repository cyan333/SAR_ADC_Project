clear;


% Fs = 1000;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 2001;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% t1 = max(t);

Fs= 2*pi*6e6/13;
T = 1/Fs;             % Sampling period       
L = 29000;             % Length of signal
t = (0:L-1)*T;        % Time vector
t1 = max(t);

% S = sin(2*pi*120*t);

% X = S + 2*randn(size(t));
% X = S;
load('SAR_output_z.mat','z');
figure(1)
plot(z);


figure(3)
plot(100000*t(1:50),z(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')

Y=fft(z);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*2*P1(2:end-1);
figure(2)
f = Fs*(0:(L/2))/L;
plot(f,P1) 

%% Initialization for parameters
% Avth_pmos=2.5e-3; %V.um
% Abita_pmos=0.5;  % %um
% Avth_nmos=3.3e-3 % v.um
% Abita_nmos=0.6; % %.um
% 
% td = 2^16-1;
% %clock frequency
% fclk =  2*pi*6e6;
% Tclk = 1/fclk;
% bw = 20e3*1;    % Bandwidth
% OSR = fclk/(2*bw);
% 
% %sampling frequency
% Fs= 2*pi*6e6/13;
% Ts=1/Fs;
% 
% % fin = 2*pi*200000;
% fin=Fs./2;
% 
% %also defined in DAC
% W = 0.13; %comp input
% L=0.065;
% Cp = W*L*6e-15;
% C_input_cancel = 100000*Cp;
% input_cancel_gain = C_input_cancel./(C_input_cancel + Cp);
% Vos = Avth_pmos*sqrt(W*L);
% 
% Ru = 200;
% Cu = 5.75e-14;
% sigma_C = 0.001322917989585;
% sigma_R = 0.001841423909340;
% % R_mismatch = normrnd(Ru,Ru*sigma_R,1,64);
% % C_mismatch = normrnd(5.75e-14,Cu*sigma_C,1,32);
% R_mismatch = normrnd(200,0,1,64);
% C_mismatch = normrnd(5.75e-14,0,1,32);
% 
% sim('SAR_ADC_ideal')
% z=sar_adc_out;
% L = length(z);   
% 
% figure(2)
% plot(z)
% 
% figure(1)
% Y=fft(z);
% 
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% f = Fs*(0:(L/2))/L;
% plot(f,P1) 
% 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')


