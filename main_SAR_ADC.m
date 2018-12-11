td = 2^16-1;
fclk = 6e6;
Tclk = 1/fclk;
bw = 20e3*1;    % Bandwidth
OSR = fclk/(2*bw);

Fs=6e6/13;
Ts=1/Fs;



sim('SAR_ADC_ideal')
z=sar_adc_out;
L = length(z);             % Length of signal
% M=2^1;
% NFFT = 2^(nextpow2(L/M)-1);% Next power of 2 from length of y
% win = ds_hann(NFFT)';
% win_rms=rms(win);
% win=win/win_rms; % if I do not normalize, to make peak 
% xw=z(1:NFFT).*win;
% Y = fft(xw);
% fb=ceil(NFFT/(OSR*2));
% y=abs(fft(xw/NFFT,NFFT));
% y(1:2)=y(1:2)*1e-0;
% si = find(max((y(1:fb+1)))==(y(1:fb+1)));
% SNDR_ramin=calculateSNR(y(1:fb),si-1);

% figure(3)
% plot(z)
Yfft=fft(z);
P2 = abs(Yfft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure(4)
plot(f*1e-3,P1) 

title('Single-Sided Amplitude Spectrum of SAR Output')
xlabel('f (kHz)')
ylabel('|P1(f)|')
snr(z,Fs)
thd(z,Fs,6)
sinad(z,Fs)
% figure(4)
% snr(z,Fs)
