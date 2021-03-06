%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2018/12/3
% Purpose:  This file is to module R C DAC
%          
%   Copyright (c) 2018 by Shanshan Xie
%   for SAR ADC project in ADC course
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

sigma_deltaC_C = 3.307294973963294e-04;
sigma_deltaR_R = 4.639952116535146e-04;

N_CDAC = 5;
N_RDAC = 6;
N_DAC = N_CDAC+N_RDAC;
Cu = 230e-15;
Ru = 180;

Vref = 1;
Vin_neg = 1;

%% Switch for C DAC
%  D5 D4 D3 D2 D1
% D = [0 0 0 0 1]; 
D = dec2bin(0:(2^N_CDAC-1)) - '0';
%% Switch for R DAC
% thisSwitch = 1; %swich from 0 to 3 Vref0  = 0, Vref63 = 63/64Vref

S = [zeros(1,63),1];

for i=2:64
    S(i,:) = [S(i-1, 2:64),0];
end

%% Generate C, R with mismatch 32 number of C and 64 number of R
C = normrnd(1,sigma_deltaC_C,[1,(2^N_CDAC)]);
% C = normrnd(1,sigma_deltaC_C,[1,(2^N_CDAC)]);

R = normrnd(1,sigma_deltaR_R,[1,(2^N_RDAC)]);
% R = normrnd(1,sigma_deltaR_R,[1,(2^N_RDAC)]);

C_binary = C(1);
N_b = 2;

%% generate binary Cap value for cap DAC
for i = 1:N_CDAC-1
    C_binary = [sum(C(N_b:(N_b+N_b-1))), C_binary];
    N_b = N_b*2;
    i=i+1;
end
C0 = C(2^N_CDAC);
Ctotal = sum(C);

C_array = repmat(C_binary,2^N_DAC,1);

%% Generate Thermometer R DAC
i = 1;
while i<length(R)
   R_seq(:,length(R)-i) = [sum(R(1:i))];
   i = i+1;
end
R_seq = [R_seq, 0];
R_seq = repmat(R_seq,2^N_DAC,1);

Rtotal = sum(R);

S = repmat(S,2^N_CDAC,1);

SR = S.*R_seq;

%% Generate Cap voltage

D = kron(D,ones(2^N_RDAC,1));

CD = C_array.*D;

Vx_cap = (sum(CD,2)*Vref) ./ Ctotal;

%% Generate Resistor Voltage


Vx_res = (Vref .* C0 .* sum(SR,2)) ./ (Rtotal .* Ctotal);
% 
% Vx = Vx_cap + Vx_res + Vin_neg;

Vx = Vx_cap + Vx_res;


step_avg = Vref/2^N_DAC;

%% Calculate DNL

for i=2:2^N_DAC
   step(i-1) = Vx(i) - Vx(i-1);
end

DNL = (step-step_avg)./step_avg;

%% Calculate INL
vout_ideal = 0;

for i=2:2^N_DAC
   vout_ideal(i,1) = (i-1)*step_avg;
end

INL = (Vx-vout_ideal)./step_avg;


%% Plot

figure(1)
subplot(2,1,1);
plot(DNL,'DisplayName','DNL','LineWidth',2);
ylabel('DNL[LSB]','FontSize',12,'FontWeight','bold');
xlabel('Output Code','FontSize',12,'FontWeight','bold');
grid on
legend('show');
xlim([0,2^N_DAC]);

subplot(2,1,2);
plot(INL,'DisplayName','INL','LineWidth',2);
ylabel('INL[LSB]','FontSize',12,'FontWeight','bold');
xlabel('Output Code','FontSize',12,'FontWeight','bold');
grid on
legend('show');
xlim([0,2^N_DAC]);





