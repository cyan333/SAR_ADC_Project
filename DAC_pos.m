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

Cu = 230e-15;
Ru = 180;

Vref = 1;
Vin_pos = 1;

%% Switch for C DAC
%  D5 D4 D3 D2 D1
D = [0 0 0 0 1]; 

%% Switch for R DAC
thisSwitch = 1; %swich from 0 to 3 Vref0  = 0, Vref63 = 63/64Vref

S = [zeros(1,63),1];
for i=2:64
    S(i,:) = [S(i-1, 2:64),0];
end

%% Generate C, R with mismatch 32 number of C and 64 number of R
% C = normrnd(Cu,Cu*sigma_deltaC_C,[1,(2^N_CDAC)]);
C = normrnd(1,sigma_deltaC_C,[1,(2^N_CDAC)]);
% R = normrnd(Ru,Ru*sigma_deltaR_R,[1,(2^N_RDAC)]);
R = normrnd(1,sigma_deltaR_R,[1,(2^N_RDAC)]);


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

%% Generate Thermometer R DAC
i = 1;
while i<length(R)
   R_seq(:,i) = [sum(R(1:i))];
   i = i+1;
end
R_seq = [0, R_seq];
% R_seq = fliplr(R_seq);
% R_seq = repmat(R_seq(1:2^N_RDAC),2^N_RDAC,1);

Rtotal = sum(R);

%% Generate Cap voltage
for i = 1:5
   CD(:,i) = C_binary(i)*D(i);
   i = i+1;
end

Vx_cap = (sum(CD)*Vref) / Ctotal;

%% Generate Resistor Voltage
Vx_res = (Vref * C0 * R_seq(thisSwitch+1)) / (Rtotal * Ctotal);

Vx = Vx_cap + Vx_res + Vin_pos;












