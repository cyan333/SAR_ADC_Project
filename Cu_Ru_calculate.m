%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2018/12/3
% Purpose:  This file is to calculate needed Cu and Ru for DAC
%          
%   Copyright (c) 2018 by Shanshan Xie
%   for SAR ADC project in ADC course
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

N_RDAC = 6; %R-DAC resolution
N_CDAC = 5; %C-DAC resolution
N_DAC = 12;

N_sample = 1;

Linearity = 11;

Cap_num = 2^N_CDAC;
Res_num = 2^N_RDAC;

Ab_R = 0.0125;
Rsq = 726;

Ab_C = 0.012;
Kc = 1.4e-15;

Iload = 90e-6; %static current max
Vref = 1;

fclk = 6e6; %MHz
fconv = fclk/(N_DAC+N_sample);

%% R-DAC Mismatch Calculation
%if INL dominate
% sigma_deltaR_R = (2^((0.5*N_RDAC)+0.5) * 2^(N_DAC-Linearity)) / (3*2^N_DAC);

%if DNL dominate
% sigma_deltaR_R = (2^(N_RDAC-N_DAC) * 2^(N_DAC-Linearity)) / (3 * sqrt(2^(N_RDAC+1) - 2));
sigma_deltaR_R = (2^(N_RDAC-N_DAC)) / (3 * sqrt(2^(N_RDAC+1) - 2));
Ru = 180; %ohm

Iload = 1/(Ru*Res_num);

% Ru = Vref/(Iload * 2^N_RDAC);

Rtotal = Ru * Res_num;

L_W_R = Ru/Rsq;

WL_R = 1e-12*(Ab_R^2 / (2 * sigma_deltaR_R^2)); % unit area, unit um^2

L_R = sqrt(L_W_R*WL_R); % unit: um

W_R = sqrt(WL_R/L_W_R);  % unit: um

%total area unit um^2
Area_R = Res_num * WL_R


%% C-DAC Mismatch Calculation
% DNL Dominate
% sigma_deltaC_C = (2^(N_CDAC-N_DAC) * 2^(N_DAC-Linearity)) / (3 * sqrt(2^(N_CDAC+1) - 2))
sigma_deltaC_C = (2^(N_CDAC-N_DAC)) / (3 * sqrt(2^(N_CDAC+1) - 2))
% min cap value to satisfy 99.73% yeild
Cu = 9 * Ab_C^2 * Kc * (2^N_CDAC-1) * 2^(2*N_DAC-2*N_CDAC)
Cu = Cu / 2^(2*N_DAC-2*Linearity);

% Cu = 500e-15; 

% WL_C = 1e-12*(Cu/Kc);  %cap unit size um^2
WL_C = 1e-12*(Ab_C^2 / (2 * sigma_deltaC_C^2)); % unit area, unit um^2


Area_C = Cap_num * WL_C


%% Input Impedance

Ctotal = 50e-15 * Cap_num;

Rin = 1/(2*fconv*Ctotal);


%% System INL
CDAC_DNL = sqrt(2^(N_CDAC+1)-2) * sigma_deltaC_C * (Vref/2^N_CDAC);
CDAC_INL = sqrt(2^(N_CDAC-1)) * sigma_deltaC_C * (Vref/2^N_CDAC);

RDAC_DNL = sqrt(2^(N_RDAC+1)-2) * sigma_deltaR_R * (Vref/2^N_RDAC);
RDAC_INL = sqrt(2^(N_RDAC-1)) * sigma_deltaR_R * (Vref/2^N_RDAC);

sysINL = sqrt( CDAC_INL^2 + RDAC_INL^2 )
sysINL_yield = sysINL * 3
sysDNL = sqrt( CDAC_DNL^2 + RDAC_DNL^2 );
sysDNL_yield = sysDNL * 3

Max_DNL_INL = 2/2^(N_DAC)















