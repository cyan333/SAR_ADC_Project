
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2018/12/3
% Purpose:  This file is to calculate needed Cu and Ru for DAC
%          
%   Copyright (c) 2018 by Shanshan Xie
%   for SAR ADC project in ADC course
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
Vref = 1;

N_CDAC = 6; %C-DAC resolution
N_DAC = 12;
N_RDAC = N_DAC-N_CDAC; %R-DAC resolution

LSB_C = Vref./2.^(N_CDAC);
LSB_N = Vref./2.^N_DAC;
LSB_R = Vref./2.^N_RDAC;

N_sample = 1;

Linearity = 11;

Cap_num = 2.^(N_CDAC-1);
Res_num = 2.^N_RDAC;

Ab_R = 0.0125;
Rsq = 726;

Ab_C = 0.012;
Kc = 1.4e-15;

% Iload = 90e-6; %static current max


fclk = 6e6; %MHz
fconv = fclk/(N_DAC+N_sample);

%% R-DAC Mismatch Calculation
%if INL dominate
sigma_R = (0.5.*2.*2.^(N_DAC-Linearity) .* LSB_N ) ./ (sqrt(2).* 3.*sqrt(2.^(N_RDAC-2)) .* LSB_R);

%if DNL dominate
% sigma_R = (2^(N_RDAC-N_DAC) * 2^(N_DAC-Linearity)) / (3 * sqrt(2^(N_RDAC+1) - 2));

Ru = 200; %ohm

% Ru = Vref./(Iload * 2.^N_RDAC);

Rtotal = Ru .* Res_num;
Iload = 1e6./(Rtotal);

L_W_R = Ru/Rsq;

WL_R = 1e-12.*(Ab_R^2 ./ (2.*sigma_R.^2)); % unit area, unit um^2

L_R = sqrt(L_W_R.*WL_R); % unit: um

W_R = sqrt(WL_R./L_W_R);  % unit: um

%total area unit um^2
Area_R = Res_num .* WL_R;


%% C-DAC Mismatch Calculation
% DNL Dominate
% sigma_C = (0.5 .* 2 .* 2.^(N_DAC-Linearity) .* LSB_N) ./ (sqrt(2) .* 3 .* sqrt(2.^N_CDAC - 1) .* LSB_C);
sigma_C = (0.5.*2.*2.^(N_DAC-Linearity) .* LSB_N) ./ (sqrt(2) .* 3 .* sqrt(2.^(N_CDAC-1) - 1) .* LSB_C);

% min cap value to satisfy 99.73% yeild
Cu = 9 .* Ab_C^2 .* Kc .* (2.^(N_CDAC-1)-1) .* 2.^(2*N_DAC-2*N_CDAC);
Cu = Cu / 2^(2*N_DAC-2*Linearity);

% Cu = 500e-15; 

% WL_C = 1e-12*(Cu/Kc);  %cap unit size um^2
WL_C = 1e-12.*(Ab_C^2 ./ (2 .* sigma_C.^2)); % unit area, unit um^2


Area_C = (Cap_num .* WL_C);


%% Input Impedance

Ctotal = Cu .* Cap_num;

Rin = 2e-3./(fconv*Ctotal);

% figure(1)
% subplot(4,1,[1,2]);
% hold off
% plot(N_CDAC, Area_C.*1e12,'DisplayName','C AREA','LineWidth',2);
% 
% hold on
% plot(N_CDAC, Area_R.*1e12,'DisplayName','R AREA','LineWidth',2);
% 
% ylabel('Area','FontSize',12,'FontWeight','bold');
% xlabel('C DAC Bit','FontSize',12,'FontWeight','bold');
% grid on
% legend('show');
% xlim([0,12]);
% ylim([0,1600]);
% 
% subplot(4,1,3);
% hold off
% semilogy(N_CDAC, Iload,'DisplayName','Static Current','LineWidth',2);
% ylabel('Static Current','FontSize',12,'FontWeight','bold');
% xlabel('C DAC Bit','FontSize',12,'FontWeight','bold');
% grid on
% legend('show');
% xlim([0,12]);
% 
% subplot(4,1,4);
% hold off
% plot(N_CDAC, Rin,'DisplayName','Rin','LineWidth',2);
% xlim([0,12]);
% ylim([400, 2500]);
% 
% ylabel('R_{in}','FontSize',12,'FontWeight','bold');
% xlabel('C DAC Bit','FontSize',12,'FontWeight','bold');
% grid on
% legend('show');
% xlim([1,12]);


%% System INL
% CDAC_DNL = sqrt(2^(N_CDAC+1)-2) * sigma_C * (Vref/2^N_CDAC);
% CDAC_INL = sqrt(2^(N_CDAC-1)) * sigma_C * (Vref/2^N_CDAC);
% 
% RDAC_DNL = sqrt(2^(N_RDAC+1)-2) * sigma_R * (Vref/2^N_RDAC);
% RDAC_INL = sqrt(2^(N_RDAC-1)) * sigma_R * (Vref/2^N_RDAC);
% 
% sysINL = sqrt( CDAC_INL^2 + RDAC_INL^2 )
% sysINL_yield = sysINL * 3
% sysDNL = sqrt( CDAC_DNL^2 + RDAC_DNL^2 );
% sysDNL_yield = sysDNL * 3
% 
% LSB=1/2^N_DAC
% 
% Max_DNL_INL = 2/2^(N_DAC)
