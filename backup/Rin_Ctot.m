%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2018/12/1
% Purpose:  This file is to find input resistance and R for R-DAC
%          
%   Copyright (c) 2018 by Shanshan Xie
%   for SAR ADC project in ADC course
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
Cu = 60e-15;  %fF

Vref = 1;

fclk = 6e6; %MHz

N_bit = 12;
N_sample = 1;

N_cap = 6;

fconv = fclk/(N_bit+N_sample);

Ctotal = Cu*2^(N_cap-1);

Rin = 1/(fconv*Ctotal)

N_res = N_bit-N_cap;

%Maximum static current from Vref is 100uA
I_static = 78.125e-6;
Ru = Vref/(I_static*2^N_res)

R = 200;
I = Vref/(R*2^N_res)




