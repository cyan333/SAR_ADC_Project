%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2018/12/1
% Purpose:  This file is to find sample frequency
%          
%   Copyright (c) 2018 by Shanshan Xie
%   for SAR ADC project in ADC course
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

fclk = 6e6; %Hz
N = 14;
fs = fclk/(1+N);