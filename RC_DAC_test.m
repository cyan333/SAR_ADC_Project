clear;

Vin=[1 2 3];
D1=1;
D2=0;
D3=0;
D4=0;
D5=0;

S5=0;
S4=0;
S3=0;
S2=0;
S1=0;
S0=0;

sigma_C = 0.001322917989585;
sigma_R = 0.001841423909340;

N_CDAC = 5;
N_RDAC = 6;
% N_DAC = N_CDAC+N_RDAC;

% Cu = 57.596e-15;
% Ru = 200;

Vref = 1;

%% Switch for C DAC
%  D5 D4 D3 D2 D1
D = [D5 D4 D3 D2 D1]; 

%% Switch for R DAC
thisSwitch = S5*2^5 + S4*2^4 + S3*2^3 + S2*2^2 + S1*2^1 + S0*2^0;


%% Generate C, R with mismatch 32 number of C and 64 number of R
C = normrnd(1,sigma_C,[1,(2^5)]);
% C = normrnd(1,sigma_C,[1,(2^N_CDAC)]);

R = normrnd(1,sigma_R,[1,(2^N_RDAC)]);
% R = normrnd(1,sigma_R,[1,(2^N_RDAC)]);


C_binary = C(1);
N_b = 2;

%% generate binary Cap value for cap DAC
for i = 1:N_CDAC-1
    C_binary = [sum(C(N_b:(N_b+N_b-1))), C_binary];
    N_b = N_b*2;
end
C0 = C(2^N_CDAC);
Ctotal = sum(C);

%% Generate Thermometer R DAC

R_seq = zeros(1,length(R));
for i=2:length(R)
   R_seq(:,i) = [sum(R(1:i-1))];

end

% R_seq = fliplr(R_seq);
% R_seq = repmat(R_seq(1:2^N_RDAC),2^N_RDAC,1);

Rtotal = sum(R);

%% Generate Cap voltage
CD = zeros(1,5);
for i = 1:5
   CD(:,i) = C_binary(i)*D(i);
end

Vout_cap = (sum(CD)*Vref) / Ctotal;

%% Generate Resistor Voltage
Vout_res = (Vref * C0 * R_seq(thisSwitch+1)) / (Rtotal * Ctotal);

Vout = Vout_cap + Vout_res + Vin;