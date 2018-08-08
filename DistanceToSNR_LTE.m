function SNR = DistanceToSNR_LTE(d, vFD, nScPerSlot)
% function SNR = DistanceToSNR(x, y, vFD, nScPerSlot)
% function SNR = DistanceToSNR(d, angle, vFD, nScPerSlot)
% Channel Model for LTE
% 
% inputs:
%   d: link distance in km
%   vFD: vector of small scale fading gain from database
% outputs:
%   SNR: sub-channel SNR vector

%% System parameters
OperationFrequency = 2000;          % 2500 MHz
hb = 30;                            % BS antenna height 32 m
hr = 1.5;                           % MS antenna height 1.5 m
Nf = 9;                             % MS noise figure 7 dB
Psc = 43-10*log10(12*100);             % power per sub-carrier dBm
%Psc = 100;                         % power per sub-carrier dBm
Nth = -174;                         % thermal noise dBm/Hz
SDsigma = 8;                        % log-normal shadowing SD sigam 8 dB
BWslot = nScPerSlot*15*1000;     % slot bandwidth in Hz
% nSc = 48;                           % number of sub-carriers per sub-channel
% nCh = 16;                           % number of sub-channels
% d = sqrt(x*x+y*y);
% angle = 0;                           % antenna angle
% antennaAnglePower = 0;              % antenna angle power

%% Calculation of antenna angle
% sinValue = y/d;
% 
% if sinValue >= 0
%     if sinValue <= (sqrt(3.0)/2.0)
%         if x > 0
%             angle = 2*90*asin(sinValue)/pi;
%         else
%             angle = 60 - 2*90*asin(sinValue)/pi;
%         end
%     elseif sinValue > (sqrt(3.0)/2.0)
%         if x > 0
%             angle = 30 + (90 - 2*90*asin(sinValue)/pi);
%         else
%             angle = 2*90*asin(sinValue)/pi - 60;
%         end
%     end
% elseif sinValue < 0
%     if sinValue >= (sqrt(3.0)/2.0)*(-1.0)
%         if x > 0
%             angle = (-1)*2*90*asin(sinValue)/pi;
%         else
%             angle = 60 - (-1)*2*90*asin(sinValue)/pi;
%         end
%     elseif sinValue < (sqrt(3.0)/2.0)*(-1.0)
%         if x > 0
%             angle = 30 + (90 - (-1)*2*90*asin(sinValue)/pi);
%         else
%             angle = (-1)*2*90*asin(sinValue)/pi - 60;
%         end
%     end  
% end
angle = 120*rand - 60;

%% COST 231 Hata Model (suburban)
% f = OperationFrequency;             % carrier frequency in MHz
% ahm = (1.1*log10(f)-0.7)*hr - (1.56*log10(f)-0.8);      % MS antenna height correction factor
% C = 0;                              % 0 dB:medium cities and suburban areas, 3 dB: metropolitan areas
% L = 46.3 + 33.9*log10(f) - 13.82*log10(hb) - ahm + (44.9-6.55*log10(hb))*log10(d) + C;      % path loss

f = OperationFrequency;             % carrier frequency in MHz
ahm = (1.1*log10(f)-0.7)*hr - (1.56*log10(f)-0.8);      % MS antenna height correction factor
C = 0;                              % 0 dB:medium cities and suburban areas, 3 dB: metropolitan areas
L0 = 46.3 + 33.9*log10(f) - 13.82*log10(hb) - ahm + C;      % path loss - (44.9-6.55*log10(hb))*log10(d)
%L1 = L0 + (44.9-6.55*log10(hb))*log10(d);
L1 = 128.1 + 37.6*log10(d);

%% shadowing
SD = SDsigma*randn;

%% small scale fading
FD = vFD;

%% mobile station SNR
% Pt = Psc + 10*log10(nScPerSlot);               % transmission
% % Gt = 15 - 9*rand;                       % transmitter antenna gain
% angle = 120*rand-60;
% antennaAnglePower = (-1.0)*min(12.0*(angle/70.0)*(angle/70.0), 20);
% Gt = 15 + antennaAnglePower;
% Gr = -1;                                % receiver antenna gain
% Pr = Pt + Gt + Gr - L - SD + FD;        % power
% N = Nth + 10*log10(BWslot) + Nf;          % noise
% SNR = Pr - N;
% % SNR = SNRtemp(1:nCh);

Pt = Psc + 10*log10(nScPerSlot);               % transmission
% Gt = 15 - 9*rand;                       % transmitter antenna gain
% antennaAnglePower = (-1.0)*min(12.0*(angle/70.0)*(angle/70.0), 20);
% Gt1a = 15 + antennaAnglePower;
Gr = -1;                                % receiver antenna gain
% Pr = Pt + Gt + Gr - L - SD + FD;        % power formula
Pr0 = Pt + Gr + FD;                  % power.....base


Pr1 = Pr0 - L1 - SD(1);                   % power.....for the base station


% affecton of all the three sectors
% Gt = 15 - 9*rand;                       % transmitter antenna gain
antennaAnglePower = (-1.0)*min(12.0*(angle/70.0)*(angle/70.0), 20);
Gt1a = 15 + antennaAnglePower;
Pr1a = Pr1 + Gt1a;        % power, sector a
antennaAnglePower = (-1.0)*min(12.0*((120-angle)/70.0)*((120-angle)/70.0), 20);
Gt1b = 15 + antennaAnglePower;
Pr1b = Pr1 + Gt1b;        % power, sector b
Gt1c = 15 - 20;
Pr1c = Pr1 + Gt1c;        % power, sector c

N = Nth + 10*log10(BWslot) + Nf;          % noise
SNR = 10*log10(power(10, Pr1a/10)+power(10, Pr1b/10)+power(10, Pr1c/10)) - N;
