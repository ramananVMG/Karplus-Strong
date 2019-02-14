% Karplus-Strong Drum
% Mac Porter
% February 11, 2019

% This code implements a basic drum model using the Karplus-Strong drum
% algorithm. The blend factor 'b' controls the probability of the sign 
% within the K-S loop. For b = 0, the output is the same as the string
% model. For b = 0.5, the sound is like a snare drum. And for b = 1, the 
% output sounds like a 'plucked bottle' for high f0 and a harp for low f0.
% -------------------------------------------------------------------------
clear all;

% Parameters/initial values
% -------------------------------------------------------------------------
f0 = 60;                % Fundamental freq (Hz)
rho = 0.98;             % Loss factor
R = 0.2;                % Dynamics filter coefficient
b = 0.5;                % Blend factor
G = 10;                 % Gain factor for distortion
tEnd = 2.0;             % Duration of output (s)

Fs = 44100;             % Sample rate
N = round(Fs/f0-0.5);   % Delay line length
halfRho = rho/2;
M = tEnd*Fs-N;   
L = M+N;                % Duration of output (samples)
% -------------------------------------------------------------------------

% Pre-filtering
% -------------------------------------------------------------------------
v = 2*rand(N,1)-1;          % Random noise signal
u = zeros(N,1);             % Initialize pre-filtered signal
ulast = 0;
for n = 1:N
    u(n) = (1-R)*v(n)+R*ulast;     % Dynamics filter
    ulast = u(n);
end
% -------------------------------------------------------------------------

% K-S Algorithm
% -------------------------------------------------------------------------
y = [u;zeros(M,1)];             % Initialize y
y(N+1) = halfRho*y(1);          % N+1 sample
for n = 1:M-1
    p = rand <= b;      % p = 1 with probabilty b, and p = 0 otherwise
    % If p = 1, first term is removed. If p = 0, second term is removed.
    y(n+N+1) = halfRho*((1-p)*(y(n+1)+y(n)) - p*(y(n+1)+y(n))); % K-S
end

y = y/max(abs(y));      % Normalize
% -------------------------------------------------------------------------

% Audio Output
% -------------------------------------------------------------------------
soundsc(y,Fs);
filename = sprintf('KS_Drum_%0.1fb_%0.1fs.wav',b,tEnd);
audiowrite(filename,y,Fs);
% -------------------------------------------------------------------------