% Karplus-Strong Acoustic Guitar
% Mac Porter
% February 11, 2019

% This code implements a basic acoustic guitar model using the tuning 
% corrected Karplus-Strong algorithm. The noise signal is convolved with
% a recorded pluck, which is then convolved with the IR of a Taylor guitar 
% body. The properties of convolution allow this to be done before being 
% input to the filter. The noise signal could be excluded, but I found that 
% it produces a more natural sound and allows for more variation. This
% aggregate signal is then pre-filtered with the same dynamics filter as 
% well as a comb filter which simulates pluck position. This position is 
% set by parameter 'mu'.
% -------------------------------------------------------------------------
clear all;

% Parameters/initial values
% -------------------------------------------------------------------------
f0 = 110;               % Fundamental freq of string (Hz)
rho = 0.98;             % Loss factor
R = 0.9;                % Dynamics filter coefficient (strength of pluck)
mu = 0.25;              % Fraction of dist between bridge and pluck point
tEnd = 5.0;             % Duration of output (s)

Fs = 44100;             % Sample rate
Nexact = Fs/f0-0.5;     % Ideal delay line length (exact)
N = floor(Nexact);      % Delay line length (rounded)
P = Nexact-N;           % Fractional delay
C = (1-P)/(1+P);        % Allpass filter coefficient
D = round(mu*N);        % For use in comb filter
halfRho = rho/2;
M = tEnd*Fs-N;             
L = M+N;                % Duration of output (samples)
% -------------------------------------------------------------------------

% Pre-filtering
% -------------------------------------------------------------------------
noise = 2*rand(N,1)-1;              % Random noise signal
pluck = audioread('pluck.wav');     % Recorded pluck
gtrIR = audioread('Guitar_IR.wav'); % Guitar body IR
v = conv(noise,pluck);              % Convolution of noise with pluck
v = conv(v,gtrIR);                  % Convolution with body IR
Lu = length(v);
w = zeros(Lu,1);                    % Initialize vectors
u = zeros(Lu,1);                    
wlast = 0;
for n = 1:Lu
    i = n > D;                      % i=1 if n > D, i=0 otherwise
    w(n) = (1-R)*v(n)+R*wlast;      % Dynamics filter
    u(n) = w(n)-i*w(n-D*i);         % Comb filter to simulate position
    wlast = w(n);
end
% -------------------------------------------------------------------------

% K-S Algorithm
% -------------------------------------------------------------------------
x = [u;zeros(L-Lu,1)];              % Input
y = zeros(L,1);                     % Initialize y
KS(1:N) = x(1:N);                   % First N samples of KS
KS(N+1) = halfRho*KS(1);            % N+1 sample of KS
y1(1) = C*KS(1);                    % 1st sample of allpass filter
y1(2:N+1) = C*KS(2:N+1)+KS(1:N);    % N+1 samples of allpass w/o last term
y(2:N+1) = y1(2:N+1)-C*y1(1:N);     % N+1 samples of allpass filter
ylast = KS(N+1);                    % Initialize ylast                          
for n = N+2:L
    ynew = x(n)+halfRho*(y(n-N)+y(n-(N+1)));    % K-S starting at N+2
    y(n) = C*ynew+ylast-C*y(n-1);               % Allpass filter
    ylast = ynew;                               % Update ylast
end

y = y/max(abs(y));      % Normalize
% -------------------------------------------------------------------------

% Audio Output
% -------------------------------------------------------------------------
soundsc(y,Fs);
filename = sprintf('KS_AcGuitar_%dHz_%0.1fs.wav',f0,tEnd);
audiowrite(filename,y,Fs);
% -------------------------------------------------------------------------
