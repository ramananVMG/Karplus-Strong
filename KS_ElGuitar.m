% Karplus-Strong Electric Guitar
% Mac Porter
% February 12, 2019

% This code implements a basic electric guitar model using the Karplus-
% Strong algorithm. It is very similar to the acoustic guitar model but 
% without the guitar body IR. Feedback is added so that the output of the
% string is used as input to the delay line. Afterwards, soft clipping
% distortion is added.
% -------------------------------------------------------------------------
clear all;

% Parameters/initial values
% -------------------------------------------------------------------------
f0 = 440;               % Fundamental freq of string (Hz)
rho = 0.97;             % Loss factor
R = 0.9;                % Dynamics filter coefficient (strength of pluck)
mu = 0.25;              % Fraction of dist between bridge and pluck point
fb_gain = 0.02;         % Feedback gain (0 to 1)
tEnd = 15.0;            % Duration of output (s)

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
v = conv(noise,pluck);              % Convolution of noise with pluck
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
y = zeros(L,1);                 % Initialize output
x = [u;zeros(M,1)];
ylast = 0;
y2last = 0;
z = zeros(N,1);
for k = 1:N:L-N
    fb(k:k+N-1) = z;            % Feedback term
    for n = k:k+N-1
        i = n>N;
        j = n>N+1;
        ynew = x(n)+halfRho*(i*y(n-N*i)+j*y(n-(N+1)*j))+fb_gain*fb(n);% K-S
        y(n) = C*ynew+ylast-C*y2last;           % Allpass filter
        ylast = ynew;
        y2last = y(n);
    end
    z = y(k:k+N-1);             % Input = output
end

% Soft clipping distortion
y_dist = zeros(L,1);
y_dist_last = 0;
for n = 1:L
    x = y(n);
    if x>=1
        y_dist(n) = 2/3;
    elseif x>-1 && x<1
        y_dist(n) = x - (x.^3)/3;
    elseif x <=-1
        y_dist(n) = -2/3;
    end
    y_dist(n) = y_dist(n) + 0.9*y_dist_last;
    y_dist_last = y_dist(n);
end
y_dist = y_dist/max(abs(y_dist));     % Normalize
% -------------------------------------------------------------------------

% Audio Output
% -------------------------------------------------------------------------
soundsc(y_dist,Fs);
filename = sprintf('KS_ElGuitar_%dHz_%0.1fs.wav',f0,tEnd);
audiowrite(filename,y,Fs);
% -------------------------------------------------------------------------
