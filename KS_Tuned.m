% Karplus-Strong with tuning correction
% Mac Porter
% February 6, 2019

% This code implements a basic plucked string model using the Karplus-
% Strong algorithm. The tuning is corrected using an allpass filter.
% The plot shows that the pitch produced is exactly in tune with
% f0 and its harmonics. This is obvious when compared to the version 
% without tuning correction, especially at high f0.
% -------------------------------------------------------------------------
clear all;

% Parameters/initial values
% -------------------------------------------------------------------------
f0 = 440;               % Fundamental freq of string (Hz)
rho = 0.98;             % Loss factor
R = 0.95;               % Dynamics filter coefficient (strength of pluck)
tEnd = 2.0;             % Duration of output (s)

Fs = 44100;             % Sample rate
Nexact = Fs/f0-0.5;     % Ideal delay line length (exact)
N = floor(Nexact);      % Delay line length (rounded)
P = Nexact-N;           % Fractional delay
C = (1-P)/(1+P);        % Allpass filter coefficient
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
x = [u;zeros(M,1)];                 % Input
y = zeros(L,1);                     % Initialize output
KS(1:N) = u;                        % First N samples of KS
KS(N+1) = halfRho*KS(1);            % N+1 sample of KS
y1(1) = C*KS(1);                    % 1st sample of allpass filter
y1(2:N+1) = C*KS(2:N+1)+KS(1:N);    % N+1 samples of allpass w/o last term
y(2:N+1) = y1(2:N+1)-C*y1(1:N);     % N+1 samples of allpass filter
ylast = KS(N+1);                    % Initialize ylast                          
for n = N+2:L
    ynew = halfRho*(y(n-N)+y(n-(N+1)));     % K-S starting at N+2
    y(n) = C*ynew+ylast-C*y(n-1);           % Allpass filter
    ylast = ynew;                           % Update ylast
end

y = y/max(abs(y));     % Normalize
% -------------------------------------------------------------------------

% Audio Output
% -------------------------------------------------------------------------
soundsc(y,Fs);
filename = sprintf('KS_TunedString_%dHz_%0.1fs.wav',f0,tEnd);
audiowrite(filename,y,Fs);
% -------------------------------------------------------------------------

% Plotting
% -------------------------------------------------------------------------
Y = abs(fft(y));                    % FFT magnitude of ouput
Y = Y/max(abs(Y));                  % Normalized FFT
f_exact = (f0:f0:Fs);               % Exact frequencies of harmonics
lines = ones(length(f_exact),1);    % Ones at exact frequencies for plot
f_axis = (0:L-1)*Fs/L;              % Frequency axis in Hz

subplot(2,1,2);
plot(f_axis,Y,'LineWidth',1);
hold on;
stem(f_exact,lines,'Marker','none','LineStyle','--','Color','k');
xlim([0 10*f0]);
legend('Output signal','Exact frequency of harmonics');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of output compared to exact frequencies of first 10 harmonics');

subplot(2,1,1);
plot(f_axis,Y,'LineWidth',1);
hold on;
stem(f_exact,lines,'Marker','none','LineStyle','--','Color','k');
xlim([f0-f0/20 f0+f0/20]);
legend('Output signal','Exact fundamental frequency');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of output compared to exact fundamental frequency');
% -------------------------------------------------------------------------
