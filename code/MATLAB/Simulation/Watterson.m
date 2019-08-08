%% Header 
% Run our transmitted file through a filter which provides the Watterson
% model
% https://www.mathworks.com/examples/matlab-communications/mw/comm_product-chandemo_hf-hf-ionospheric-channel-models#3


%% Initializations
% M = 4;                        % Modulation order
% qpskMod = comm.QPSKModulator(0); % 4-PSK modulator object with 0 phase offset

% Rsym = 1200;                  % Input symbol rate
% Rbit = Rsym * log2(M);        % Input bit rate
% Nos = 4;                      % Oversampling factor
% ts = (1/Rbit) / Nos;          % Input sample period

sample_rate = 200e3;            % sampling rate of 200kHz
ts = 1/sample_rate;             % Input sample period
txtime = 0.5;                   % in s. Depends on the TxFile

%% Create doppler for first magneto-ionic component
fd = 10; % Chosen maximum Doppler shift for simulation
sGauss1 = 2.0;
fGauss1 = -5.0;
dop1 = doppler('BiGaussian', ...
               'NormalizedStandardDeviations', [sGauss1/fd 1/sqrt(2)], ...
               'NormalizedCenterFrequencies',  [fGauss1/fd 0], ...
               'PowerGains',                   [0.5        0])
           
%% Simulate first magneto-ionic component (single-path Rayleigh channel)
chan1 = comm.RayleighChannel( ...
    'SampleRate',          1/ts, ...
    'MaximumDopplerShift', fd, ...
    'DopplerSpectrum',     dop1, ...
    'RandomStream',        'mt19937ar with seed', ...
    'Seed',                99, ...
    'PathGainsOutputPort', true);

%% Create doppler for second magneto-ionic component
sGauss2 = 1.0;
fGauss2 = 4.0;
dop2 = doppler('BiGaussian', ...
               'NormalizedStandardDeviations', [sGauss2/fd 1/sqrt(2)], ...
               'NormalizedCenterFrequencies',  [fGauss2/fd 0], ...
               'PowerGains',                   [0.5        0])
           
%% Simulate second magneto-ionic component (single-path Rayleigh channel) 
chan2 = comm.RayleighChannel( ...
    'SampleRate',          1/ts, ...
    'MaximumDopplerShift', fd, ...
    'DopplerSpectrum',     dop2, ...
    'RandomStream',        'mt19937ar with seed', ...
    'Seed',                999, ...
    'PathGainsOutputPort', true);

%% Channels complex gains are found by summing
gGauss1 = 1.2;
gGauss2 = 0.25;

% This needs attention
% Nsamp = 2e6;                % Total number of channel samples
% Nsamp_f = 1000;             % Number of samples per frame
% Nframes = Nsamp/Nsamp_f;    % Number of frames

Nsamp = sample_rate*txtime + 1; % Total number of channel samples
Nsamp_f = Nsamp;             % Number of samples per frame
Nframes = Nsamp/Nsamp_f;    % Number of frames

s = zeros(Nsamp, 1);  y = zeros(Nsamp, 1);
for iFrames = 1:Nframes
%     inputSig = qpskMod(randi([0 M-1], Nsamp_f, 1));
    [inputSig, fs] = audioread('../Transmission-files/100kHz_0.5s_fs200k.wav');
    inputSig = inputSig(:,1);
    [s1, y1] = chan1(inputSig);
    [s2, y2] = chan2(inputSig);
    s((1:Nsamp_f)+(iFrames-1)*Nsamp_f) = sqrt(gGauss1) * s1 ...
                                       + sqrt(gGauss2) * s2;
    y((1:Nsamp_f)+(iFrames-1)*Nsamp_f) = sqrt(gGauss1) * y1 ...
                                       + sqrt(gGauss2) * y2;
end

%% Doppler spectrum is estimated from the complex path gains and plotted
hFig = figure;
pwelch(y, hamming(Nsamp/100), [], [], 1/ts, 'centered');
axis([-0.1 0.1 -80 0]);
legend('Simulation');

%% Theoretical bi-Gaussian Doppler spectrum is overlaid
f = -1/(2*ts): 0.1 :1/(2*ts);
Sd = gGauss1 * 1/sqrt(2*pi*sGauss1^2) * exp(-(f-fGauss1).^2/(2*sGauss1^2)) ...
   + gGauss2 * 1/sqrt(2*pi*sGauss2^2) * exp(-(f-fGauss2).^2/(2*sGauss2^2));

hold on;
plot(f(Sd>0)/1e3, 10*log10(Sd(Sd>0)), 'k--');
legend('Simulation', 'Theory');