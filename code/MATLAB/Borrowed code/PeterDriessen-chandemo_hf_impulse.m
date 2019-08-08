%% HF Ionospheric Channel Models
% This example shows how to simulate High-Frequency (HF) ionospheric
% channels, based on the models described in Recommendation ITU-R F.1487.
% In particular, it shows how to simulate the general Watterson channel
% model, and other simplified channel models used in the quantitative
% testing of HF modems. It makes use of the Rayleigh multipath fading
% channel object and the Gaussian and bi-Gaussian Doppler structures from
% Communications System Toolbox(TM).

% Copyright 2007-2016 The MathWorks, Inc.


%% ITU-R HF Channel Models: Overview
% In HF ionospheric radio communications, the transmitted signal can bounce
% off several times from the E and F layers of the ionosphere, which
% results in several propagation paths, also called modes [1]. Typically,
% the multipath delay spreads are large, as compared to mobile radio.
% Also, the signal can suffer from Doppler spread due to the turbulence of
% the ionosphere.  However, the fading rate is usually smaller than for
% mobile radio.
%
% Recommendation ITU-R F.1487 [1] proposes a general Gaussian scatter model
% for the simulation of HF ionospheric channels.  This model is based on
% Watterson's channel model [2].  Simpler models are also proposed in [1]
% for use in HF modem tests, with specified parameters. 


%% Initialization of Simulation-Specific Parameters 
% The simulation sampling rate is specified, and kept the same for the
% remainder of the example.  The input to the channel simulator is
% oversampled by a factor of four.

M = 4;                        % Modulation order
qpskMod = comm.QPSKModulator(0); % 4-PSK modulator object with 0 phase offset
%inputSig = qpskMod(randi([0 M-1], Nsamp_f, 1));

Rsym = 1200;                  % Input symbol rate
Rbit = Rsym * log2(M);        % Input bit rate
Nos = 4;                      % Oversampling factor
ts = (1/Rbit) / Nos;          % Input sample period


%% Watterson Channel Model
% The Watterson channel model consists of a tapped delay line, where each
% tap corresponds to a resolvable propagation path.  On each tap, two
% magneto-ionic components are present: each one is modeled as a complex
% Gaussian random process with a given gain and frequency shift, and whose 
% Doppler spectrum is Gaussian with a given standard deviation [2].
% Hence, each tap is characterized by a bi-Gaussian Doppler spectrum, which
% consists of two Gaussian functions in the frequency domain, each one with
% its own set of parameters (power gain, frequency shift, and standard
% deviation).

%%
% In this example, we follow the Watterson simulation model specified in
% [1], in which the complex fading process on each tap is obtained by
% adding two independent frequency-shifted complex Gaussian random
% processes (with Gaussian Doppler spectra) corresponding to the two
% magneto-ionic components.  This simulation model leads to a complex
% fading process whose envelope is in general _not_ Rayleigh distributed.
% Hence, to be faithful to the simulation model, we cannot simply generate
% a Rayleigh channel with a bi-Gaussian Doppler spectrum.  Instead, we
% generate two independent Rayleigh channels, each with a frequency-shifted
% Gaussian Doppler spectrum, gain-scale them, and add them together to
% obtain the Watterson channel model with a bi-Gaussian Doppler spectrum.
% For simplicity, we simulate a Watterson channel with only one tap.

%%
% A frequency-shifted Gaussian Doppler spectrum can be seen as a
% bi-Gaussian Doppler spectrum in which only one Gaussian function is
% present (the second one having a zero power gain).  Hence, to emulate the
% frequency-shifted Gaussian Doppler spectrum of each magneto-ionic
% component, we construct a bi-Gaussian Doppler structure such that one of
% the two Gaussian functions has the specified frequency shift and standard
% deviation, while the other has a zero power gain.

%%
% The first magneto-ionic component has a Gaussian Doppler spectrum with
% standard deviation |sGauss1|, frequency shift |fGauss1|, and power gain
% |gGauss1|.  A bi-Gaussian Doppler structure |dop1| is constructed such
% that the second Gaussian function has a zero power gain (its standard
% deviation and center frequency are hence irrelevant, and take on default
% values), while the first Gaussian function has a normalized standard
% deviation |sGauss1/fd| and a normalized frequency shift |fGauss1/fd|,
% where the normalization factor |fd| is the maximum Doppler shift of the
% corresponding channel.  In this example, since the gain of the second
% Gaussian function is zero, the value assigned to the gain of the first
% Gaussian function is irrelevant (we leave it to its default value of
% 0.5), because the associated channel System object created later
% normalizes the Doppler spectrum to have a total power of 1.
%
% Type |help doppler| for more information on how to construct a
% bi-Gaussian Doppler structure.

fd = 10; % Chosen maximum Doppler shift for simulation
sGauss1 = 2.0;
fGauss1 = -5.0;
dop1 = doppler('BiGaussian', ...
               'NormalizedStandardDeviations', [sGauss1/fd 1/sqrt(2)], ...
               'NormalizedCenterFrequencies',  [fGauss1/fd 0], ...
               'PowerGains',                   [0.5        0])

%%                     
% To simulate the first magneto-ionic component, we construct a single-path
% Rayleigh channel System object |h1| with a frequency-shifted Gaussian
% Doppler spectrum specified by the Doppler structure |dop1|. The average
% path power gain of the channel is 1 (0 dB).
   
chan1 = comm.RayleighChannel( ...
    'SampleRate',          1/ts, ...
    'MaximumDopplerShift', fd, ...
    'DopplerSpectrum',     dop1, ...
    'RandomStream',        'mt19937ar with seed', ...
    'Seed',                99, ...
    'PathGainsOutputPort', true);
           
%%
% Similarly, the second magneto-ionic component has a Gaussian Doppler
% spectrum with standard deviation |sGauss2|, frequency shift |fGauss2|,
% and power gain |gGauss2|.  A bi-Gaussian Doppler structure |dop2| is
% constructed such that the second Gaussian function has a zero power gain
% (its standard deviation and center frequency are hence irrelevant, and
% take on default values), while the first Gaussian function has a
% normalized standard deviation |sGauss2/fd| and a normalized frequency
% shift |fGauss2/fd| (again its power gain is irrelevant).

sGauss2 = 1.0;
fGauss2 = 4.0;
dop2 = doppler('BiGaussian', ...
               'NormalizedStandardDeviations', [sGauss2/fd 1/sqrt(2)], ...
               'NormalizedCenterFrequencies',  [fGauss2/fd 0], ...
               'PowerGains',                   [0.5        0])

%%                      
% To simulate the second magneto-ionic component, we construct a
% single-path Rayleigh channel System object |h2| with a frequency-shifted
% Gaussian Doppler spectrum specified by the Doppler structure |dop2|. 

chan2 = comm.RayleighChannel( ...
    'SampleRate',          1/ts, ...
    'MaximumDopplerShift', fd, ...
    'DopplerSpectrum',     dop2, ...
    'RandomStream',        'mt19937ar with seed', ...
    'Seed',                999, ...
    'PathGainsOutputPort', true);

%%
% We compute in the loop below the output to the Watterson channel in
% response to an input signal, and store it in |s|. In obtaining |s|, the
% first filter operation emulates the effect of the first magneto-ionic
% component, while the second filter operation emulates the effect of the
% second component.
%
% To obtain the desired power gains, |gGauss1| and |gGauss2|, of each
% magneto-ionic component,  we need to scale the output signal for each
% magneto-ionic component by their corresponding amplitude gains,
% |sqrt(gGauss1)| and |sqrt(gGauss2)|.
%
% Due to the low Doppler shifts found in HF environments and the fact that
% the bi-Gaussian Doppler spectrum is combined from two objects, obtaining
% measurements for the Doppler spectrum using the built-in visualization of
% the System objects is not appropriate. Instead, we store the channel's
% complex path gains and later compute the Doppler spectrum for each path
% at the command line.  In the loop below, the channel's complex path gains
% are obtained by summing (after scaling by the corresponding amplitude
% gains) the complex path gains associated with each magneto-ionic
% component, and then stored in |y|.

gGauss1 = 1.2;
gGauss2 = 0.25;

Nsamp = 2e6;                % Total number of channel samples
Nsamp_f = 1000;             % Number of samples per frame
Nframes = Nsamp/Nsamp_f;    % Number of frames

s = zeros(Nsamp, 1);  y = zeros(Nsamp, 1);
for iFrames = 1:Nframes
    inputSig = qpskMod(randi([0 M-1], Nsamp_f, 1));
    [s1, y1] = chan1(inputSig);
    [s2, y2] = chan2(inputSig);
    s((1:Nsamp_f)+(iFrames-1)*Nsamp_f) = sqrt(gGauss1) * s1 ...
                                       + sqrt(gGauss2) * s2;
    y((1:Nsamp_f)+(iFrames-1)*Nsamp_f) = sqrt(gGauss1) * y1 ...
                                       + sqrt(gGauss2) * y2;
end

%%
% The Doppler spectrum is estimated from the complex path gains and plotted.

hFig = figure;
pwelch(y, hamming(Nsamp/100), [], [], 1/ts, 'centered');
axis([-0.1 0.1 -80 0]);
legend('Simulation');

%%
% The theoretical bi-Gaussian Doppler spectrum is overlaid to the estimated
% Doppler spectrum.  We observe a good fit between both.

f = -1/(2*ts): 0.1 :1/(2*ts);
Sd = gGauss1 * 1/sqrt(2*pi*sGauss1^2) * exp(-(f-fGauss1).^2/(2*sGauss1^2)) ...
   + gGauss2 * 1/sqrt(2*pi*sGauss2^2) * exp(-(f-fGauss2).^2/(2*sGauss2^2));

hold on;
plot(f(Sd>0)/1e3, 10*log10(Sd(Sd>0)), 'k--');
legend('Simulation', 'Theory');


%% ITU-R F.1487 Low Latitudes, Moderate Conditions (LM) Channel Model
% Recommendation ITU-R F.1487 specifies simplified channel models used in
% the quantitative testing of HF modems.  These models consist of two
% independently fading paths with equal power.  On each path, the two
% magneto-ionic components are assumed to have zero frequency shift and
% equal variance: hence the bi-Gaussian Doppler spectrum on each tap
% reduces to a single Gaussian Doppler spectrum, and the envelope of the
% complex fading process is Rayleigh-distributed.
%
% Below, we construct a channel object according to the Low Latitudes,
% Moderate Conditions (LM) channel model specified in Annex 3 of ITU-R
% F.1487.  The path delays are 0 and 2 ms.  The frequency spread, defined
% as twice the standard deviation of the Gaussian Doppler spectrum, is 1.5
% Hz. The Gaussian Doppler spectrum structure is hence constructed with a
% normalized standard deviation of (1.5/2)/ |fd|, where |fd| is chosen as 1
% for simplicity (type |help doppler| for more information).

close(hFig);

fd = 1;
chan3 = comm.RayleighChannel( ...
    'SampleRate',          1/ts, ...
    'PathDelays',          [0 0.002], ...
    'AveragePathGains',    [0 0], ...
    'MaximumDopplerShift', fd, ...
    'DopplerSpectrum',     doppler('Gaussian', (1.5/2)/fd), ...
    'RandomStream',        'mt19937ar with seed', ...
    'Seed',                9999, ...
    'PathGainsOutputPort', true, ...
    'Visualization',       'Impulse response')

%%
% We have turned on the impulse response visualization in the Rayleigh
% channel System object. The code below simulates the LM channel and
% visualizes its bandlimited impulse response. By default, the channel
% responses for one of every four samples are visualized for faster
% simulation. In other words, for a frame of length 1000, the responses for
% the 1st, 5th, 9th, ..., 997th samples are shown. To observe the response
% for every sample, set the |SamplesToDisplay| property of |hLMChan| to
% |'100%'|.

Nsamp_f = 1e3;      % Number of samples per frame
Nframes = 100;       % Number of frames
for iFrames = 1:Nframes
   %inputSig = qpskMod(randi([0 M-1], Nsamp_f, 1));
   inputSig = [1; zeros(Nsamp_f-1,1)];
   outputSig(:,iFrames)=chan3(inputSig);
   %figure(iFrames)
   %plot(abs(outputSig(:,iFrames)))
   %axis([0 25 0 1.5]);
end
figure(1)
plot(abs(outputSig));
axis([0 25 0 1.5]);
figure(2)
plot(abs(outputSig(1,:)));
hold on
plot(abs(outputSig(20,:)),'r-');
legend('delay 0ms','delay 2ms');
hold off
figure(3)
plot(20*log10(abs(outputSig(1,:))));
hold on
plot(20*log10(abs(outputSig(20,:))),'r-');
legend('delay 0ms','delay 2ms');
hold off
figure(4)
plot(angle(outputSig(1,:)));
hold on
plot(angle(outputSig(20,:)),'r-');
legend('delay 0ms','delay 2ms');
hold off

figure(5)
plot(unwrap(angle(outputSig(1,:))));
hold on
plot(unwrap(angle(outputSig(20,:))),'r-');
legend('delay 0ms','delay 2ms');
hold off

%%
% We now turn on the Doppler spectrum visualization for the channel object
% to observe the theoretical and empirical Gaussian Doppler spectra for the
% first discrete path. Due to the very low Doppler shift, it may take a
% while to have the empirical spectrum converge to the theoretical
% spectrum.

% If there is no motion, then there is no Doppler shift, and
% all paths at a given delay add together to give
% some constant amplitude and phase.  
% In this case, signals spaced
% too close in time to be resolved wtihin the specified channel bandwidth
% are combined together in one channel filter tap at some constant delay,
% amplitude and phase.

% When there is motion, then 
% Signals traveling along different paths of different lengths
% can have different Doppler shifts, 
% corresponding to different rates of change in phase. Signals spaced
% too close in time to be resolved wtihin the specified channel bandwidth
% are combined together in one fading channel filter tap.
% The difference in Doppler shifts between different signal components 
% contributing to a single fading channel tap is known as the Doppler spread. 
% Channels with a large Doppler spread have signal components that are each 
% changing independently in phase over time.

% Because of the different incident angle, not only will the longer path 
% signal arrive later, but it will have a arrive with a different incident angle, 
% and therefore (due to the Doppler Effect) it will have a different frequency. 
% So the Doppler spread would be the difference of the two frequencies received 
% (even though there is only a single fixed frequency being transmitted).

release(chan3);
chan3.Visualization = 'Doppler spectrum';

Nsamp_f = 2e6;      % Number of samples per frame
Nframes = 80;       % Number of frames
for iFrames = 1:Nframes
   inputSig = qpskMod(randi([0 M-1], Nsamp_f, 1));
   chan3(inputSig);
end

%% ITU-R F.1487 High Latitudes, Disturbed Conditions (HD) Channel Model
% We use the function |stdchan| to construct a channel object according to
% the HD model.  For a list of standardized channel models supported by
% |stdchan|, type |help stdchan|. When using |stdchan| to construct ITU-R
% HF channel models, the maximum Doppler shift must be set to 1 Hz: this
% ensures that the Gaussian Doppler spectrum of the constructed channel has
% the correct standard deviation.

clear hLMChan;

fd = 1;
hdchan = stdchan(ts, fd, 'iturHFHD')

%%
% As before, we obtain the path gains by processing random data through the
% channel.  These path gains are stored in |y| for post-processing.

hdchan.StorePathGains = 1;
hdchan.ResetBeforeFiltering = 0;
hdchan.NormalizePathGains = 1;

Nsamp = 2e6;                    % Total number of channel samples
Nsamp_f = 1000;                 % Number of samples per frame
Nframes = Nsamp/Nsamp_f;        % Number of frames

s = zeros(Nsamp, 1);  y = zeros(Nsamp, 2);
for iFrames = 1:Nframes
    inputSig = qpskMod(randi([0 M-1], Nsamp_f, 1));
    s( (1:Nsamp_f) + (iFrames-1)*Nsamp_f ) = filter(hdchan, inputSig);
    y( (1:Nsamp_f) + (iFrames-1)*Nsamp_f, :) = hdchan.PathGains;
end

%%
% The Doppler spectrum for each path is estimated and plotted, alongside
% the theoretical Gaussian Doppler spectrum.

sigmaGaussian = 30/2;
f = -1/(2*ts) : 0.1 : 1/(2*ts);
Sd = 0.5 * 1/sqrt(2*pi*(sigmaGaussian*fd)^2) ...
         * exp(-f.^2/(2*(sigmaGaussian*fd)^2));

figure; 

hs1 = subplot(2, 1, 1); hold on;
pwelch(y(:,1), hamming(Nsamp/100), [], [], 1/ts, 'centered');
axis([-0.1 0.1 -80 0]);

plot(hs1, f(Sd>0)/1e3, 10*log10(Sd(Sd>0)), 'k--');
legend('Simulation', 'Theory');
title('Welch Power Spectral Density Estimate for Path 1');

hs2 = subplot(2, 1, 2); hold on;
pwelch(y(:,2), hamming(Nsamp/100), [], [], 1/ts, 'centered');
axis([-0.1 0.1 -80 0]);

plot(hs2, f(Sd>0)/1e3, 10*log10(Sd(Sd>0)), 'k--');
legend('Simulation', 'Theory');
title('Welch Power Spectral Density Estimate for Path 2');

        
%%
% References:
%
%  [1] Recommendation ITU-R F.1487, "Testing of HF modems with bandwidths
%      of up to about 12 kHz using ionospheric channel simulators," 2000.
%  [2] C. C. Watterson, J. R. Juroshek, and W. D. Bensema, "Experimental
%      confirmation of an HF channel model," IEEE(R) Trans. Commun. Technol.,
%      vol. COM-18, no. 6, Dec. 1970.

displayEndOfDemoMessage(mfilename)
