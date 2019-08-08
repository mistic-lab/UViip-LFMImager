% sweep5a

% manually find one cycle of sweep file

clear variables; close all;

fs=2e6; % sampling rate
T=0.25; %wavefile time in sec
%fc=256; % carrier frequency

% linear change in frequency;
fc=(fs/4)*(1:fs*T)/(fs*T); % 0 to fs/4
% fs/4 is the sweep bandwidth = bandwidth of signal
% want to write cos (2 pi fc t)
% but in sampled systems t = n*ts = n/fs

% reference sweep waveform
%twopi_fc_t=(1:fs*T)*2*pi*fc/fs; from amfmdem3
twopi_fc_t=fc.*(1:fs*T)*2*pi/fs; 

% equation for reference sweep waveform
c_t=0.5*cos(twopi_fc_t)';
c2_t = [c_t;c_t];

N=length(c_t);
c_t=c_t(1:N)';

% figure(1)
% plot(c_t(1:N));

% multipath
delay=3; % #samples plus 1
h_t=[1 zeros(1,delay) 1 zeros(1,fs*T-delay-1)];

m_t=conv(c_t,h_t);


r1=xcorr(c_t,[c_t c_t]);
r2=xcorr(m_t,[c_t c_t]);
mc2=20*log10(abs(r2/max(abs(r2))));

figure(2)
plot(mc2);
% axis([0 length(r2) -100 0]);
title('mc2')
xlabel('Samples')
ylabel('Magnitude of correlation [dB]')

figure(3)
plot(fc)
ylabel('Samples')
xlabel('Frequency')
title('fc')



% time = sample number/sampling rate
% distance = velocity * time = 340 m/sec * sec
% find max value of reference value
% refsample = find(max(cc))