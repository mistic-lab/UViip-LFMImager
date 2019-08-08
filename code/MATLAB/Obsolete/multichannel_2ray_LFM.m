%%Constants
fs = 1e6;                  %sampling frequency
T = 1;                      %length of sweep
t = 0:1/fs:T;                       %initial frequency
B = 5e5;                   %final frequency
%% Generate double sin sweep
sweep = exp((i*pi*B/T).*t.^2);    %generate sine sweep
sweep_x2 = [sweep sweep];    %put 2 of them back to back
%% Create two ray received signal with noise
delay = round(6.56e-4*fs);
delayed_sweep = [sweep_x2(delay+1:length(sweep_x2)) sweep_x2(1:delay)];
weight = 0.158;        %relative weight of weaker signal
received_signal = sweep_x2 + weight * delayed_sweep;
noisy_signal = awgn(received_signal, 10);    %add noise to the received signal
%% Time reversed conjugate signal
transmitted_conjuage_rev = exp((-i*pi*B/T).*(T-t).^2);
%% Matched Filter
pulse_comp_freq = fft(noisy_signal).*fft([transmitted_conjuage_rev zeros(1, length(sweep))]);
pulse_comp_time = ifft(pulse_comp_freq);
% % 
plot([t t], 20*log((abs(pulse_comp_time))/max(abs(pulse_comp_time))));
%% 



% spectrogram(received_signal, 1024, 840, 512, fs);
% plot(abs(fft(sweep_x2)));


% f2 = figure;
% f2.Position = [650, 500, 550, 600];
% spectrogram(noisy_signal, 1024, 840, 512, fs);
% title('2-ray multichannel simulation received signal');
% 
% %correlate received signal with transmitted signal
% %and normalize on the dB scale
% sweep_correlation = xcorr(noisy_signal, sweep);
% normalized_correlation = sweep_correlation/max(sweep_correlation);
% sweep_correlation_db = 10*log(abs(normalized_correlation));
% 
% %generate time vector to plot correlation
% correlation_time = transpose(0:length(sweep_correlation_db)-1);
% correlation_time = correlation_time./fs;
% 
% f3 = figure;
% f3.Position = [0, 1320, 1200, 600];
% plot(correlation_time, sweep_correlation_db);
% title('2-ray multichannel simulation correlation');
% ylabel('normalized magnitude (dB)');
% xlabel('time (s)');
% % xlim([0 40]);
% xlim([9.99 10.01]);
% ylim([-75 0]);