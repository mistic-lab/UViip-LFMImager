e% read in transmitted and received signals
[transmitted_sweep, fs] = audioread('double_sweep_nbruce.wav');
[received_sweep, fs] = audioread('sweep_vertical_tx_inverted_rx_170428');
% [received_sweep, fs] = audioread('received_sweep_virtual_noise.wav');

sweep_correlation = xcorr(received_sweep, transmitted_sweep);   % cross-correlate the signals
sweep_correlation = sweep_correlation(length(received_sweep)-length(transmitted_sweep)+1:...
    end-length(transmitted_sweep)+1);   % trim the correlation signal to remove zero padding
% normalized_correlation = sweep_correlation/max(sweep_correlation);  % normalize the correlation to 1
sweep_correlation_db = 10*log(abs(normalized_correlation));     % convert correlation to dB

%generate time vector to plot correlation
correlation_time = transpose(0:length(sweep_correlation_db)-1);
correlation_time = correlation_time./fs;

% plot and format correlation
plot(correlation_time, sweep_correlation_db);   
title('2-ray multichannel simulation correlation');
ylabel('normalized magnitude (dB)');
xlabel('time (s)');
ylim([-110 10]);