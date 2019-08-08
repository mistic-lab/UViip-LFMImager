function [correlation_db,correlation,t_correlation] = correlate(fs, tx_symbol,rx_data)
%CORRELATE takes two binary files and correlates them
%   Arguments are tx_sweep vector and rx_binary vector
%   Returns correlation vector and time vector


fprintf('CORRELATE: ');

%% Correlate the signal
correlation = xcorr(rx_data,tx_symbol);   % cross-correlate the signals
correlation = correlation(length(rx_data)-length(tx_symbol)+1:end-length(tx_symbol)+1);   % trim the correlation signal to remove zero padding
correlation_db = 10*log(abs(correlation));     % convert correlation to dB

%% Generate time vector the length of the correlation
t_correlation = transpose(0:length(correlation_db)-1);
t_correlation = t_correlation./fs;

fprintf('done\n');
end

