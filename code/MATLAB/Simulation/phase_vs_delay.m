% clear all;

fs = 1e6;                  % sampling frequency
T = 0.5;                      % length of sweep
t = 0:1/fs:T;   
f1 = 0;                     % initial frequency
f2 = 1e5;                   % final frequency

sweep = exp(j*2*pi*(f1.*t + ((f2-f1)/(2*T)).*t.^2));
% spectrogram(sweep, 2048, 840, 512, fs);
% plot(abs(fft(sweep)));

transmitted_sweep = transpose([real(sweep); imag(sweep)]);
% received_sweep = [zeros(100000, 2); transmitted_sweep];
received_sweep = [transmitted_sweep; 0.5.*transmitted_sweep; zeros(100000, 2)];

%% Correlate the signals
% sweep_correlation = xcorr(transmitted_sweep(:, 1) + 1i.*transmitted_sweep(:, 2),...
%     transmitted_sweep(:, 1) + 1i.*transmitted_sweep(:, 2));   % cross-correlate the signals
sweep_correlation = xcorr(received_sweep(:, 1) + 1i.*received_sweep(:, 2),...
    transmitted_sweep(:, 1) + 1i.*transmitted_sweep(:, 2));   % cross-correlate the signals
% sweep_correlation = xcorr(transmitted_sweep(:, 1) + 1i.*transmitted_sweep(:, 2),...
%     received_sweep(:, 1) + 1i.*received_sweep(:, 2));   % cross-correlate the signals
% sweep_correlation = sweep_correlation(length(received_sweep)-length(transmitted_sweep)+1:...
%     end-length(transmitted_sweep)+1);   % trim the correlation signal to remove zero padding
normalized_correlation = sweep_correlation/1;%max(sweep_correlation);  % normalize the correlation to 1
sweep_correlation_db = 10*log(abs(normalized_correlation));     % convert correlation to dB

%% Generate time vector for plotting
correlation_time = transpose(0:length(sweep_correlation_db)-1);
correlation_time = correlation_time./fs;

figure
plot(correlation_time, sweep_correlation_db);
ylim([0 200])

% saveas(1, 'full_half_half_second_delay', 'png');
