%% Header
% Correlates and plots the transmitted and received signals
% Default is to correlate with complex_sweep_0.5s.wav

%% Read in transmitted and received signals
% <<<<<<< Updated upstream
[transmitted_sweep, fs] = audioread(Tx_file);
% =======
% [transmitted_sweep, fs] = audioread('/home/peter/UViip/MATLAB/Transmission-files/complex_sweep_0.5s.wav');
% [transmitted_sweep, fs] = audioread('/home/peter/UViip/MATLAB/Transmission-files/100kHz_0.5s_fs200k.wav');
% >>>>>>> Stashed changes

disp('Use the UI to select the received WAV file');
[received_filename, received_path] = uigetfile('*.wav','Select received WAV file')
[received_sweep, fs] = audioread(strcat(received_path, received_filename));

%% Correlate the signal
sweep_correlation = xcorr(received_sweep(:, 1) + 1i.*received_sweep(:, 2),...
    transmitted_sweep(:, 1) + 1i.*transmitted_sweep(:, 2));   % cross-correlate the signals
sweep_correlation = sweep_correlation(length(received_sweep)-length(transmitted_sweep)+1:...
    end-length(transmitted_sweep)+1);   % trim the correlation signal to remove zero padding
normalized_correlation = sweep_correlation/1;%max(sweep_correlation);  % normalize the correlation to 1
sweep_correlation_db = 10*log(abs(normalized_correlation));     % convert correlation to dB

%% Generate time vector the length of the correlation
correlation_time = transpose(0:length(sweep_correlation_db)-1);
correlation_time = correlation_time./fs;



