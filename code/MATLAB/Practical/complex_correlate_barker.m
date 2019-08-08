%% Header
% Correlates and plots the transmitted and received signals
clear all;

%% Read in transmitted and received signals
% [transmitted_sweep, fs] = audioread('../Transmission-files/complex_sweep_0.5s.wav');
[transmitted_sweep, fs] = audioread('../Transmission-files/barker.wav');

[received_filename, received_path] = uigetfile('*.wav','Select Received WAV file');
[received_sweep, fs] = audioread(strcat(received_path, received_filename));

%% Correlate the signal
sweep_correlation = xcorr(received_sweep(:, 1),...
    transmitted_sweep(:, 1));   % cross-correlate the signals
sweep_correlation = sweep_correlation(length(received_sweep)-length(transmitted_sweep)+1:...
    end-length(transmitted_sweep)+1);   % trim the correlation signal to remove zero padding
normalized_correlation = sweep_correlation/1;%max(sweep_correlation);  % normalize the correlation to 1
sweep_correlation_db = 10*log(abs(normalized_correlation));     % convert correlation to dB

%% Generate time vector for plotting
correlation_time = transpose(0:length(sweep_correlation_db)-1);
correlation_time = correlation_time./fs;

%% Plot and save correlations
received_filename = received_filename(1:end-4); %this gets rid of the ".wav" for figure/file titling

% Plot of full experiment
figure(1); clf;
plot(correlation_time, sweep_correlation_db);
title([received_filename ' full experiment']);
ylabel('magnitude (dB)');
xlabel('time (s)');
ylim([-10 60]);
% Save plot of full experiment
% saveas(1, fullfile(received_path,[received_filename '_full']), 'png'); %Saves it in png format in the folder where the Rx file is

% Find correlation points
[peaks, locations] = findpeaks(sweep_correlation_db, 'MINPEAKHEIGHT', 15); %index all points of at least 15dB
[val, idx] = max(peaks); %save amplitude and time of the point with highest SNR
sample = locations(idx); %locations is the vector from two lines up
time = correlation_time(sample); %create a vector to plot time
t_ref = floor(time*100)/100; %rounds time elements down
grid

% Plot the peak with maximum correlation
figure(2); clf;
plot(correlation_time, sweep_correlation_db);
title([received_filename ' max peak']);
ylabel('normalized magnitude (dB)');
xlabel('time (s)');
ylim([-10 60]);
% xlim([t_ref t_ref + 0.005]);
xlim([time-0.001 time+0.008]);
grid minor
% Save plot with maximum correlation zoomed in
% saveas(2, fullfile(received_path, [received_filename '_max-peak']), 'png'); %Saves it in png format in the folder where the Rx file is

hold on
yyaxis right

plot(correlation_time, angle(sweep_correlation));



