%% Header
% Calls complex_correlate.m which requests a Rx file to correlate with
% then plots both the full experiment and the max peak

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
saveas(1, fullfile(received_path,[received_filename '_full']), 'png'); %Saves it in png format in the folder where the Rx file is

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
xlim([time-0.001 time+0.004]);
grid minor
% Save plot with maximum correlation zoomed in
saveas(2, fullfile(received_path, [received_filename '_max-peak']), 'png'); %Saves it in png format in the folder where the Rx file is

