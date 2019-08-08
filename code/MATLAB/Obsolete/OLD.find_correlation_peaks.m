%% Header
% Calls complex_correlate.m which requests a Rx file to correlate with
% then generates a vector of skywave and groundwave points


%% Find skywave and groundwave points
[peaks, locations] = findpeaks(sweep_correlation_db,correlation_time, 'MinPeakHeight', 15, 'MinPeakDistance', 0.001 ); %index all points of at least 15dB


%% Plotting
figure(1); clf;
plot(correlation_time, sweep_correlation_db);
title('full experiment');
ylabel('magnitude (dB)');
xlabel('time (s)');
ylim([-10 60]);

hold
plot(locations, peaks,'x')





