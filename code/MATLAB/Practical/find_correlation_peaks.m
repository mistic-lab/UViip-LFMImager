
%% Header
% This script finds all of the ground wave points in the correlation
% complex_correlate.m must be run first

%% Initializations
if exist('fignum','var') == 0 fignum = 1; end

%% Plot experiment (to be overlaid with x's
figure(fignum); clf;
hold on
plot(correlation_time, sweep_correlation_db);
title('Are the x''s on the peaks?');
ylabel('magnitude (dB)');
xlabel('time (s)');
ylim([-40 40]);

%% Find groundwave peaks (ask for user confirmation)
count = 1;
correct_peaks = 0;
while correct_peaks == 0
    % Find pwr level which isolates all peaks
    if count == 1 %only run this once
        non_neg_sweep_correlation_db(sweep_correlation_db<0)=0; %makes all negative values 0
        minPWR = ((max(sweep_correlation_db)-mean(non_neg_sweep_correlation_db))/2) + mean(non_neg_sweep_correlation_db)
    end
    count = count + 1; % ensures that the default minPWR isn't used again if it wasn't right the first time
    
    % Locate peaks based on minPWR
    [peaks, locations] = findpeaks(sweep_correlation_db,correlation_time, 'MinPeakHeight', minPWR, 'MinPeakDistance', 0.2 );

    %Plot them over the main plot
    figure(fignum)
    peaksplot(count) = plot(locations, peaks,'x'); %passed into a variable in order to later hide it if it was wrong
    
    % Check whether the peaks have been appropriately found
    prompt = 'Are the peaks correctly indicated? [y/n] \n~>';
    correct_peaks = input(prompt,'s');
    
    if correct_peaks == 'n'
        set(peaksplot(count),'Visible','off') % Hiding the incorrect peaks
        prompt1 = 'Pick a new minPWR which should isolate all of the peaks \n~>';
        minPWR = input(prompt1);
    end
    
    figure(fignum);
    clf; %clear the figure used for these checks
end

%% Store the groundwave magnitudes and times
groundwave_peakvalue = peaks(2:end-1); % eliminating the first and last because they could have less that (-1000:4000) samples surrounding them
groundwave_timevalue = locations(2:end-1);
for i = 1:length(groundwave_peakvalue) 
    groundwave_timeindex(i) = find(correlation_time == groundwave_timevalue(i));      %finds index of the peaktime in correlation_time
    groundwave_peakindex(i) = find(sweep_correlation_db == groundwave_peakvalue(i)); %finds index of the peakvalue in sweep_correlation_db
end




