%% Header
% This script plots both the full experiment and the max peak
% complex_correlate must be run first!

%% Initializations
if exist('fignum','var') == 0 fignum = 1; end
received_filename = received_filename(1:end-4); %this gets rid of the ".wav" for figure/file titling
time_prebuffer      = -1e-3;  % how many s before each groundwave to plot
time_postbuffer     = 8e-3; % how many s after each groundwave to plot
samples_prebuffer   = time_prebuffer*fs;    % puts it in ms and then finds equivalent number of samples
samples_postbuffer  = time_postbuffer*fs;

%% Plot full experiment
figure(fignum); clf;
plot(correlation_time, sweep_correlation_db);
title([received_filename ' full experiment']);
ylabel('magnitude (dB)');
xlabel('time (s)');
ylim([-40 40]);
fullexpfignum = fignum;
fignum = fignum + 1;

%% Plot all ground waves overlaid
figure(fignum); clf;
grid
plot_xtime = time_prebuffer:1/fs:time_postbuffer;
plot_ysep = 0:1/length(groundwave_peakindex):1; % would be nice for this to show correlation_time

for i = 1:length(groundwave_peakvalue)
    plot_magnitudes(i,:) = sweep_correlation_db( groundwave_peakindex(i)+samples_prebuffer:groundwave_peakindex(i)+samples_postbuffer);
    plot_yMAT(i,1:length(plot_xtime)) = plot_ysep(i);
end

plot3(plot_xtime,plot_yMAT(1,:),plot_magnitudes(1,:));
hold on
for i = 2:length(plot_ysep)-1
    plot3(plot_xtime,plot_yMAT(i,:),plot_magnitudes(i,:));
end



title([received_filename ' max peaks overlaid']);
xlabel('time normalized around the highest peak (s)');
ylabel('No meaningful separation value')
zlabel('correlation magnitude (dB)');
xlim([time_prebuffer time_postbuffer]);
ylim([0 1]);
zlim([-40 40]);

grid minor

peaksfignum = fignum;
fignum = fignum + 1;

%% Save plots
if saveplots == 1
    % Save plot of full experiment
    saveas(fullexpfignum, fullfile(received_path,[received_filename '_full']), 'png'); %Saves it in png format in the folder where the Rx file is
    % Save plot of all peaks overlaid
    saveas(peaksfignum, fullfile(received_path, [received_filename '_max-peak']), 'fig'); %Saves it in fig format in the folder where the Rx file is
end
