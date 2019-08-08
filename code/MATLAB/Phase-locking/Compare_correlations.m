%% Header
% Compare correlations

clear variables

num_corrs = 8;
plot_ysep = 0:num_corrs; % number of data vectors
fignum=1;

for i=1:num_corrs
    ant = strcat('A',num2str(i));
    correlations(i,:) = read_complex_binary(strcat('Binaries/',ant));
    plot_yMAT(i,1:length(correlations)) = plot_ysep(i);
end
correlations(i+1,:) = read_complex_binary('Binaries/Time');

%% Plot all ground waves overlaid
grid
figure(fignum)
hold on
for i = 1:num_corrs
    plot3(correlations(end,:),plot_yMAT(i,:),correlations(i,:));
end
title('All test transmissions');
xlabel('Time (s)');
ylabel('Antenna #')
zlabel('Correlation magnitude (dB)');
zlim([-15 40]);
ylim([1 8]);
xlim([1.986,1.993])
view([0 270 0])

hold off
fignum = fignum+1;

for i = 1:num_corrs
    figure(fignum)
    plot(correlations(end,:),correlations(i,:));
    title(strcat('A',num2str(i)))
    ylim([-15 40]);
    xlabel('Time (s)')
    ylabel('Correlation magnitude (dB)')
    fignum = fignum+1;
end




