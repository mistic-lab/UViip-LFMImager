%% Header
% data_condenser.m
% Prepares data already in workspace (from data_analysis_main.m) for
% sharing

%% Save to .MAT file (large file size)
save 20170603-0119-DATA.mat correlation_time plot_magnitudes plot_xtime received_sweep sweep_correlation transmitted_sweep

%% Save to text file
% make everything column data
plot_magnitudes = plot_magnitudes';
plot_xtime = plot_xtime';

% write to matrix
data_matrix = [correlation_time; plot_magnitudes; plot_xtime; received_sweep; sweep_correlation; transmitted_sweep;

% create file
fileID = fopen('20170603-0119-DATA.txt','w');
fprintf(fileID,'%6s %12s %18s %24s %30s %36s\n',...
    'correlation_time',...
    'plot_magnitudes',...
    'plot_xtime',...
    'received_sweep',...
    'sweep_correlation',...
    'transmitted_sweep');
fprintf(fileID,'%6.2f %12.8f\n',data_matrix); %not done this yet
fclose(fileID);