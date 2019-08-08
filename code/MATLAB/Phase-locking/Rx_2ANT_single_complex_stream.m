%% Header
% Take Complex stream from binary and correlate 

clear variables
fs = 4e6;
Ax1 = 'A7';
Ax2 = 'A8';

%% Main

%Read in rx_data
%A:AB

rx_data=read_complex_binary(strcat('Binaries/rx_data_complex_',Ax1,Ax2));

R1 = real(rx_data)';
R2 = imag(rx_data)';

%% Read in transmitted and received signals
Tx=read_complex_binary('Binaries/TX_300k_1s_4e6fs_lfm')';

%% Correlate R1
fprintf('Correlating R1\n')
sweep_correlation_R1 = xcorr(Tx,R1);
sweep_correlation_R1 = sweep_correlation_R1(1:length(Tx)+length(R1));
sweep_correlation_db_R1 = 10*log(abs(sweep_correlation_R1));     % convert correlation to dB

%% Correlate R2
fprintf('Correlating R2\n')
sweep_correlation_R2 = xcorr(Tx,R2);
sweep_correlation_R2 = sweep_correlation_R2(1:length(Tx)+length(R2));
sweep_correlation_db_R2 = 10*log(abs(sweep_correlation_R2));     % convert correlation to dB

%% Generate time vector the length of the correlations
correlation_time = transpose(0:length(sweep_correlation_db_R1)-1);
correlation_time = correlation_time./fs;

write_complex_binary(sweep_correlation_db_R1,strcat('Binaries/',Ax1))
write_complex_binary(sweep_correlation_db_R2,strcat('Binaries/',Ax2))
% write_complex_binary(correlation_time,'Binaries/Time')

% %%Plotting
% figure(1)
% plot(correlation_time,sweep_correlation_db_R1)
% title(Ax1)
% axis([0 3 -20 30])
% 
% figure(2)
% plot(correlation_time,sweep_correlation_db_R2)
% title(Ax2)
% axis([0 3 -20 30])
% 
% figure(3)
% hold on
% plot(correlation_time,sweep_correlation_db_R1)
% plot(correlation_time,sweep_correlation_db_R2)
% title(strcat('Overlaid ', Ax1,' & ',Ax2))
% legend(Ax1,Ax2)
% axis([0 3 -20 30])
