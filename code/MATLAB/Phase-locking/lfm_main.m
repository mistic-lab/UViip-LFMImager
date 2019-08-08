%% Header
% Run phase-locking test for N210s

clear variables
% delete ./tx_symbol* rx_data

%% Initializations
fs = 4e6;
f1 = 0;
f2 = .3e6;
T = 1;


%% Main
%Generate sweep
[tx_symbol, t] = generate_chirp(fs,f1,f2,T);

%Write tx to binary
write_complex_binary(tx_symbol,'Binaries/TX_300k_1s_4e6fs_lfm');

spectrogram(tx_symbol,[],[],[],fs) % check sweep

%Pause while experiment is run
fprintf('paused\n');
pause

%Read in rx_data
rx_data=read_complex_binary('rx_data');

%Correlate
[correlation_db,correlation,t_correlation] = correlate(fs,tx_symbol,rx_data);

%Find peaks
[FirstTwoPeaks,DiffPeakSample] = DiffPeak(correlation);
time_between_peaks = [DiffPeakSample,0.5]./fs;
speed_in_coax=5.05e-9; %RG6
length_between_peaks = time_between_peaks./speed_in_coax;
fprintf('--> %.0f m +/- %.0f m (%.0f ft +/- %.0f ft) between peaks\n',...,
    length_between_peaks(1),length_between_peaks(2),...,
    3.28084*length_between_peaks(1),3.28084*length_between_peaks(2))

%Plot
plot(abs(correlation))
xlabel('Samples')
ylabel('Impulse Response [dB]')