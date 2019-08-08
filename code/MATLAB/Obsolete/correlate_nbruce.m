
    % read in transmitted and received signals
    [transmitted_sweep, fs] = audioread('double_sweep_nbruce.wav');
    %[received_sweep, fs] = audioread('received_sweep_coax_nbruce_60.wav');
    [received_sweep, fs] = audioread('../Results/2017_05_11_2346_PHTx_ELWRx/20170511_2346_PHTx_ELWRx.wav');
    
    sweep_correlation = xcorr(received_sweep, transmitted_sweep);   % cross-correlate the signals
    sweep_correlation = sweep_correlation(length(received_sweep)-length(transmitted_sweep)+1:...
        end-length(transmitted_sweep)+1);   % trim the correlation signal to remove zero padding
    sweep_correlation_db = 10*log(abs(sweep_correlation));     % convert correlation to dB
    
    %generate time vector to plot correlation
    correlation_time = transpose(0:length(sweep_correlation_db)-1);
    correlation_time = correlation_time./fs;
    
    % plot and format correlation
    figure
    plot(correlation_time, sweep_correlation_db);
    title('Complex correlation (same file for Tx and Rx)');
    ylabel('normalized magnitude (dB)');
    xlabel('time (s)');
    ylim([-20 140]);



