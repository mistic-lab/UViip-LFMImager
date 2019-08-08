%% Header
%   This file calls a series of functions and scripts to analyze, plot and
%   save data from the UViip experiments.

clear all % start fresh

%% Correlate Tx and Rx signals
default_Tx_file = '../Transmission-files/100kHz_0.5s_fs200k.wav';

%Check whether to use the normal sweep
prompt = ['Use default Tx file? [y/n] (default = ',default_Tx_file,') \n~>'];
usedefaultTx = input(prompt,'s');

if usedefaultTx == 'y'
    Tx_file = default_Tx_file;
elseif usedefaultTx == 'n'
    [transmitted_filename, transmitted_path] = uigetfile('*.wav','Select transmitted WAV file');
    Tx_file = [transmitted_path transmitted_filename];
end

disp('Launching correlation ...');
complex_correlate % Call script
disp('Correlation complete.');

%% Find groundwave peaks
disp('Calculating minimum power to isolate peaks ...');
find_correlation_peaks
disp('Peak data stored.');

%% Plot the correlation
prompt = 'Save plots? [y/n] \n~>';
saveplotsbin = input(prompt,'s');

if saveplotsbin == 'y'
    saveplots = 1;
elseif saveplotsbin == 'n'
    saveplots = 0;
end


plot_correlation % Call script
disp('Data analysis complete.');