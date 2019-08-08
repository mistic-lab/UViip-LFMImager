%% Header
% Simulate phase-locking test for N210s

clear variables

%% Initializations
fs = 2e6;
f1 = 0;
f2 = fs/4;
T = .25;

length_of_short_path=1.5;   %meters
length_of_long_path=400;    %meters
time_of_short_path=length_of_short_path/3e8;    %seconds
time_of_long_path=length_of_long_path/3e8;      %seconds

time_delay_start_end=0.25;

t_start_to_first = 10*T;
t_start_to_second = 15*T;
t_open_close = t_start_to_first;
t_in_between = t_start_to_second-t_start_to_first;

%% Main
%Generate sweep
[tx_symbol, t] = generate_chirp(fs,f1,f2,T);

%Create rx_data
open_close_padding = complex(zeros(1,time_delay_start_end*fs),0);
in_between_padding = complex(zeros(1,t_in_between*fs),0);

rx_data = cat(2,open_close_padding,tx_symbol,in_between_padding,tx_symbol,open_close_padding);

%Correlate
[correlation_db,t_correlation] = correlate(fs,tx_symbol,rx_data);

%Plot
plot(t_correlation,correlation_db)
