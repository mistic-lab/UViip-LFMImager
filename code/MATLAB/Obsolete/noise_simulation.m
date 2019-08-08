fs = 10e5;                  %sampling frequency
T = 5;                      %length of sweep
t = 0:1/fs:T;   
f1 = 0;                     %initial frequency
f2 = 5e5;                   %final frequency
sweep = chirp(t, f1, T, f2);    %generate sine sweep
sweep_x2 = [sweep sweep];

noise = wgn(1, 10000002, -100);
signal_weight = 1e-8
sig_noise_ratio = snr(signal_weight*sweep_x2, noise);

noisy_signal = signal_weight*sweep_x2 + noise;
% plot(abs(fft(noisy_signal)));

correlation = xcorr(sweep, noisy_signal);
norm_corr_db = 10*log(abs(correlation ./ max(correlation)));
figure;
plot(norm_corr_db);

title(sprintf('Correlation with snr = %f', sig_noise_ratio));

xlabel('sample');
ylabel('mag');