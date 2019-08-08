fs = 2e5;                  % sampling frequency
T = 0.5;                      % length of sweep
t = 0:1/fs:T;   
% f1 = 0;                     % initial frequency
f1 = -5e4;
% f2 = 1e5;                   % final frequency
f2 = 5e4
% sweep = exp(j*(pi*((f2-f1)/T).*t.^2));
sweep = exp(1j*2*pi*(f1.*t + ((f2-f1)/(2*T)).*t.^2));
% spectrogram(sweep, 2048, 840, 512, fs);
% figure;
% plot(abs(fft(sweep)));


complex_sweep = transpose([real(sweep); imag(sweep)]);
audiowrite('./100kHz_0.5s_fs200k.wav', complex_sweep, fs);


