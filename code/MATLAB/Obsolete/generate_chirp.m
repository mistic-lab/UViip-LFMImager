fs = 1e6;                  % sampling frequency
T = 5;                      % length of sweep
t = 0:1/fs:T;   
f1 = 0;                     % initial frequency
f2 = 5e5;                   % final frequency
sweep = chirp(t, f1, T, f2);%./2;    % generate sine sweep
sweep_x2 = [sweep sweep]; % put 2 of them back to back
sweep_x2 = sweep_x2./(max(abs(sweep_x2))+.001);

plot(sweep_x2)
% spectrogram(sweep_x2, 1024, 840, 512, fs);

% delay = round(6.56e-4*fs);  % calculate 0.656ms delay in samples
% delayed_sweep = [sweep_x2(delay+1:length(sweep_x2)) sweep_x2(1:delay)]; % generate delayed (skywave) signal
% weight = 0.0158;        % relative weight of weaker signal
% two_ray_signal = sweep_x2 + weight * delayed_sweep;% generate 2-ray simulated received signal
% two_ray_signal = two_ray_signal./(max(abs(two_ray_signal))+.001);
% max(two_ray_signal)
% 
% 
% % write received signal and original sin sweep to .wav file
% %audiowrite('./2_ray.wav', two_ray_signal, fs);
audiowrite('./double_sweep.wav', sweep_x2, fs);
% audiowrite('./sweep.wav', sweep, fs);

