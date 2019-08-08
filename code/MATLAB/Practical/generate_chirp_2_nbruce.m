clear all; close all;

fs = 4e6;                  % sampling frequency
T = 0.5;                      % length of sweep
t = 0:1/fs:T;
f1 = 5e6;                     % initial frequency
f2 = 6e6;                   % final frequency


initial_phase_shift = 1;% exp(1j*pi/4);
shift_entire_spectrum =1;% exp(-j*(pi*((f2-f1)/T).*t));
sweeping_spectrum = exp(1j*(pi*((f2-f1)/T).*t.^2));

sweep  = initial_phase_shift.*sweeping_spectrum.*shift_entire_spectrum;
% sweep = exp(1j*2*pi*(f1.*t + ((f2-f1)/(2*T)).*t.^2))
%% Sanity check
% sweepR = cos((pi*((f2-f1)/T).*t.^2));
% sweepI = sin((pi*((f2-f1)/T).*t.^2));
% sweepC=sweepR+j*sweepI;

%% Make double sweep
% sweep_x2 = [sweep sweep]; % put 2 of them back to back

%% Don't understand yet
% complex_sweep = transpose([real(sweep); imag(sweep)]); I don't quite get
% this one yet. Why is there a semicolon in there?

%% Write to file
% audiowrite('../Transmission-files/double_sweep_nbruce.dat', sweep_x2, fs);
% 
% audiowrite('../Transmission-files/sweep_nbruce.dat', sweep, fs);

%% Plotting
%  x1=sweep; %just because I copied all this stuff in from the 517
% assignment
% 
% figure(1)
% plot(real(x1))
% 
% figure(2)
% spectrogram(x1, 2024, 2023, 512, fs,'yaxis');
% 
% figure(3)
% plot3(t, real(x1), imag(x1), 'LineWidth',2)
% hold on
% plot3(t, real(x1), zeros(size(x1))-1.5,'r-')
% plot3(t, zeros(size(x1))-2, imag(x1),'g-')
% hold off
% grid on
% axis([0.24  0.26    -2  2    -1.5  1.5])
% % to view the flipping from pos to neg (set fs=1e4)
% % axis([0.24  0.26    -2  2    -1.5  1.5])
% view([-125  30])
% xlabel('Time', 'Rotation',-30)
% ylabel('Real Axis', 'Rotation',10)
% zlabel('Imag Axis')
