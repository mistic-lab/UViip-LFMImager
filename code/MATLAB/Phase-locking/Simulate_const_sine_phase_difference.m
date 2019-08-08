%% Header
% Find phase difference between two complex constant sine waves
clear variables


%% Initializations
fs = 200;
T=1;
fc=5;

% fs = 4e6;
% T=0.5;
% fc=.1e6;

% Time vector
t = 0:1/fs:T-1/fs;

%Signal generation FIX THIS
R1=exp(-1j*(2*pi*((fc)/T).*t));
R2_phase_shift = exp(-1j*(pi/4));
R2=R1*R2_phase_shift;

phase_difference = angle(R1./R2);

figure(1)
plot(R1)
hold on
plot(R2)
title Complex
hold off

figure(2)
plot(real(R1))
hold on
plot(real(R2))
title Real
hold off

figure(3)
plot(imag(R1))
hold on
plot(imag(R2))
title Imaginary
hold off

figure(4)
plot(phase_difference)
title Phase-Difference