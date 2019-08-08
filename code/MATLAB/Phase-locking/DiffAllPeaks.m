function [samples_between_peaks]=DiffAllPeaks(R1,R2, fs)

fprintf('DIFFALLPEAKS: ');

% Rx_Signal_1=abs(Rx_Signal_1);
% Rx_Signal_2=abs(Rx_Signal_2);

% Difference between delayed and advance recieved signal

% diff1=(Rx_Signal_1)-[(Rx_Signal_1(2:end)); 0];

% diff2=(Rx_Signal_1)-[0; (Rx_Signal_1(1:end-1))];

%Hilbert method
R1_h = hilbert(R1);
R2_h = hilbert(R2);
phase_rad = angle( R1_h ./ R2_h);

figure(1)
plot(phase_rad);
title Hilbert-method

%Cross-correlation method
[c,lag]=xcorr(R1,R2);
[maxC,I]=max(c);
lag = abs(lag(I));
t_lag = lag/fs;
speed_in_coax = (5.05e-9)^-1; %m/s
cable_difference = t_lag*speed_in_coax;


% FFT method
f1 = fft(R1);
f2 = fft(R2);
phase_rad = angle(f1(3e6)/f2(3e6));

figure(2)
plot(phase_rad)



% 1 means local peak

Binary=and(sign(real(R1))>0,sign(real(R2))>0);

% Determine the local maximum 

max1=find(Binary>0);

% Find difference between each peak in samples
samples_between_peaks=max1(2:end)-max1(1:end-1);

figure(3)
plot(samples_between_peaks)

% Result
fprintf('done\n');
end