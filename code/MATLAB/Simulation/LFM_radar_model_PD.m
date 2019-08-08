% *************************************************************************
% * Program to Simulate LFM Radar Model                                   *                                           *
% * Programmed by:  Ahmed Youssef                                         *
% * Revision 0:  02-11-2016                                               *
% * Revision 1:  10-12-2016 adding GO-CFAR anf fixed threshold  
% * modified by Peter Driessen to show correlation of complex sweeps
% *************************************************************************
% Parameter	                       Units              Description
%   taw:                           second           Pulse Duration
%   Tr:			                   second	     Pulse Repition Interval
%   fs:			                  	Hertz	        Sampling Frequency
%   Ts:			                  	Second	        Sampling Time
%   fc:			                  	Hertz	        Carrier Frequency
%   f0:			                  	Hertz	        Starting Frequency
%   B:			                  	Hertz	           Bandwidth
%   N_PRI:			                Scaler	       Number of PRI's per one scan
%   phi_az:			                Radian	           Target Azimuth
%   R:			                  	meter	           Target Range
%   fd:			                  	Hertz	           Target Doppler
%% clearing
close all
clear 
clc
%% Parameter
taw=100e-6;
fs=2^26;
Ts=1/fs;
Tr=.5e-3;
f0=0;
B=15e6; %bandwidth
N_PRI=16;
u=B/(taw);
wo=2*pi*f0;
A=1;
A1=1;
Cell_no=9;
R=35.763e3;
fd=7000;
CFAR=1;
%% Basic Data and dB conversions
c=3e8;
fc=3e9;
lamda=c/fc;
%% Basic Calculation
t_taw=(0:Ts:taw-Ts)';
%Number of samples in Chirp
samp_taw=floor(taw*fs);
%Number of samples in PRI
Tr=floor(Tr*fs)/fs; % floor(Tr)=0.Thus,we have to multiply with factor
samp_Tr=floor(Tr*fs);
%Round trip Time
td=2*R/c;
%Time samples for single and entire CPI
t_cpi=linspace(0,(N_PRI*Tr)-Ts,((N_PRI*Tr))/Ts)';
t_pri=linspace(0,Tr-Ts,Tr/Ts);
%% Creating Chirp signal
S_chirp=(A*exp(1j*pi*u.*((t_taw).^2)));
%% creating complex baseband chirp signal t_taw=(0:Ts:taw-Ts)';
S_cplx=S_chirp.*exp(-1j*pi*B*t_taw).*exp(-1j*pi/4);

% lowpass not needed since shifted signal was complex to begin with
%lowpass = fdesign.lowpass('Fp,Fst,Ap,Ast',freq/(Fs/2),2*freq/(Fs/2),1,80);
%H_d = design(lowpass,'butter');


%% creating CPI Chirp signal
CPI_chirp = zeros(length(t_cpi),1);
for i=0:N_PRI-1
    n = find( (t_cpi >= i*Tr+td) & (t_cpi < i*Tr+td+taw) );
    CPI_chirp(n)= exp(1j*pi*u*(t_cpi(n)-td-i*Tr).^2);
end
%% Creating Doppler effect
D_eff = (exp(1j*2*pi*fd*(t_cpi-td)));
CPI=D_eff.*(CPI_chirp);
%% creating conj. and reverse of signal add delay to be causal
S=(A*exp(-1j*pi*u.*((taw-t_taw).^2))).*exp(1j*pi*B*t_taw);
%S=(A*exp(-1j*pi*u.*((taw-t_taw).^2)));
%S_chirp=(A*exp(1j*pi*u.*((t_taw).^2)));
%S_cplx=S_chirp.*exp(-1j*pi*B*t_taw).*exp(-1j*pi/4);
S_spec=(fft((S),length(t_pri)));
S_cplx_spec=(fft((S_cplx),length(t_pri)));
yr=ifft(S_spec.*S_cplx_spec);
%% Matched filter
op_pc=[];
for i=1:N_PRI
    match_op=1.176948729352729*(fft(CPI((i-1)*samp_Tr+1:i*samp_Tr),length(t_pri))).*(S_spec/max(S_spec));
    op_pc1=((ifft((match_op))));
    op_pc=[(op_pc) ; op_pc1];
    if i==1
        op_pc_max=max(op_pc); % to grab the first max for simulation seek
    end
end



% instantaneous frequency

z=S_cplx;
zd=filter([0 1],1,z); % delay one sample
fi=angle(conj(zd).*z)*fs/(2*pi);





%% Plotting
figure(1)
subplot(221),plot(t_taw,(real((S_chirp)))),title('Single Chirp')
xlabel({'Time','(a)'},'fontweight','bold'),ylabel('Normalized Amplitude')
subplot(222),plot(t_cpi,real(CPI)),title('one CPI (16 pulses)')
xlabel({'Time','(b)'},'fontweight','bold'),ylabel('Normalized Amplitude')
subplot(223),plot(abs(S_spec)/max(abs(S_spec))),
axis([0 3e4 0 1]),grid,title('Spectrum of single chirp')
xlabel({'Frequency','(c)'},'fontweight','bold'),ylabel('Normalized Amplitude')
subplot(224),plot(abs(fft(CPI))/max(abs(fft(CPI)))),
axis([0 5e5 0 1]),grid,title('Specturm of one CPI (16 pulses)')
xlabel({'Frequency','(d)'},'fontweight','bold'),ylabel('Normalized Amplitude')


figure(2) % plot chirp signal complex baseband
plot3(t_taw, real(S_cplx), imag(S_cplx), 'LineWidth',2)
hold on
plot3(t_taw, real(S_cplx), zeros(size(S_cplx))-1.5,'r-')
plot3(t_taw, zeros(size(S_cplx))-2, imag(S_cplx),'g-')
hold off

cn=samp_taw/2;
Npts=cn/5;
cm=cn-Npts;
cp=cn+Npts;
figure(3) % plot chirp signal complex baseband zoom in near zero frequency
plot3(t_taw(cm:cp), real(S_cplx(cm:cp)), imag(S_cplx(cm:cp)), 'LineWidth',2)
hold on
plot3(t_taw(cm:cp), real(S_cplx(cm:cp)), zeros(size(S_cplx(cm:cp)))-1.5,'r-')
plot3(t_taw(cm:cp), zeros(size(S_cplx(cm:cp)))-2, imag(S_cplx(cm:cp)),'g-')
hold off

figure(4)
plot(t_taw,fi)
xlabel('time (seconds)');
ylabel('instantaneous frequency (Hz)');
title('instantaneous frequency of complex baseband chirp signal');

figure(5)
plot(abs(yr));

figure(6)
cn=13300;
Npts=cn/100;
cm=cn-Npts;
cp=cn+Npts;
plot(abs(yr(cm:cp)));

figure(7)
plot(angle(yr(cm:cp)));

figure(8)
% plot correlation complex baseband zoom in near peak
plot3((cm:cp), real(yr(cm:cp)), imag(yr(cm:cp)), 'LineWidth',2)
hold on
plot3((cm:cp), real(yr(cm:cp)), zeros(size(yr(cm:cp)))-1.5,'r-')
plot3((cm:cp), zeros(size(yr(cm:cp)))-2, imag(yr(cm:cp)),'g-')
hold off

% figure(2)
% subplot(121),plot(t_cpi,abs(op_pc/max(op_pc))),grid,
% title('MF output for one CPI')
% xlabel('Time'),ylabel('Normalized Amplitude')
% subplot(122),plot(t_cpi,abs(op_pc/max(op_pc))),
% axis([3.2e-4 3.5e-4 0 1]),grid,title('MF output for one pulse')
% xlabel('Time'),ylabel('Normalized Amplitude')
% figure(3)
% sig_db=20*log10((abs(op_pc)/(abs((op_pc_max)))));
% plot(t_cpi,sig_db),axis([3.36e-4 3.41e-4 -100 0]),
% title('MF Output for one pulse in dB'),
% xlabel('Time'),ylabel('Amplitude in dB')
% figure(4)
% X_range = (c/(2*1e3))*(0:Ts:Tr-Ts);
% Y_freq = (1/1e3)*((0:(N_PRI)-1)*(1/Tr)/(N_PRI));
% mesh(X_range,Y_freq,abs(Y_FF.^2)/max(max(abs(Y_FF).^2))),%view([-90,0])
% xlabel('Range(Km)'),
% ylabel('Frequency(Khz)'),
% zlabel('Normalized Amplitude'),title('MTD Output')
% figure(5)
% plot(t_pri,(abs(sig_d(Cell_no,:)))/Max,t_pri,Threshold/Max),
% axis([3.36e-4 3.42e-4 0 1]),grid,
% title('CFAR output'),xlabel('Time'),ylabel('Normalized Amplitude')
% 




