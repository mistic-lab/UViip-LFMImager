% *************************************************************************
% * Program to Simulate LFM Radar Model                                   *                                           *
% * Programmed by:  Ahmed Youssef                                         *
% * Revision 0:  02-11-2016                                               *
% * Revision 1:  10-12-2016 adding GO-CFAR anf fixed threshold            *
% *************************************************************************
% Parameter	                       Units              Description
%   taw:                           second           Pulse Duration
%   Tr:			                   second	     Pulse Repition Interval
%   fs:			                  	Hertz	        Sampling Frequency
%   Ts:			                  	Second	        Sampling Time
%   fc:			                  	Hertz	        Carrier Frequency
%   f0:			                  	Hertz	        Starting Frequency
%   B:			                  	Hertz	           Bandwdith
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
B=15e6;
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
S=(A*exp(-1j*pi*u.*((taw-t_taw).^2)));
S_spec=(fft((S),length(t_pri)));
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
%% Moving Target Detection - MTD
x1=[];
for m=1:N_PRI
    x=(transpose((op_pc((m-1)*samp_Tr+1:m*samp_Tr))));
    x1=[(x1) ; (x)];
end
Y_FF =(1/sqrt(N_PRI))*(fft(x1,N_PRI));
%%  Constant False Alarm Rate - Cfar
cell = 4; % Resoultion of compressed signal=  floor(fs/B)
a1=1;
b1=ones(1,8*cell);
a2=1;
b2=[zeros(1,(11*cell)) ones(1,(8*cell))];
for i=1:N_PRI
    y1(i,:)=filter(b1,a1,abs(Y_FF(i,:)));
    y2(i,:)=filter(b2,a2,abs(Y_FF(i,:)));
    if CFAR==1 %CA-CFAR
        thr(i,:) = ((y1(i,:))+(y2(i,:)))/(2*cell*8);
    elseif CFAR==2 %GO
        thr(i,:)= max(y1(i,:),y2(i,:));
    elseif CFAR==3 %fixed
        thr(i,:)=ones(1,length(y1));
    end
    sig_d(i,:)=[zeros(1,(9*cell)) (Y_FF(i,(1:(samp_Tr)-(9*cell)),:))];
end
if CFAR==1 %CA-CFAR
    Threshold=5.3223*abs(thr(Cell_no,:));%old-6.442902029133024
elseif CFAR==2
    Threshold=4.9452*abs(thr(Cell_no,:));%4.8263 CAGO
elseif CFAR==3
    Threshold=1.7619*abs(thr(Cell_no,:));
end
Max=max(max(abs(sig_d(Cell_no,:))),max(Threshold));
%% Calculation For MF output
pc_peak=findpeaks(abs(op_pc));
s_sort=sort(pc_peak,'descend');
x=(s_sort(1));
y=(s_sort(18));
sll=pow2db((y/x)^2)
Improvment_factor_MF_dB=20*log10(x)
Ideal_Improvment_factor_dB=pow2db(B*taw)
max_amp_exact=sqrt(B*taw);
max_real=max(abs(op_pc));
error_MF=((max_amp_exact-max_real)/max_amp_exact)*100
%% Calculation for MTD ouput
max_Y_FF=max(max(abs(Y_FF)));
Improvment_factor_MTD=(max_Y_FF/x).^2
error_MTD=((16-15.999972745241624)/16)*100
Improvment_factor_MTD_dB=pow2db(max(max(abs(Y_FF.^2))))-Improvment_factor_MF_dB
%% Calculation for all system
Improvment_factor_all_system=pow2db(max_Y_FF^2)+6.5
%% Ploting
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




