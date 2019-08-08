% *************************************************************************
% * Program to Simulate LFM Radar Model                                   *                                           *
% * Programmed by:  Ahmed Youssef                                         *
% * Revision 0:  02-11-2016                                               *
% * Revision 1:  10-12-2016 adding GO-CFAR anf fixed threshold            *
% * Revision 2:  05-02-2017 adding window function                        *
% * Revision 3:  10-03-2017 put jamming technique (Convolution noise)     *
% * Revision 4:  10-05-2017 SO-CAFR and OS-CFAR                           *
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
B1=75/15;
N_PRI=16;
u=B/(taw);
wo=2*pi*f0;
A=1;
A1=1;
Cell_no=9;
R=35.763e3;
fd=7000; %max test fd=99500 max fd in the world 50000--> x15-north-America
CFAR=1;
%% Basic Data and dB conversions
c=3e8;
fc=3e9;
lamda=c/fc;
%% Signal to noise ratio (SNR), and jamming to signal ratio (JSR)
SNR=0;%dB
JSR=0;%dB
%% Basic Calculation
t_taw=(0:Ts:taw-Ts)';
t_Tr=0:Ts:Tr-Ts;
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
noise=wgn(samp_Tr*N_PRI,1,0,'dBW','complex');
CPI=D_eff.*(CPI_chirp);
%% Convlution noise jamming
if JSR==0
    xj=0;
else
    res_jam=100e-9; % Resonpse of jammer to detect Radar signal and retransmite
    chirp_jam=exp(1j*pi*u.*((t_taw).^2));
    S= (taw/2)+td+res_jam:Tr:Tr*N_PRI+td+res_jam;
    jam_pulses =pulstran(t_cpi,S,'rectpuls',taw);
    %---------------------------------------------------
    %Note for using Pulse train: pulsetran(x,y,'z',d);
    %x: Number of point. Example:t=0:.0001:1
    %y: 'very important' it contains (start:step:end)
    %start:must be (end/2)
    %step & end : how many pulse do you want
    %z: shape of signal. Example 'rectpuls'
    %d: width of the pulse
    %--------------------------------------------------
    %xj_f= wgn(samp_Tr*N_PRI,1,0,'dBW','complex').*jam_pulses;
    xj_f= wgn(samp_Tr*N_PRI,1,0,'dBW','complex').*abs(CPI_chirp);
    
    xj=ifft(conj(fft((chirp_jam),samp_Tr*N_PRI)).*fft(xj_f));
    xj= db2pow((SNR+JSR)/2)*(xj/sqrt(var(xj)));
end
CPI=db2pow(SNR/2)*CPI+xj+noise;
%% window function
W_ID=1;
if W_ID == 1
    w = ones(length(t_taw),1);
    disp('The window used is rectangular window')
elseif W_ID == 2
    w = hanning(length(t_taw));
    disp('The window used is hanning window')
elseif W_ID == 3
    w = hamming(length(t_taw),'periodic');
    disp('The window used is hamming window')
elseif W_ID == 4
    w = bartlett(length(t_taw));
    disp('The window used is bartlett window')
elseif W_ID == 5
    w = blackman(length(t_taw),'periodic');
    disp('The window used is blackman window')
elseif W_ID == 6
    w = blackmanharris(length(t_taw),'periodic');
    disp('The window used is blackmanharris window')
elseif W_ID == 7
    w = chebwin(length(t_taw),150);
    disp('The window used is cheby window')
elseif W_ID == 8
    w = barthannwin(length(t_taw));
    disp('The window used is barthann window')
elseif W_ID == 9
    w = bohmanwin(length(t_taw));
    disp('The window used is bohman window')
elseif W_ID == 10
    w = kaiser(length(t_taw),15.5713);
    disp('The window used is kaiser window')
end
%% creating conj.and reverse of signal and then add filter delay to be causal
S=(A*exp(-1j*pi*u.*((taw-t_taw).^2)));
S=S.*w;
S_spec=(fft((S),length(t_pri)));
%% Matched filter
op_pc=[];
for i=1:N_PRI
    match_op=1.1769*(fft(CPI((i-1)*samp_Tr+1:i*samp_Tr),length(t_pri))).*(S_spec/max(S_spec));
    op_pc1=((ifft((match_op))));
    op_pc=[(op_pc) ; op_pc1];
    if i==1
        op_pc_max=max(op_pc); % to grab the first max for simulation seek
    end
end
% %% Matched filter
% %------creating conj. and reverse of signal add delay to be causal-------
% f=linspace(-fs/2,fs/2,samp_Tr);
% f = -fs/2+(0:samp_Tr-1)*fs/samp_Tr;
% tp = (fftshift(exp(-1j*2*pi*f*taw)))'; % effect of filter delay
% sig2=(A*exp(1j*pi*u.*((t_taw).^2)));
% sig1=sig2.*w;
% S_spec=(fft((sig1),samp_Tr));
% S_spec=conj(S_spec).*tp;
% op_pc=[];
% for i=1:N_PRI
%     match_op=1.1769*(fft(CPI((i-1)*samp_Tr+1:i*samp_Tr),length(t_pri))).*(S_spec)/max((S_spec));
%     op_pc1=((ifft((match_op))));
%     op_pc=[op_pc; op_pc1];
%     if i==1
%         op_pc_max=max(op_pc); % to grab the first max for simulation seek
%     end
% end
%% Moving Target Detection - MTD
x1=[];
for m=1:N_PRI
    x=(transpose((op_pc((m-1)*samp_Tr+1:m*samp_Tr))));
    x1=[(x1) ; (x)];
end
Y_FF =(1/sqrt(N_PRI))*(fft(x1,N_PRI));
%%  Constant False Alarm Rate - Cfar
for cell =4 % Resoultion of compressed signal=  floor(fs/B)
    window=8;
    window_length=window*cell;
    full_window=2*window_length+3;
    a1=1;
    in=0;
    b1=ones(1,window*cell)/window_length;
    a2=1;
    thr=[];
    b2=[zeros(1,((window_length+3))) ones(1,(window*cell))/window_length];
    for i=Cell_no;%1:N_PRI
        in=in+1;
        thr_del=zeros(1,full_window);
        if CFAR==4  %OS
            Y_FF(i,:)=abs(Y_FF(i,:));
            optimi_os_cell=floor((2*window_length)*(3/4));
            %TC(1:window_length+1)=0; %delay because of w_L and GC
            %TC(window_length+2:length(Y_FF(i,:)))=Y_FF(i,1:length(Y_FF(i,:))-(window_length+1)); %OS parameters
            TC_del=zeros(1,window_length+1) ; %delay because of w_L and GC
            TC=[TC_del Y_FF(i,1:length(Y_FF(i,:))-(window_length+1))] ;
            for ii=1:(length(Y_FF(i,:))-full_window)
                w_l=Y_FF(i,ii:ii+window_length-1);
                w_r=Y_FF(i,ii+window_length+3:ii+2*window_length+2);
                window_os=[w_l w_r] ;
                K_win=sort(window_os);
                T_pf_os=4.7119;
                OS_point(ii)=K_win(optimi_os_cell)*T_pf_os;
            end
            thr(in,:)=[thr_del OS_point];
            sig_d(i,:)=[zeros(1,((window_length+1))) (Y_FF(i,(1:(samp_Tr)-((window_length+1))),:))];
            
        else
            y1(i,:)=filter(b1,a1,abs(Y_FF(i,:)));
            y2(i,:)=filter(b2,a2,abs(Y_FF(i,:)));
            if CFAR==1 %CA-CFAR
                thr(i,:) = (y1(i,:)+(y2(i,:)))/2;
            elseif CFAR==2 %GO
                thr(i,:)= max(y1(i,:),y2(i,:));
            elseif CFAR==3 %GO
                thr(i,:)= min(y1(i,:),y2(i,:));
            elseif CFAR==5 %fixed
                thr(i,:)=ones(1,length(y1));
            end
            sig_d(i,:)=[zeros(1,((window_length+1))) (Y_FF(i,(1:(samp_Tr)-((window_length+1))),:))];
        end
    end
    if CFAR==1 %CA-CFAR
        Threshold=5.3223*abs(thr(Cell_no,:));%5.1656  for  9  4.6570
    elseif CFAR==2 %GO
        Threshold=4.9452*abs(thr(Cell_no,:));
    elseif CFAR==3 %SO
        Threshold=7.1349*abs(thr(Cell_no,:));
    elseif CFAR==4 %OS
        Threshold=thr(in,:);
    elseif CFAR==5 % Fixed
        Threshold=1.7619*abs(thr(Cell_no,:));
    end
    Max=max(max(abs(sig_d(Cell_no,:))),max(Threshold));
    Max1=max(abs(sig_d(Cell_no,:)));
    Max2=max(Threshold);
    Max1=1;
    Max2=1;
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
    %     figure(1)
    %     subplot(221),plot(t_taw,(real((S_chirp)))),title('Single Chirp')
    %     xlabel({'Time [Sec]','(a)'},'fontweight','bold'),ylabel('Normalized Amplitude')
    %     subplot(222),plot(t_cpi,real(CPI)),title('one CPI (16 pulses)')
    %     xlabel({'Time [Sec]','(b)'},'fontweight','bold'),ylabel('Normalized Amplitude')
    %     freq=0:1/taw:(length(S_chirp)/taw)-1;
    %     subplot(223),plot(freq,abs(fft(S_chirp))/max(abs(fft(S_chirp)))),
    %     axis([0 3e7 0 1]),grid,title('Spectrum of single chirp')
    %     xlabel({'Frequency [Hz]','(c)'},'fontweight','bold'),ylabel('Normalized Amplitude')
    %     freq=0:1/(16*Tr):(length(CPI_chirp)/(16*Tr))-1;
    %     subplot(224),plot(freq,abs(fft(CPI_chirp))/max(abs(fft(CPI_chirp)))),
    %     axis([0 3e7 0 1]),grid,title('Specturm of one CPI (16 pulses)')
    %     xlabel({'Frequency [Hz]','(d)'},'fontweight','bold'),ylabel('Normalized Amplitude')
    %     figure(2)
    %     subplot(121),plot(t_cpi,abs(op_pc/max(op_pc))),grid,
    %     title('MF output for one CPI')
    %     xlabel('Time'),ylabel('Normalized Amplitude')
    %     subplot(122),plot(t_cpi,abs(op_pc/max(op_pc))),
    %     axis([3.2e-4 3.5e-4 0 1]),grid,title('MF output for one pulse')
    %     xlabel('Time'),ylabel('Normalized Amplitude')
    %     figure(3)
    %     sig_db=20*log10((abs(op_pc)/(abs((op_pc_max)))));
    %     plot(t_cpi,sig_db),axis([3.36e-4 3.41e-4 -100 0]),
    %     title('MF Output for one pulse in dB'),
    %     xlabel('Time'),ylabel('Amplitude in dB')
    %     figure(4)
    %     X_range = (c/(2*1e3))*(0:Ts:Tr-Ts);
    %     Y_freq = (1/1e3)*(0:N_PRI-1)/(Tr*N_PRI);
    %     mesh(X_range,Y_freq,abs(Y_FF.^2)/max(max(abs(Y_FF).^2))),%view([-90,0])
    %     xlabel('Range(Km)'),
    %     ylabel('Frequency(Khz)'),
    %     zlabel('Normalized Amplitude'),title('MTD Output')
    figure(5)
    plot(t_pri,(abs(sig_d(Cell_no,:)))/Max1,t_pri,Threshold/Max2),
    %axis([3.36e-4 3.42e-4 0 1]),grid,
    xlim([3.36e-4 3.42e-4]),grid,
    title('CFAR output'),xlabel('Time'),ylabel('Normalized Amplitude')
    %     hold on
    %% test
    %       figure(6)
    %     plot(thr(1,:))
    %     hold on
end
%% plot in paper
% figure(1)
% t_taw=(-taw/2:Ts:(taw/2)-Ts)';
% S_chirp=(A*exp(1j*pi*u.*((t_taw).^2)));
% subplot(211),plot(t_taw*1e6,(real((S_chirp)))),title('Single Chirp')
% xlabel({'Time [\mu sec]','(a)'},'fontweight','bold'),ylabel('Normalized Amplitude')
% l_f=(length(S_chirp)/taw);
% freq=(-l_f/2):1/taw:(l_f/2)-(1/taw);
% subplot(212),plot(freq*1e-6,abs(fftshift(fft(S_chirp)))/max(abs(fft(S_chirp)))),
% axis([-30 30 0 1]),grid,title('Spectrum of single chirp')
% xlabel({'Frequency [MHz]','(b)'},'fontweight','bold'),ylabel('Normalized Amplitude')
% figure(2)
% subplot(211),plot(t_cpi,abs(op_pc/max(op_pc))),
% axis([3.36e-4 3.41e-4 0 1]),grid,title('MF output for one pulse')
% xlabel({'Time [sec]','(a)'}),ylabel('Normalized Amplitude')
% sig_db=20*log10((abs(op_pc)/(abs((op_pc_max)))));
% subplot(212),plot(t_cpi,sig_db),axis([3.36e-4 3.41e-4 -100 0]),
% title('MF Output for one pulse in dB'),
% xlabel({'Time [sec]','(b)'}),ylabel('Amplitude in dB'),grid
% figure(3)
% X_range = (c/(2*1e3))*(0:Ts:Tr-Ts);
% Y_freq = (1/1e3)*(0:N_PRI-1)/(Tr*N_PRI);
% mesh(X_range,Y_freq,abs(Y_FF.^2)/max(max(abs(Y_FF).^2))),%view([-90,0])
% xlabel('Range(Km)'),
% ylabel('Frequency(Khz)'),
% zlabel('Normalized Amplitude'),title('MTD Output')
% figure(5)
% plot(t_pri,(abs(sig_d(Cell_no,:)))/Max,t_pri,Threshold/Max),
% axis([3.36e-4 3.42e-4 0 1]),grid,
% title('CFAR output'),xlabel('Time [sec]'),ylabel('Normalized Amplitude')
% legend('Radar Signal','CFAR Detector')
% figure(6)
%  subplot(211),plot(t_cpi,real(CPI)),grid
%  xlabel({'Time [sec]','(a)'}),ylabel('Amplitude')
%  xlim([1e-4 4e-4])
%  subplot(212),plot(t_cpi,abs(op_pc)),grid
%  xlabel({'Time [sec]','(b)'}),ylabel('Amplitude')
%  xlim([1e-4 5e-4])
%%  sepctrogram jamming
% samp_td=floor(fs*td);
% %[S,F,t,P] =
% spectrogram(real(xj(samp_td:samp_td+samp_taw)),512,500);
%         t=taw*(0:2/length(t):(2-2/length(t)));
%         F=B*(0:1/length(F):(1-1/length(F)));
%         surf(t,2.25*F,10*log10(P/max(max(P))),'edgecolor','none');% axis tight;
%         xlabel('Time [sec]','FontSize',12,'FontWeight','bold');
%         ylabel('Frequency [Hz]','FontSize',12,'FontWeight','bold');
%         view(0,90);
%         axis([0 200e-6 0 15e6]);
%         colorbar
%% spectrogram signal
% samp_td=floor(fs*td);
% [S,F,tt,P] = spectrogram(real(S_chirp(1:samp_taw,1)),512,500,1*samp_taw);
%         tt=taw*(0:1/length(tt):(1-1/length(tt)));
%         F=B*(0:1/length(F):(1-1/length(F)));
%         surf(tt,2.25*F,10*log10(P/max(max(P))),'edgecolor','none'); axis tight;
%         xlabel('Time [sec]','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
%         ylabel('Frequency [Hz]','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
%         view(0,90);
%         axis([0 100e-6 0 15e6]);
%         colorbar
%%
% plot(t_taw,(real((S_chirp)))),title('Single Chirp')
% xlabel({'Time [Sec]','(a)'},'fontweight','bold'),
% ylabel('Normalized Amplitude')
% xlim([0 10e-6])













