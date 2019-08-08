%% Header
% Run phase-locking test for N210s

clear variables
% possibilities={'rx_data', 'rx_data_RXA', 'rx_data_RXB'};
% for i=1:numel(possibilities)
%     if exist(possibilities{i},'file')==2
%         delete(possibilities{i})
%     end
% end

%% Main

%Read in rx_data
%A:AB
rx_data=read_complex_binary('Binaries/rx_data_closest');
% rx_phase=read_complex_binary('Binaries/rx_data_phase_latest');

R1 = real(rx_data);
R2 = imag(rx_data);
% rx_phase = real(rx_phase);

%A:A A:B
% R1=read_complex_binary('Binaries/rx_data_RXA_latest');
% R2=read_complex_binary('Binaries/rx_data_RXB_latest');
% rx_mult = read_complex_binary('Binaries/rx_data_mult_latest');

phase_difference = angle(R1./R2);
% phase_difference = angle(rx_data);

% find_value = find(phase_difference > 0);
% find_value = find_value(1);
% 
% l = length(phase_difference);
% num_0 = nnz(phase_difference==0);
% percent_0=num_0*100/l
% num_not = nnz(phase_difference==phase_difference(find_value));
% percent_not = num_not*100/l
% any_others = l - num_0 - num_not



figure(1)
plot(R1)
hold on
plot(R2)
% title Signals
% legend('R1','R2')
% hold off

% figure(2)
% plot(real(R1))
% hold on
% plot(real(R2))
% title Real
% legend('R1','R2')
% hold off
% 
% figure(3)
% plot(imag(R1))
% hold on
% plot(imag(R2))
% title Imaginary
% legend('R1','R2')
% hold off

% figure(4)
plot(phase_difference)
% legend ('R1, 'R2','angle(R1./R2)')
% title Phase-Difference
axis([0 length(R1) -4 4])

% figure(5)
% plot(rx_phase)
% legend('Phase from GR')
% title Phase-Difference

% figure(6)
% plot(rx_mult)
% legend('Mult from GR')
% title Phase-Difference