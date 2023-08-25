function output_modu = modulation(input_frame, index)
% Input_modu: input bit stream (0,1)
% index:  modulation index
%		1---bpsk
%		2---qpsk
%		4---16qam
%		6---64qam
% 	else is error
f_length = length(input_frame)/index;
QAM_input_I = zeros(1,f_length);
QAM_input_Q = zeros(1,f_length);
% note: Matlab index starts from 1
switch index
case 1,
    BPSK_I = [-1 1];    % refer to Table82 on page21 of IEEE802.11a
    QAM_input_I = BPSK_I(input_frame+1);
case 2,
    QPSK_IQ = [-1 1];   % refer to Table83 on page21 of IEEE802.11a
    QAM_input_I = QPSK_IQ(input_frame(1:2:end)+1)/sqrt(2);
    QAM_input_Q = QPSK_IQ(input_frame(2:2:end)+1)/sqrt(2);
case 4,
    QAM_16_IQ = [-3 -1 3 1];    % refer to Table84 on page21 of IEEE802.11a
    QAM_input_I = QAM_16_IQ(input_frame(1:4:end)*2+input_frame(2:4:end)+1)/sqrt(10);
    QAM_input_Q = QAM_16_IQ(input_frame(3:4:end)*2+input_frame(4:4:end)+1)/sqrt(10);
case 6,
    QAM_64_IQ = [-7 -5 -1 -3 7 5 1 3];  % refer to Table85 on page21 of IEEE802.11a
    QAM_input_I = QAM_64_IQ(input_frame(1:6:end)*4+input_frame(2:6:end)*2+input_frame(3:6:end)+1);
    QAM_input_Q = QAM_64_IQ(input_frame(4:6:end)*4+input_frame(5:6:end)*2+input_frame(6:6:end)+1);
end
output_modu = QAM_input_I + j * QAM_input_Q;
