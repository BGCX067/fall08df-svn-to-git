% Digital Filtering Project - FIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A heterodyne communication system is designed that oversample samples the
% IF signal at a rate of fs Sa/s.  You are designing a baseband
% communication filter (FIR) that is used to recover a signal whose
% bandwidth is 0.15fs.  Specifically you are to design a set of 31st order
% lowpass linear phase FIRs that has a normalized passband frequency edge
% at wp=0.3pi, stopband beginning at ws=0.4pi, and an arbitrary transition
% band (don�t care).  The design specifications are graphically interpreted
% below on a normalize frequency axis.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all;

filter_order_matlab = 30; %30th order filter for 'rest' of the world
N = filter_order_matlab;
N_fft = 512;
w = [0:2/N_fft: 1 - 1/N_fft];
wp = 0.3;
ws = 0.4;

Hd = w<wp;
mask_pb = w<wp;
mask_sb = w>ws;
mask = mask_pb | mask_sb;

F = [0 wp ws 1];
A = [1 1 0 0];
W = [1 1];
h_uni = firpm(N, F, A, W);

N_symbols = 10000; %total number of symbols to test;

fs      = 1;
fc      = 0.125;         % taking fc/fs = 1/4;
fQPSK   = 0.125;        % taking fQPSK/fs = 0.125;

x_k = [];               % x_k is the time series. x_k = W(k - n N_s)[Acos(2pi fc/fs k + phi_1n)
                        %                                          + Asin(2pi fc/fs k + phi_2n)]
                        % N_s is the number of samples in one QPSK symbol =
                        % N_s = (1/fQPSK)/(1/fs) = fs/fQPSK;

                        
%
N_s = fs/fQPSK;
k = [1:N_s];
A = 1;
re = round(rand(1,N_symbols));
im = round(rand(1,N_symbols));
for i = 1:N_symbols
    phi_1n = pi*re;
    phi_2n = pi*im;
    temp_sym = A*cos(2*pi*fc/fs *k + ones(1,length(k))*phi_1n(i)) +...
        A*sin(2*pi*fc/fs *k +ones(1,length(k))*phi_2n(i));
    x_k = [x_k temp_sym];
end

x_k_tx =x_k;

% stem(x_k)
% h_uni = [1 0 0 0 ];
ber=[]; snr =[];
for SNR_dB = 0:15
%     SNR_dB = 10;
    SNR_lin = 10^(SNR_dB/10);
    signal_power = sum(x_k_tx.^2)/length(x_k_tx);
    noise_variance = signal_power/SNR_lin;

    n = wgn(1, length(x_k), 1, 'linear')*sqrt(noise_variance);
    x_k = x_k_tx + n;

    y_k_unsync = conv(h_uni, x_k);


    offset = 0; % debugging parameter
    y_k = y_k_unsync(16 + offset : length(y_k_unsync) );

    re_hat = []; im_hat = [];%demodulated symbols go here

    for i = 1:N_symbols
        temp_sym = y_k(1 + N_s*(i-1): i*N_s);
        re_hat_temp =  (1 - sign(sum(temp_sym.*cos(2*pi*fc/fs *k ))))/2;
        im_hat_temp =  (1 - sign(sum(temp_sym.*sin(2*pi*fc/fs *k ))))/2;
        re_hat = [re_hat re_hat_temp];
        im_hat = [im_hat im_hat_temp];
    end


    total_bit_errs = sum(abs(re-re_hat)) + sum(abs(im-im_hat));
    SNR_dB
    BER = total_bit_errs/N_symbols/2

    ber = [ber BER];
    snr = [snr SNR_dB];
end

semilogy(snr, ber)
grid on
title('BER vs SNR for QPSK transmitted here');
xlabel('SNR in dB')
ylabel('BER');