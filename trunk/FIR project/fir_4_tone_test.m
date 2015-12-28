% Digital Filtering Project - FIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A heterodyne communication system is designed that oversample samples the
% IF signal at a rate of fs Sa/s.  You are designing a baseband
% communication filter (FIR) that is used to recover a signal whose
% bandwidth is 0.15fs.  Specifically you are to design a set of 31st order
% lowpass linear phase FIRs that has a normalized passband frequency edge
% at wp=0.3pi, stopband beginning at ws=0.4pi, and an arbitrary transition
% band (don’t care).  The design specifications are graphically interpreted
% below on a normalize frequency axis.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 4:  Testing (a)	Test your equiripple (Part 3.a) output response to
% sinusoidal inputs at frequencies f0=0.05fs and f0=0.15fs.  Measure the
% input-output phase delay. (b)	Test your equiripple Part 3.a) filter’s
% output response to a unit bound worst case input.  Compare the outcome to
% the worst case gain prediction.


clear; clc; close all

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

fs = 1
f01 = 0.05*fs;
f02 = 0.15*fs;

n = [0:0.1:25];
x1 = cos(2*pi*f01/fs *n);
x2 = cos(2*pi*f02/fs *n);

y1 = conv(h_uni, x1);
y2 = conv(h_uni, x2);

subplot(2,2,1);
plot(x1);
grid on;
hold on
plot(y1, '-rs');
legend('input','output');
title('f_0 = 0.05f_s');
text(1.5,-1.25,'output appears after 15 samples');

subplot(2,2,3)
plot(x2);
grid on;
hold on;
plot(y2,'-rs');
legend('input','output');
title('f_0 = 0.15f_s');
text(1.5,-1.25,'output appears after 15 samples');


%%%%%%%%%%%%%%%
% Worst case  %
%    input    %
%%%%%%%%%%%%%%%
subplot(2,2,2)
x_worst_case_input = sign(h_uni);
worst_case_gain = sum(abs(h_uni))

y_worst_case = conv(h_uni,x_worst_case_input);
stem(x_worst_case_input);
hold on
stem(y_worst_case, 'rs')
Worst_case_max_output = max(y_worst_case)
legend('worst case input', 'output');
text(20, 1.3,'Maximum value of the output is 1.742');
text(20, 1.4,'Worst case gain = 1.742');