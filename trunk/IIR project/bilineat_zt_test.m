clear;close all;clc

fp = 1e3;
fa = 2e3;

fs = 10e3;

wp = 2*pi*fp/fs;
pw_Omega_p = 2*fs*tan(wp/2);

wa = 2*pi*fa/fs;
pw_Omega_a = 2*fs*tan(wa/2);

Ap = 3;
Aa = 10;

N_freq = 1024*64
W = 0: pi/N_freq: pi-pi/N_freq;

butter_order = buttord(pw_Omega_p, pw_Omega_a, Ap, Aa, 's');
[B,A] = butter(butter_order, pw_Omega_p, 's');

[Bz, Az] = bilinear(B, A, fs);

% H  = freqz(B ,A , 1024*64);
Hz = freqz(Bz,Az, W);

% subplot 211 
% plot(abs(H))

% subplot 212
plot(W*fs/2/pi,abs(Hz))
