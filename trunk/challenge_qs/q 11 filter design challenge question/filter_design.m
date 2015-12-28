clear;clc;close all
N=101;
w = [0:2/N:1.998]*pi;

H_abs = w<0.4*pi | w>2*pi-0.4*pi;
H_ang = exp(-(N-1)/2*j*w);
H = H_abs.*H_ang;

h = ifft(H);
subplot(2,2,1)
h_trunc = h(N/2 - 2: N/2 + 2)
stem(real(h_trunc))

subplot(2,2,2)
semilogy(w, abs(fft(real(h_trunc), length(w))));
hold on
semilogy(w, H_abs, 'r-');

subplot(2,2,3)
plot(w, angle(fft(real(h_trunc), length(w))));
grid on;
hold on;
plot(w, angle(H_ang),'r.');

h_4  = [.1,-.3,.4,-.3,.1];
h_3  = [.1,.3,.4,.3,.1];
h_2  = [.1,-.2,.4,-.2,.1];
h_1  = [.1,.2,.4,.2,.1];
h_tool = [0.07, 0.25, 0.33, 0.25, 0.07];
subplot(2,2,4);semilogy(w,abs(fft(h_1, length(w))))
grid on