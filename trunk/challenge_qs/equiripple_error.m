clear;clc;close all;

N = 1000;
w = [0:2/N:2-1/N];

wc = 0.5
f = [0 wc-0.035 wc+0.035 1];
a = f<wc;

W = [1,1];

h = firgr(30, f, a, W);
Hd = w<wc | w>2-wc;
H = fft(h,N);

plot(w,abs(H));
hold on
plot(w,Hd,'r');
H = abs(fft(h,length(w)));

