clc;close all;clear
x = ones(11,1);

blk_win = blackman(length(x));

x_win = x.*blk_win;

X = fft(x,31);

X_WIN = fft(x_win, 31);

max(X)
max(X_WIN)

subplot(2,2,1);
stem(x);
subplot(2,2,2);
plot(abs(X));

subplot(2,2,3);
stem(x_win);
subplot(2,2,4);
plot(abs(X_WIN));

