clear;clc;close all;

%%% 8th order chebyshev II , l_inf norm

k = 0.0888;
DEN = [1 1.10  1.97  1.55  1.22   .61   .24   .061 0.008 ];
NUM = [1 4.43 10.76 17.46 20.48 17.46 10.76  4.43  1     ];


[A, b, c, d] = tf2ss(NUM,DEN);

c = [0 0 0 1 0 0 0 0];

[NUM_i, DEN_i] = ss2tf(A,b,c,d);

H = freqz(NUM_i, DEN_i,1024*32);

% max(abs(H))
% 
% plot(abs(H))


x0 = k*[1 0 0 0 0 0 0 0]';
x = zeros(8, 1024*32);
x(:,1) = x0;
for i = 2:1024*32
    x(:,i) = A*x(:,i-1);
end

h = x(1,:);

H = fft(h);
plot(abs(H));
max(abs(H))