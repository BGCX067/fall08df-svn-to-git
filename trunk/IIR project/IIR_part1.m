clear;clc;close all;
warning off

% center frequency and sampling frequency
f   =   6*1e3;    % in Hz
fs  =   44.1*1e3; % in Hz

%cutoff frequencies, in Hz
fp  =   [f-0.2*1e3 f+0.2*1e3]
fa  =   [f-0.5*1e3 f+0.5*1e3]

%cutoff frequencies in digital domain (0-2pi). 
wp = 2*pi*fp/fs;
wa = 2*pi*fa/fs;

%pre-warping, rad/s. 
pw_omega_p = 2*fs*tan(wp/2);
pw_omega_s = 2*fs*tan(wa/2);

%allowed error deviation in dB
Ap  =   1.5; 
Aa  =   40;

% define frequency points
N_freq = 1024*32; %number of points in freq domain
W = 0:pi/N_freq:pi - pi/N_freq;
W_Hz = W*fs/2/pi;
W_pass = wa(1): pi/N_freq:wa(2)-pi/N_freq;  %includes transistion band


% define passband and stop-band masks
pass_mask = (wp(1) < W) & (wp(2) > W);

stop_mask = (W < wa(1)) | (W > wa(2));
lower_stop_mask = W<wa(1);
upper_stop_mask = W>wa(2);

H_desired = pass_mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  part 1.1                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        butterworth          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% design an analog filter matching the pre-warped frequencies
butter_order    =   buttord(pw_omega_p, pw_omega_s, Ap, Aa,'s');
[B, A] = butter(butter_order, pw_omega_p,'s'); % s for analog


% substitute s=2*fs*(z-1)/(z+1). warping should take place here
[Bz, Az] = bilinear(B, A, fs);
[Z_butter, P_butter, K] = tf2zp(Bz, Az);



subplot(2,2,1)
H = freqz(Bz,Az,W);
H_butter = H;
plot(W*fs/2/pi,abs(H));

xlim([f-1e3 f+1e3]);
ylim([0 1.1]);

title('Butterworth');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        chebyshev I          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cheby1_order    =   cheb1ord(pw_omega_p, pw_omega_s, Ap, Aa,'s');
[B, A] = cheby1(cheby1_order, Ap, pw_omega_p, 's');

% substitute s=2*fs*(z-1)/(z+1). warping should take place here
[Bz, Az] = bilinear(B, A, fs);
[Z_cheby1, P_cheby1, K] = tf2zp(Bz, Az);
subplot(2,2,2)
H = freqz(Bz,Az,W);
H_cheby1 = H;
plot(W*fs/2/pi,abs(H));



xlim([f-1e3 f+1e3]);
ylim([0 1.1]);

title('Chebyshev I');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       chebyshev II          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cheby2_order    =   cheb2ord(pw_omega_p, pw_omega_s, Ap, Aa,'s');
[B, A] = cheby2(cheby2_order, Aa, pw_omega_s, 's');


[Bz, Az] = bilinear(B, A, fs);
[Z_cheby2, P_cheby2, K] = tf2zp(Bz, Az);

subplot(2,2,3)
H = freqz(Bz,Az,W);
H_cheby2 = H;

plot(W*fs/2/pi,abs(H));
xlim([f-1e3 f+1e3]);
ylim([0 1.1]);
title('Chebyshev II');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         elliptic            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elliptic_order    =   ellipord(pw_omega_p, pw_omega_s, Ap, Aa,'s');
[B_ellip, A_ellip] = ellip(elliptic_order, Ap, Aa, pw_omega_p, 's');

[Bz_ellip, Az_ellip] = bilinear(B_ellip, A_ellip, fs);
[Z_ellip, P_ellip, K_ellip] = tf2zp(Bz_ellip, Az_ellip);


subplot(2,2,4)
H_ellip = freqz(Bz_ellip,Az_ellip,W);
plot(W*fs/2/pi,abs(H_ellip));
xlim([f-1e3 f+1e3]);
ylim([0 1.1]);
title('Elliptic');




%%%%%%%%%%%%%%%%%%%%%%%%%
%     1.1 error         %
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%error calculation
error_butter = abs(abs(H_butter) - H_desired);
max_passband_error_dB_butter = 20*log10(max(pass_mask.*error_butter));
max_lower_stopband_error_dB_butter = 20*log10(max(lower_stop_mask.*error_butter));
max_upper_stopband_error_dB_butter = 20*log10(max(upper_stop_mask.*error_butter));

%%%error calculation
error_chebyI = abs(abs(H_cheby1) - H_desired);
max_passband_error_dB_chebyI = 20*log10(max(pass_mask.*error_chebyI));
max_lower_stopband_error_dB_chebyI = 20*log10(max(lower_stop_mask.*error_chebyI));
max_upper_stopband_error_dB_chebyI = 20*log10(max(upper_stop_mask.*error_chebyI));

%%%error calculation
error_chebyII = abs(abs(H_cheby2) - H_desired);
max_passband_error_dB_chebyII = 20*log10(max(pass_mask.*error_chebyII));
max_lower_stopband_error_dB_chebyII = 20*log10(max(lower_stop_mask.*error_chebyII));
max_upper_stopband_error_dB_chebyII = 20*log10(max(upper_stop_mask.*error_chebyII));

%%%error calculation
error_elliptic = abs(abs(H_ellip) - H_desired);
max_passband_error_dB_elliptic = 20*log10(max(pass_mask.*error_elliptic));
max_lower_stopband_error_dB_elliptic = 20*log10(max(lower_stop_mask.*error_elliptic));
max_upper_stopband_error_dB_elliptic = 20*log10(max(upper_stop_mask.*error_elliptic));


fprintf('Filter\t Order\t Max. P Err\t Max. LS Err\t Max. US Err\n')
fprintf('%s\t %d\t %d\t %d\t %d\n','butter', butter_order, max_passband_error_dB_butter, max_lower_stopband_error_dB_butter, max_upper_stopband_error_dB_butter)
fprintf('%s\t %d\t %d\t %d\t %d\n','cheby1', cheby1_order, max_passband_error_dB_chebyI, max_lower_stopband_error_dB_chebyI, max_upper_stopband_error_dB_chebyI)
fprintf('%s\t %d\t %d\t %d\t %d\n','cheby2', cheby2_order, max_passband_error_dB_chebyII, max_lower_stopband_error_dB_chebyII, max_upper_stopband_error_dB_chebyII)
fprintf('%s\t %d\t %d\t %d\t %d\n','elliptic', elliptic_order, max_passband_error_dB_elliptic, max_lower_stopband_error_dB_elliptic, max_upper_stopband_error_dB_elliptic)




%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  1.1  Pole Zero   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(1,2,1)
zplane(Z_butter,P_butter);
hold on; 
zplane(Z_cheby1, P_cheby1,'r');
zplane(Z_cheby2, P_cheby2,'k');
zplane(Z_ellip, P_ellip, 'g');
axis([-1.5 1.5 -1.5 1.5]);

subplot(1,2,2)
zplane(Z_butter,P_butter);
hold on; 
zplane(Z_cheby1, P_cheby1,'r');
zplane(Z_cheby2, P_cheby2,'k');
zplane(Z_ellip, P_ellip, 'g');
axis([0.58 0.72 0.65 0.82]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 tf, mag resp, group delay  
%1.2  	Choose an elliptic IIR for continued study. If the minimum order of
%the elliptic designed in Part 1.1 exceeds order N=6, then use (re-design)
%a 6th order elliptic IIR that best meets the specification (i.e., relax
%the stopband attenuation). Display the 6th order filter’s transfer
%function, magnitude frequency response (linear and dB), and group delay.%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defined earlier;
[Bz_ellip, Az_ellip];
elliptic_order;
[Z_ellip, P_ellip];

H = freqz(Bz_ellip,Az_ellip,W);
figure;
title('Elliptic filter, order 6');

subplot 221 
plot(W_Hz, abs(H));
title('Magnitude frequency response, linear');
subplot 223
plot(W_Hz, 20*log10(abs(H)));
title('Magnitude frequency response, dB');

subplot 222
plot(W_Hz, angle(H));
title('Phase response');
subplot 224
plot(W_Hz, grpdelay(Bz_ellip, Az_ellip, W));
title('group delay');



%%%%%%%%%%%%%%%%%%%%%%%%%
%  1.3 pre-emphasize    %
%%%%%%%%%%%%%%%%%%%%%%%%%
% You are motivated to design a pre-emphasized IIR having a slightly
% increasing the gain in the upper range of the passband relative to the
% lower range of the passband by 15%.  Achieve this design goal by
% manipulating the poles and/or zeros of your 6th order filter (see Figure
% 2 for an example of increasing the lower range of the passband by 15%).

P_ellip_emph = P_ellip.*[1.00065 1.00065 1 1  1 1]';
[B_emph, A_emph] = zp2tf(Z_ellip, P_ellip_emph, K_ellip);

figure
subplot 221
zplane(Z_ellip, P_ellip);
axis([0.58 0.72 0.65 0.82]);

subplot 222
zplane(Z_ellip, P_ellip_emph);
axis([0.58 0.72 0.65 0.82]);


pass_index = find(pass_mask>0);
H = freqz(Bz_ellip,Az_ellip,W);

subplot 223
plot(W_Hz(pass_index), abs(H(pass_index)));


subplot 224
H_emph = freqz(B_emph, A_emph, W);
plot(W_Hz(pass_index), abs(H_emph(pass_index)));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4 test with tone and BB noise %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1.4  	Using MATLAB and the 6th order
% IIR from Part 1.2, simulate the filter’s performance in the presence of a
% tone and broadband noise.  The broadband audio process is to be modeled
% as a long zero mean Gaussian random time.  Adjust the background noise to
% have SNR = {INFTY, 10, 0, -10}, compute and display the input and output
% power spectral densities (overlay).  What is the SNR measured at the
% filter’s output?  Does the outcome agree with your intuitive and
% mathematical understanding on how the filter should work? A simple
% example of power spectrum analysis, shown below, was conducted using
% MATLAB’s pmtm function.  You may or may not want to go down this path.
% You can use any technique that is mathematically consistent with the
% studies goals.  
% 
% >> Fs=1000; t=0:1/Fs:.5; 
% >> x=cos(2*pi*t*100)+randn(size(t)); 
% >> [Pxx,Pxxc,F]=pmtm(x,4,[],Fs,'adapt'); 
% >> plot(Pxx)

figure;
t = 0: 1/fs: 100/f-1/fs;
x = cos(2*pi*t*f);
plot(x)

energy_in_x = sum(x.^2)/length(x);
iter = 1;
for SNR_dB = [1/0.000001 10 1 -10]
    snr = 10^(SNR_dB/10);
    
    noise_power = energy_in_x/snr;
    
    noise = randn(length(x), 1)*sqrt(noise_power);
    
    x_noise = x + noise';
    
    y = filter(Bz_ellip, Az_ellip, x_noise);
    
    
    subplot(2,2,iter)
    plot(abs(fft(x_noise)));
    hold on;
    plot(abs(fft(y)),'r');
    
    iter = iter+1;
    
end