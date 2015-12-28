clear;clc;close all; warning off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Part II: IIR Architecture  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
%%%         elliptic            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elliptic_order    =   ellipord(pw_omega_p, pw_omega_s, Ap, Aa,'s');
[NUM_nobz, DEN_nobz] = ellip(elliptic_order, Ap, Aa, pw_omega_p, 's');


[NUM, DEN] = bilinear(NUM_nobz, DEN_nobz, fs);
[Z_ellip, P_ellip, K_ellip] = tf2zp(NUM, DEN);

% K = NUM(1);
% NUM = NUM/ K;
% subplot(2,2,4)
H_ellip = freqz(NUM,DEN,W);
plot(W*fs/2/pi,abs(H_ellip));
% xlim([5000 7000]);
% ylim([0 1.1]);
% title('Elliptic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Part 4           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4.1  	Due to the long duration of the impulse response of a narrowband
% IIR filter (see Figure 5), the response time of your system is considered
% important.  In the absence of any additive background interference, how
% long (in samples) after the tone arrives would it take for the output to
% raise to 50%, 75% and 95% of its steady state value (i.e., signaling
% latency)


t = 0: 1/fs: 30/f-1/fs;
u_tone = cos(2*pi*f*t);
% u_tone = zeros(1,length(t));
% u_tone(1)=1;


y = filter(NUM, DEN, u_tone);
subplot 211
plot(u_tone);
hold on
plot(y,'-ro');
hold off
H_filter = freqz(NUM,DEN,W);
% plot(W*fs/2/pi,abs(H_ellip));
% figure
subplot 212
plot(W_Hz, grpdelay(NUM, DEN, W));
% figure 
% plot(y)

