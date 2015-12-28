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






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Part 2.1            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 2.1 	Implement, using MATLAB, your elliptic IIR as a Direct II and
% Cascade filter using MATLAB-enabled architectures. What are the resulting
% 4-tuples [A,b,c,d] per architecture, along with any associated non-unity
% input scale factor k? Illustration:  Direct II/Cascade Solutions Given an
% % IIR with the transfer function H(z)=k Num(z)/Den(z), where:
% scale factor  = 0.0027;
% num = 1.0000   -3.0136    3.2621   -0.0000   -3.2621    3.0136  -1.000;
% den = 1.0000   -4.4962    9.6788  -12.1817    9.4904   -4.3228 0.9427;
% the Direct II filter 4-tuple is (using normal state  representation)
%        A = [  0    1.0000         0         0         0         0
%          0         0    1.0000         0         0         0 0         0
%          0    1.0000         0         0 0         0         0         0
%          1.0000         0 0         0         0         0         0
%          1.0000
%    -0.9427    4.3228   -9.4904   12.1817   -9.6788    4.4962];
% b = [0;0;0;0;0;1]; 	c = [-1.9427  7.3364  -12.7524  12.1817  -6.4167
% 1.4825];	d = 1; For Cascade, the 2nd sections are individually:
% scale  factor = 0.0027;
% SOS =[1.0000   -0.0000   -1.0000    1.0000   -1.4902  0.9682
%      1.0000   -1.4151    1.0000    1.0000   -1.4660  0.9863
%      1.0000   -1.5986    1.0000    1.0000   -1.5400  0.9872];
%   For the each sos, the state 4-tuples are:
%   A1 = [0 1;-0.9682 1.4902]; b1 = [0;1]; c1 =
% [-1.9682 1.4902]; d1 = 1; A2 = [0 1;-0.9863 1.4660]; b2 = [0;1]; c2 =
% [0.0137 0.0509]; d2 = 1; A3 = [0 1;-0.9872 1.5400]; b3 = [0;1]; c3 =
% [0.0128 -0.0586]; d3 = 1;


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

[A_d2, B_d2, C_d2, D_d2] = tf2ss(NUM, DEN);

[SOS_c] = tf2sos(NUM, DEN);

for l = 1: length(SOS_c(:,1))
    [A_c(:,:,l),B_c(:,:,l), C_c(:,:,l), D_c(:,:,l)] = tf2ss(SOS_c(l,1:3), SOS_c(l, 4:6));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% part 2.2            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.2	Summarize your designs in table form (see Table 3) which lists the
% absolute value of the largest and smallest (non-zero) Direct II and
% Cascade filter coefficients or constant, and the number of non-zero,
% non-unity coefficients associated with each of the architectures. Table
% 3: Filter Coefficients Summary Architecture\coefficient	|Largest|
% |Smallest|	#of nonzero, non-unity coefficients Direct II Cascade


%%% Direct II
Max_coeff_d2 = max([max(max(A_d2)) max(D_d2)  max(C_d2) max(D_d2)]);
Min_coeff_d2 = min([min(min(abs(A_d2(find(A_d2 ~= 0))))) min(abs(C_d2(find(C_d2 ~= 0))))  ...
    min(abs(B_d2(find(B_d2 ~= 0)))) min(abs(D_d2(find(D_d2 ~= 0))))]);
zero_unity_d2 = length([find(A_d2 == 0 )' find(B_d2 == 0 )' ...
    find(C_d2 == 0 ) find(D_d2 == 0 )' ]) + ...
    length([find(A_d2 == 1 )' find(B_d2 == 1 )' ...
    find(C_d2 == 1 ) find(D_d2 == 1 )' ]);
Non_zero_unity_d2 = numel(A_d2) + numel(B_d2) + numel(C_d2) + numel(D_d2) - zero_unity_d2;

%%% Cascade
Max_coeff_c = max([max(max(max(A_c))) max(max(B_c))  max(max(C_c)) max(D_c)]);
Min_coeff_c = min([min(min(min(abs(A_c(find(A_c ~= 0)))))) min(min(abs(C_c(find(C_c ~= 0)))))  ...
    min(min(abs(B_c(find(B_c ~= 0))))) min(abs(D_c(find(D_c ~= 0))))]);
zero_unity_c = length([find(A_c == 0 )' find(B_c == 0 )' ...
    find(C_c == 0 ) find(D_c == 0 )' ]) + ...
    length([find(A_c == 1 )' find(B_c == 1 )' ...
    find(C_c == 1 ) find(D_c == 1 )' ]);
Non_zero_unity_c = numel(A_c) + numel(B_c) + numel(C_c) + numel(D_c) - zero_unity_c;


%%% Table format

fprintf('Architecture/Coeff\t Largest\t Smallest\t #non-zero-unity \n');
fprintf('%s\t %d\t %d\t %d\n','Direct II', Max_coeff_d2, Min_coeff_d2, Non_zero_unity_d2)
fprintf('%s\t %d\t %d\t %d\n','Cascade  ', Max_coeff_c, Min_coeff_c, Non_zero_unity_c)

%%%%%%%%%%%%%%%%
%%% part 2.3 %%%
%%%%%%%%%%%%%%%%

% 2.3 Since the audio application requires a relatively low sample rate
% (44.1 kSa/s) and cost is not an issue, you consider a floating-point
% design. Based on your architectural analysis, which architecture would
% you recommend and why?


fprintf('\n\npart 2.3: I  think it is not an issue. \nMax integer for direct II ~ 10 and for cascade ~ 2.\ndirect II needs more integer bits\n\n\n ')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Part III: IIR Test and Evaluation %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider again the Direct II and Cascade architectures studied in Part II.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 3.1	Determine the L2 state norm for each state and output for each %%%
%%% architecture.                                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% calculate the lyapunov matrix

%%% Direct II
B0_d2 = B_d2*B_d2';
K0_d2 = B0_d2;

for i =1 :1000
    K0_d2 = A_d2*K0_d2*A_d2' + B0_d2;
end

%%% Cascade

for section = 1 : length(D_c)
    B0_c(:,:,section) = B_c(:,:,section)*B_c(:,:,section)';
    K0_c(:,:,section) = B0_c(:,:,section);
    for i = 1:1000
        K0_c(:,:,section) = A_c(:,:,section)*K0_c(:,:,section)*A_c(:,:,section)' + B0_c(:,:,section);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 3.2  	Determine the L_inf frequency domain state norm for each state
%%% and output for each architecture.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Direct II
N_time = N_freq;
x = zeros(length(A_d2), N_time);
x0 = zeros(length(A_d2), 1); x0(1) = 1;
x(:,1) = x0;

for i = 2:N_time
    x(:,i) = A_d2*x(:,i-1);
end

h_d2 = x(1,:);
H_d2 = fft(h_d2, N_freq*2);

figure
subplot(2,2,1)
plot(W_Hz, abs(H_d2(1:length(H_d2)/2)))
title('L-inf state norm for Direct II');
Linf_d2 = max(abs(H_d2));



%%% Cascade

x_c = zeros(2, N_time,length(D_c)); % second order sections, length(D_c) in number
x0 = [1 0]';
h_c = zeros(length(D_c), N_time);
Linf_c = zeros(1,length(D_c));

for section = 1:length(D_c)
    x_c(:,1,section) = x0;
    for i = 2:N_time
        x_c(:,i,section) = A_c(:,:,section)*x_c(:,i-1,section);
    end

    h_c(section, :) = x_c(1,:,section);
    H_c = fft(h_c(section,:), N_freq*2);
    subplot(2,2,section+1)
        plot(W_Hz, abs(H_c(1:length(H_c)/2)))
        title('L-inf state norm for Cascade section %d');
    Linf_c(section) = max(abs(H_c));
    hold on

end

% 3.3  	Determine the L1 state norm for each state and output for each
% architecture.


% Direct II
L1_d2 = sum(abs(h_d2));

% Cascade
L1_c = zeros(1,length(D_c));
for section = 1:length(D_c)
    L1_c(section)=sum(abs(h_c(section,:)));
end

%%% Table format
% Table the results from Parts 3.1, 3.2, and 3.3.
% Table 4: Lp norms
% L2	x1	x2	x3	x4	x5	x6	y
% Direct II
% Cascade
% Linf	x1	x2	x3	x4	x5	x6	y
% Direct II
% Cascade
% L1	x1	x2	x3	x4	x5	x6	y
% Direct II
% Cascade
%


fprintf('Direct II:\t L2\t\t\t\t Linf\t\t\t\t L1\t\n');
fprintf('\t\t%d\t\t%d\t\t%d\t\n',sqrt(K0_d2(1,1)), Linf_d2, L1_d2);

fprintf('Cascade:\t L2\t\t\t\t Linf\t\t\t\t L1\t\n');
for i = 1:length(D_c)
    fprintf('\t\t%d\t\t%d\t\t%d\t\n',sqrt(K0_c(1,1,i)), Linf_c(i), L1_c(i));
end

close all

% % 3.4	 Conduct a finite wordlength study for a signed 16-bit word
% % architecture consisting of I integer bits and F fractional bits
% % (I+F+1=16) for both the Direct II and Cascade architecture.  Using
% % MATLAB, design and implement a fixed-point simulator (without extended
% % precision accumulators (i.e., M=0)), for each architecture. The simulator
% % will need to implement fixed-point arithmetic for the production all
% % states and outputs involving a non-unity or non-zero production rules.
% % Add to you simulation a counter that will determine the number of
% % system-wide accumulator overflows that occur run-time. Handle all
% % run-time overflows using saturating arithmetic. Use floating-point to
% % define an assumed infinite precision (ideal) response, and define the
% % error to be the difference between the ideal and fixed-point response.
% % Conduct a study for F ranging from 15 to 0-bits, using a common unit
% % bound signed random noise record of length at least 1024 samples,
% % preferably longer.  Display the outcome in graphical form, plotting
% % realized fractional precision against F (in bits) as well as the number
% % of run-time overflows.  A motivational example is shown in Figure 4.
% 
% 
% N_time = 100;
% N_B = 16; I = 0; S = 1; F =15;
% 
% x_d2    = zeros(length(A_d2),N_time);
% x_d2_pr = zeros(length(A_d2),N_time);
% x_d2_r = zeros(length(A_d2),N_time);
% 
% % u = rand(1,N_time);
% u = zeros(1,N_time);
% u(1) = 1;
% for F = 1:15
%     A_d2_r = fi(A_d2, S, N_B, F);
%     B_d2_r = fi(B_d2, S, N_B, F);
%     C_d2_r = fi(C_d2, S, N_B, F);
%     D_d2_r = fi(D_d2, S, N_B, F);
% 
%     overflow_count_x_d2 = 0;
%     overflow_count_y_d2 = 0;
% 
% 
%     for i = 2:N_time
%         x_d2(:,i) = A_d2(:,:)*x_d2(:,i-1) + B_d2*u(i);
%         x_d2_pr(:,i) = A_d2_r(:,:)*x_d2_r(:,i-1) + B_d2_r*u(i);
% 
%         if(x_d2_pr(1,i) > 1 || x_d2_pr(1,i) < -1)
%             overflow_count_x_d2 = overflow_count_x_d2 + 1;
%         end
% 
%         x_d2_r(:,i) = fi(x_d2_pr(:,i), S, N_B, F);
%         e_r_d2(i) = x_d2_r(1,i) - x_d2(1,i);
% 
%         y_d2(i) = C_d2*x_d2(:,i) + D_d2*u(i);
%         y_d2_pr(i) = C_d2_r*x_d2_r(:,i) + D_d2_r*u(i);
%         y_d2_r(i) = fi(y_d2_pr(i), S, N_B, F);
% 
% 
%         if(y_d2_pr(i)>1 || y_d2_pr(i)<-1)
%             overflow_count_y_d2 = overflow_count_y_d2 + 1;
%         end
% 
%     end
%     var_e_d2(F) = var(e_r_d2);
%     %
%     % overflow_count_y_d2
%     % overflow_count_x_d2
% 
%     % fprintf('direct II var(e_r_d2) =%d \tlog2(var(e_r_d2))= %d\n',var(e_r_d2),log2(var(e_r_d2)))
% 
%     %%%cascade
% 
% 
% 
%     x_c    = zeros(2,N_time,length(D_c));
%     x_c_pr = zeros(2,N_time,length(D_c));
%     x_c_r = zeros(2,N_time,length(D_c));
% 
% 
%     y_c    = zeros(1,N_time,length(D_c));
%     y_c_pr = zeros(1,N_time,length(D_c));
%     y_c_r = zeros(1,N_time,length(D_c));
% 
%     u = ones(1,N_time);
% 
%     A_c_r = fi(A_c, S, N_B, F);
%     B_c_r = fi(B_c, S, N_B, F);
%     C_c_r = fi(C_c, S, N_B, F);
%     D_c_r = fi(D_c, S, N_B, F);
% 
%     % A_c_r = A_c;
%     % B_c_r = B_c;
%     % C_c_r = C_c;
%     % D_c_r = D_c;
% 
%     overflow_count_x_c = 0;
%     overflow_count_y_c = 0;
% 
%     for sect = 1:length(D_c)
%         for i = 2:N_time
% 
% 
%             if(sect == 1)
%                 x_in(i) = u(i);
%                 x_in_r(i) = u(i);
%             else
%                 x_in(i) = y_c(1,i,sect-1);
%                 x_in_r(i) = y_c_r(1,i,sect-1);
%             end
%             x_c(:,i,sect) = A_c(:,:,sect)*x_c(:,i-1,sect) + B_c(:,1,1)*x_in(i);
%             x_c_pr(:,i,sect) = A_c_r(:,:,sect)*x_c_r(:,i-1,sect) + B_c_r(:,1,1)*x_in_r(i);
% 
% 
%             if(x_c_pr(1,i) > 1 || x_c_pr(1,i) < -1)
%                 overflow_count_x_c = overflow_count_x_c + 1;
%             end
% 
%             x_c_r(:,i,sect) = fi(x_c_pr(:,i,sect), S, N_B, F);
%             e_r_c(i,sect) = x_c_r(1,i,sect) - x_c(1,i,sect);
% 
%             y_c(1,i,sect) = C_c(1,:,sect)*x_c(:,i,sect) + D_c(sect)*x_in(i);
%             y_c_pr(1,i,sect) = C_c_r(1,:,sect)*x_c_r(:,i,sect) + D_c_r(sect)*x_in_r(i);
%             y_c_r(i) = fi(y_c_pr(i), S, N_B, F);
% 
% 
%             if(y_c_pr(i)>1 || y_c_pr(i)<-1)
%                 overflow_count_y_c = overflow_count_y_c + 1;
%             end
%         end
% 
%     end
%     %
%     % overflow_count_y_c
%     % overflow_count_x_c
%     %     fprintf('cascade var(e_r_c) =%d \tlog2(var(e_r_c))= %d\n',var(e_r_c),log2(var(e_r_c)))
%     var_e_cascade(:,F) = var(e_r_c);
% end
% plot([1:15],log2(var_e_cascade(1,:))/2);% I'm cheating by dividing by two. Something what engineers do - according to Dr Taylor
% ylabel('Precision bits');
% xlabel('Fractional bits');
% title('Plot of fractional bits vs Precision bits')
% % hold on
% % plot([1:15],log2(var_e_cascade(2,:))/2,'r');
% % plot([1:15],log2(var_e_cascade(3,:))/2);
% log2(L1_c)
% hold on
% % plot([1:15],log2(var_e_d2),'r');
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%         part 3.5             %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Best 16-bit solution is Cascade with 13 bits of Fractional precision
% 
% t = 0: 1/fs: 100/f-1/fs;
% x = cos(2*pi*t*f);
% % plot(x)
% F=13;
% energy_in_x = sum(x.^2)/length(x);
% iter = 1;
% for SNR_dB = [1/0.000001 10 1 -10]
%     snr = 10^(SNR_dB/10);
%     
%     noise_power = energy_in_x/snr;
%     
%     noise = randn(length(x), 1)*sqrt(noise_power);
%     
%     x_noise = x + noise';
%     
%     y_noise = filter(NUM, DEN, x_noise);
%     y_wnoise = filter(NUM, DEN, x);
%     y_r = fi(y_noise, S, N_B, F);
%     
%     filter_noise = y_wnoise - y_noise;
%     
%     SNR_filtered(iter) = var(y_wnoise)/var(filter_noise);
%     
%     
% %     subplot(2,2,iter)
% %     plot(abs(fft(x_noise)));
% %     hold on;
% %     plot(abs(fft(y)),'r');
%     
%     iter = iter+1;
%     
% end
% 
% SNR_dB
% 10*log10(SNR_filtered)
% figure
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% part 3.7           %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% for SNR_dB = [1/0.000001 ]
%     snr = 10^(SNR_dB/10);
%     
%     noise_power = energy_in_x/snr;
%     
%     noise = randn(length(x), 1)*sqrt(noise_power);
%     
%     x_noise = x + noise';
%     
%     y_noise = filter(NUM, DEN, x_noise);
%     y_wnoise = filter(NUM, DEN, x);
%     y_r = fi(y_noise, S, N_B, F);
%     
%     error = y_r - y_noise;
%     filter_noise = y_wnoise - y_noise;
%     
%     SNR_filtered(iter) = var(y_wnoise)/var(filter_noise);
%     
%     
% %     subplot(2,2,iter)
% %     plot(abs(fft(x_noise)));
% %     hold on;
% %     plot(abs(fft(y)),'r');
%     
%     iter = iter+1;
%     
% end
% 
% plot(x(1:200))
% hold on
% plot(y_r(1:200),'r')
% 
% fprintf('error mean %d\t error variance=%d\n',mean(y_r),var(y_r))
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Part 4           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4.1  	Due to the long duration of the impulse response of a narrowband
% IIR filter (see Figure 5), the response time of your system is considered
% important.  In the absence of any additive background interference, how
% long (in samples) after the tone arrives would it take for the output to
% raise to 50%, 75% and 95% of its steady state value (i.e., signaling
% latency)


t = 0: 1/fs: 20/f-1/fs;
u_tone = cos(2*pi*f*t);


y = filter(NUM, DEN, u_tone);
plot(u_tone);
hold on
plot(y,'-ro');
hold off
H_filter = freqz(NUM,DEN,W);
plot(W*fs/2/pi,abs(H_ellip));
plot(W_Hz, grpdelay(NUM, DEN, W));
