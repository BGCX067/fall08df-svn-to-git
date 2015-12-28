clear;clc;close all;

%%% 8th order chebyshev II , run_time overflow

k = 0.0888;
NUM = [1 4.43 10.76 17.46 20.48 17.46 10.76  4.43  1     ];
DEN = [1 1.10  1.97  1.55  1.22   .61   .24   .061 0.008 ];


[A, b, c, d] = tf2ss(NUM*k,DEN);

% c = [0 0 0 1 0 0 0 0];

% [NUM_i, DEN_i] = ss2tf(A,b,c,d);
%
% H = freqz(NUM_i, DEN_i,1024*32);

% max(abs(H))
%
% plot(abs(H))

N_time = 1000;
% u = ones(1, N_time);
u = randn(1,N_time);
x_true = zeros(8, N_time);
% y = zeros(1, N_time);

x_r = zeros(8,N_time);
x_pr = zeros(8,N_time);

N_B = 16; S = 1;I = 0; F = 15;%N_B = I+F+S;

x0 = b*u(1);
% x_true(:,1) = x0;
% x_r(:,1) = x0;
A_r = fi(A, S, N_B, F);
b_r = fi(b, S, N_B, F);
c_r = fi(c, S, N_B, F);
d_r = fi(d, S, N_B, F);

overflow_count_x = 0;
overflow_count_y = 0;
for i = 2:N_time
    %
    %     for m = 1:8
    %         x_true(m,i) = b(m)*u(i);
    %         x_pr(m,i) = b(m)*u(i);
    %         for n =1:8
    %             x_true(m,i) = x_true(m,i) + A(m,n)*x_true(n,i-1);
    %             temp = A_r(m,n)*x_r(n,i-1);
    %             x_pr(m,i) = x_pr(m,i) + temp;
    %
    %
    % %
    % %             if(x_pr(m,i) > 1 || x_pr(m,i) < -1)
    % %                 overflow_count_x = overflow_count_x + 1;
    % %             end
    %         end
    %     end

    x_true(:,i) = A(:,:)*x_true(:,i-1) + b(:)*u(i);
    x_pr(:,i) = A_r(:,:)*x_r(:,i-1) + b_r(:)*u(i);

    if(x_pr(1,i) > 1 || x_pr(1,i) < -1)
        overflow_count_x = overflow_count_x + 1;
    end

    x_r(:,i) = fi(x_pr(:,i), S, N_B, F);
    e_r(i) = x_r(1,i) - x_true(1,i);
    %     x_r;
    y_true(i) = c*x_true(:,i) + d*u(i);
    y_r(i) = fi(c_r*x_r(:,i) + d_r*u(i), S, N_B, F);
    y_pr(i) = c_r*x_r(:,i) + d_r*u(i);

    if(y_pr(i)>1 || y_pr(i)<-1)
        overflow_count_y = overflow_count_y + 1;
    end
end

% plot(e_r)
hold on
% plot(y_r(1,:),'r')
% plot(y_pr(1,:))
% plot(y_true(1,:),'k')
overflow_count_x
overflow_count_y

% plot(y_true)
% hold on
% plot(y_r,'r')
% plot(y_pr, 'g')


%
% %%%% cascade
SOS_c = tf2sos(NUM*k,DEN,'up','inf');
for l = 1: length(SOS_c(:,1))
    [A_c(:,:,l),B_c(:,:,l), C_c(:,:,l), D_c(:,:,l)] = tf2ss(SOS_c(l,1:3), SOS_c(l, 4:6));
end

% u = ones(1, N_time);
xc_true = zeros(2, N_time, length(D_c));
xc_r = zeros(2, N_time, length(D_c));

yc_true  = zeros(1, N_time, length(D_c));
yc_r  = zeros(1, N_time, length(D_c));

Ac_r = fi(A_c, S, N_B, F);
Bc_r = fi(B_c, S, N_B, F);
Cc_r = fi(C_c, S, N_B, F);
Dc_r = fi(D_c, S, N_B, F);
overflow_count_xc = 0;
overflow_count_yc = 0;
for sect = 1:length(D_c)
    for i = 2:N_time
        if(sect == 1)
            xc_true(:,i,sect) = A_c(:,:,sect)*xc_true(:,i-1,sect) + B_c(:,1,sect)*u(i);
            xc_r(:,i,sect) = fi(Ac_r(:,:,sect)*xc_r(:,i-1,sect) + Bc_r(:,1,sect)*u(i), S, N_B, F);
            xc_pr(:,i,sect) = Ac_r(:,:,sect)*xc_r(:,i-1,sect) + Bc_r(:,1,sect)*u(i);
            if(xc_pr(1,i) > 1 || xc_pr(1,i) < -1)
                overflow_count_xc = overflow_count_xc + 1;
            end


            yc_true(1,i,sect)  = C_c(1,:,sect)*xc_true(:,i  ,sect) + D_c(sect)*u(i);
            temp_true(:,i,sect) = C_c(1,:,sect)*xc_true(:,i  ,sect);
            yc_r(1,i,sect)  = fi(Cc_r(1,:,sect)*xc_r(:,i  ,sect) + Dc_r(sect)*u(i), S, N_B, F);
            yc_pr(1,i,sect)  = Cc_r(1,:,sect)*xc_r(:,i  ,sect) + Dc_r(sect)*u(i);
%             if(yc_pr(1,i,sect)>1 || yc_pr(1,i,sect)<-1)
%                 overflow_count_yc = overflow_count_yc + 1;
%             end

        else
            xc_true(:,i,sect) = A_c(:,:,sect)*xc_true(:,i-1,sect)  + B_c(:,1,sect)*yc_true(:,i,sect-1);
            xc_r(:,i,sect) = fi(Ac_r(:,:,sect)*xc_r(:,i-1,sect) + Bc_r(:,1,sect)*yc_r(:,i,sect-1), S, N_B, F);
            xc_pr(:,i,sect) = Ac_r(:,:,sect)*xc_r(:,i-1,sect) + Bc_r(:,1,sect)*yc_r(:,i,sect-1);
            if(xc_pr(1,i) > 1 || xc_pr(1,i) < -1)
                overflow_count_xc = overflow_count_xc + 1;
            end

            yc_true(1,i,sect)  = C_c(1,:,sect)*xc_true(:,i  ,sect) + D_c(sect)*yc_true(1,i,sect-1);
            temp_true(:,i,sect) = C_c(1,:,sect)*xc_true(:,i-1  ,sect);

            yc_r(1,i,sect)  = fi(Cc_r(1,:,sect)*xc_r(:,i  ,sect) + Dc_r(sect)*yc_r(1,i,sect-1), S, N_B, F);
            yc_pr(1,i,sect)  = Cc_r(1,:,sect)*xc_r(:,i  ,sect) + Dc_r(sect)*yc_r(1,i,sect-1);

            if((yc_pr(1,i,sect)>1 || yc_pr(1,i,sect)<-1) )
                overflow_count_yc = overflow_count_yc + 1;
            end


        end
    end
end
%
% close all;
% subplot(2,2,1)
% plot(xc_true(1,:,1))

y_1 = filter(SOS_c(1,1:3), SOS_c(1,4:6),u  );
y_2 = filter(SOS_c(2,1:3), SOS_c(2,4:6),y_1);
y_3 = filter(SOS_c(3,1:3), SOS_c(3,4:6),y_2);
y_4 = filter(SOS_c(4,1:3), SOS_c(4,4:6),y_3);

% plot(y_4)
%
% hold on
% plot(temp_true(1,:,4),'g')
% plot(u(1,:,1)*D_c(1) + temp_true(1,:,1),'k')
% plot(yc_true(1,:,1),'r')
%
% subplot(2,2,2)
% plot(xc_true(1,:,2))
% hold on
% plot(yc_true(1,:,2),'r')
% plot(yc_true(1,:,1)*D_c(2) + temp_true(1,:,2),'k')
% plot(temp_true(1,:,4),'g')
%
% subplot(2,2,3)
% plot(xc_true(1,:,3))
% hold on
% plot(yc_true(1,:,3),'r')
% plot(temp_true(1,:,4),'g')
%
% subplot(2,2,4)
% plot(xc_true(1,:,4))
% hold on
% plot(yc_true(1,:,4),'r')
% % plot(temp_true(1,:,4),'g')
plot(yc_true(1,:,4),'-b*')
hold on
plot(yc_r(1,:,4),'-ro')
% plot(yc_pr(1,:,4), 'r');
% plot(yc_r(1,:,4), 'g');
overflow_count_xc
overflow_count_yc