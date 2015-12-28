clc;close all;clear
f = fipref('loggingmode','on');
k = 0.0888;
num = [1 4.43 10.76 17.46 20.48 17.46 10.76 4.43 1]*k;
den = [1 1.10 1.97 1.55 1.22 0.61 0.24 0.061 0.008];
Hcas.arithmetic = 'fixed';


SOS_c = tf2sos(num,den);
H1 = dfilt.df2t(SOS_c(1,1:3), SOS_c(1,4:6));
H2 = dfilt.df2t(SOS_c(2,1:3), SOS_c(2,4:6));
H3 = dfilt.df2t(SOS_c(3,1:3), SOS_c(3,4:6));
H4 = dfilt.df2t(SOS_c(4,1:3), SOS_c(4,4:6));

H1.Arithmetic = 'fixed';
Hcas = dfilt.cascade(H1, H2, H3, H4);
u = ones(1,1000);
set(H1, 'OverflowMode', 'Saturate');
y = filter(H1,u);
rlog = qreport(H1)
get(H1)