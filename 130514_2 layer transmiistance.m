clear all

n=1.6;
n1=1.5;
n2=1;

t_AN100=2*n2/(n2+n1);   %0.8
t_1=2*n1/(n1+n);        %0.96
t_2=2*n/(n+n2);         %1.23

T=(t_AN100*t_1*t_2)^2   %0.9