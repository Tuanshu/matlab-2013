clear all

%pi*(w_2^2)/lambda=((pi*(w_1^2)/lambda)*f^2)/((f-d)^2+(pi*(w_1^2)/lambda)^2)

w_1=1;      %mm
f=200;       %mm
lambda=0.55;%micron
w_2=lambda*f/pi/w_1;%micron

Zr=pi*(w_2^2)/lambda/1000;%mm