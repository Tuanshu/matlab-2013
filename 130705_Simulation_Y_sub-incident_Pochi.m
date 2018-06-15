%% Options

clear all

C=3E8;

Frequency=[0:0.0002:7]*1E14;

N_f=length(Frequency);
N_t=100*N_f; %multi to N_f

n_air=1;

l_sam=10000;

l_ref=9000;

G=gaussmf(Frequency,[max(Frequency)/20 max(Frequency)/2]);    %incident power (spectral density)

%%
Time_total=1/(max(Frequency)/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;


%%

n1=1.5;
k1=0;
l1=500; %micron


n2=1.4;
k2=0;
l2=3;


n3=1.6;
k3=0;
l3=5;

%%

t= @ (n1,n2) (2*n1)./(n1+n2); 
t_r= @ (n1,n2) (2*n2)./(n1+n2);
r= @ (n1,n2) (n1-n2)./(n1+n2);
r_r= @ (n1,n2) (n2-n1)./(n1+n2);

%%
%%突然發現, 若是考慮material dispersion, 事情好會沒這麼單純,
%%因為遠處的訊號會被broaden
%%不過這次Pochi沒要我考慮dispersion
%%經過dispersion拓寬後影響也會降低, 所以應該還好吧(可以偷偷試著考慮一下)


% interface ref
s_ref=G.*exp(1i*4*pi*Frequency/C*l_ref/1E6)*r(n_air,n1);

% interface 1 (air-n2 interface)
s_1=G.*exp(1i*4*pi*Frequency/C*l_sam/1E6)*r(n_air,n1);

% interface 2(n2-n1 interface)
s_2=G.*exp(1i*4*pi*Frequency/C*l_sam/1E6)*t(n_air,n1)*t_r(n_air,n1).*exp(1i*4*pi*Frequency/C*n2*l2/1E6)*r(n2,n1);

% interface 3 (n1-n3 interface)
s_3=G.*exp(1i*4*pi*Frequency/C*l_sam/1E6)*t(n_air,n1)*t_r(n_air,n1).*exp(1i*4*pi*Frequency/C*n2*l2/1E6)*t(n2,n1)*t_r(n2,n1).*exp(1i*4*pi*Frequency/C*n1*l1/1E6)*r(n1,n3);

% interface 4 (n3-air interface)
s_4=G.*exp(1i*4*pi*Frequency/C*l_sam/1E6)*t(n_air,n1)*t_r(n_air,n1).*exp(1i*4*pi*Frequency/C*n2*l2/1E6)*t(n2,n1)*t_r(n2,n1).*exp(1i*4*pi*Frequency/C*n1*l1/1E6)*t(n1,n3)*t_r(n1,n3).*exp(1i*4*pi*Frequency/C*n3*l3/1E6)*r(n3,n_air);

S=(s_ref+s_1+s_2+s_3+s_4).*conj(s_ref+s_1+s_2+s_3+s_4);

FFT_S=fft(S,N_t);

plot(Frequency,S);

plot(Position_micron,real(FFT_S));