%% Options

clear all

C=3E8;

If_dispersive=0;
If_Bandlimited=1;

Shortest_Wavelength=0.3;        %micron
N_f=50000;
N_t=100*N_f; %multi to N_f

Center_Wavelength=0.55;          %nm
Bandwidth_Wavelength=0.09;                   

Center_Frequency=C/(Center_Wavelength*1E-6);
Bandwidth_Frequency=Bandwidth_Wavelength/Center_Wavelength*Center_Frequency;

Largest_Frequency=C/(Shortest_Wavelength*1E-6);
Spectral_Resolution_Hz=Largest_Frequency/N_f;
Frequency=[0:fix(Largest_Frequency/(N_f)):fix(Largest_Frequency/(N_f))*(N_f-1)];
Wavelength_micron=C./Frequency*1E6;
n_air=1;

l_sam=100000;

l_ref=99000;
if If_Bandlimited==1
    G=gaussmf(Frequency,[Bandwidth_Frequency/(2*(2*log(2))^0.5) Center_Frequency]);    %incident power (spectral density)
else
    G=1;
end

%%
Time_total=1/(max(Frequency)/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;

%%

n_BK7=real(1+1.03961212*(Wavelength_micron.^2)./((Wavelength_micron.^2)-0.00600069867) + 0.231792344*(Wavelength_micron.^2)./((Wavelength_micron.^2)-0.0200179144) + 1.01046945*(Wavelength_micron.^2)./((Wavelength_micron.^2)-103.560653)).^0.5;


%%
if If_dispersive==0
    n1=1.5;
else
    n1=1.5+n_BK7-n_BK7(find(Wavelength_micron<Center_Wavelength,1,'first'));
    n1(isnan(n1))=1.5;
end
k1=0;
l1=500; %micron

if If_dispersive==0
    n2=1.4;
else
    n2=1.4+n_BK7-n_BK7(find(Wavelength_micron<Center_Wavelength,1,'first'));
    n2(isnan(n2))=1.4;
end
k2=0;
l2=3;


if If_dispersive==0
    n3=1.6;
else
    n3=1.6+n_BK7-n_BK7(find(Wavelength_micron<Center_Wavelength,1,'first'));
    n3(isnan(n3))=1.6;
end
k3=0;
l3=5;


%%

t= @ (n1,n2) (2*n1)./(n1+n2); 
t_r= @ (n1,n2) (2*n2)./(n1+n2);
r= @ (n1,n2) (n1-n2)./(n1+n2);
r_r= @ (n1,n2) (n2-n1)./(n1+n2);
d= @ (n,l) exp(1i*4*pi.*Frequency/C.*n*l/1E6);

%%

% FP of 1st film

r_fp2=r(n_air,n2)+(t(n_air,n2).*t_r(n_air,n2).*r(n2,n1).*d(n2,l2))./(1-r_r(n_air,n2).*r(n2,n1).*d(n2,l2));
t_fp2=(t(n_air,n2).*t(n2,n1).*exp(1i*2*pi.*Frequency/C.*n2*l2/1E6))./(1-r_r(n_air,n2).*r(n2,n1).*d(n2,l2));


r_r_fp2=r(n1,n2)+(t(n1,n2).*t_r(n1,n2).*r(n2,n_air).*d(n2,l2))./(1-r_r(n1,n2).*r(n2,n_air).*d(n2,l2));

t_r_fp2=(t(n1,n2).*t(n2,n_air).*exp(1i*2*pi.*Frequency/C.*n2*l2/1E6))./(1-r_r(n1,n2).*r(n2,n_air).*d(n2,l2));
%%
% FP of 2nd film
r_fp3=r(n1,n3)+(t(n1,n3).*t_r(n1,n3).*r(n3,n_air).*d(n3,l3))./(1-r_r(n1,n3).*r(n3,n_air).*d(n3,l3));

t_fp3=(t(n1,n3).*t(n3,n_air).*exp(1i*2*pi.*Frequency/C.*n3*l3/1E6))./(1-r_r(n1,n3).*r(n3,n_air).*d(n3,l3));

r_r_fp3=r(n_air,n3)+(t(n_air,n3).*t_r(n_air,n3).*r(n3,n1).*d(n3,l3))./(1-r_r(n_air,n3).*r(n3,n1).*d(n3,l3));
t_r_fp3=(t(n_air,n3).*t(n3,n1).*exp(1i*2*pi.*Frequency/C.*n3*l3/1E6))./(1-r_r(n_air,n3).*r(n3,n1).*d(n3,l3));
% coupled FP of both FPs

r_total=r_fp2+(t_fp2.*t_r_fp2.*r_fp3.*d(n1,l1))./(1-r_r_fp2.*r_fp3.*d(n1,l1));
%%
%%突然發現, 若是考慮material dispersion, 事情好會沒這麼單純,
%%因為遠處的訊號會被broaden
%%不過這次Pochi沒要我考慮dispersion
%%經過dispersion拓寬後影響也會降低, 所以應該還好吧(可以偷偷試著考慮一下)


% interface ref
s_ref=((G).^0.5).*exp(1i*4*pi*Frequency/C*l_ref/1E6);
%s_ref=s_ref;

s_sam=((G).^0.5).*exp(1i*4*pi*Frequency/C*l_sam/1E6).*r_total;
%s_sam=s_sam;
% interface 1 (air-n2 interface)
%s_1=G.*exp(1i*4*pi*Frequency/C*l_sam/1E6).*r(n_air,n1);

% interface 2(n2-n1 interface)
%s_2=G.*exp(1i*4*pi*Frequency/C*l_sam/1E6).*t(n_air,n1).*t_r(n_air,n1).*exp(1i*4*pi*Frequency/C.*n2*l2/1E6).*r(n2,n1);

% interface 3 (n1-n3 interface)
s_3=G.*exp(1i*4*pi*Frequency/C*l_sam/1E6).*t(n_air,n1).*t_r(n_air,n1).*exp(1i*4*pi*Frequency/C.*n2*l2/1E6).*t(n2,n1).*t_r(n2,n1).*exp(1i*4*pi*Frequency/C.*n1*l1/1E6).*r(n1,n3);

% interface 4 (n3-air interface)
s_4=G.*exp(1i*4*pi*Frequency/C*l_sam/1E6).*t(n_air,n1).*t_r(n_air,n1).*exp(1i*4*pi*Frequency/C.*n2*l2/1E6).*t(n2,n1).*t_r(n2,n1).*exp(1i*4*pi*Frequency/C.*n1*l1/1E6).*t(n1,n3).*t_r(n1,n3).*exp(1i*4*pi*Frequency/C.*n3*l3/1E6).*r(n3,n_air);

s_target=s_3+s_4;
%S=(s_ref+s_1+s_2+s_3+s_4).*conj(s_ref+s_1+s_2+s_3+s_4);

S=(s_ref+s_sam).*conj(s_ref+s_sam);
%S_res=(s_ref+s_sam-s_target).*conj(s_ref+s_sam-s_target);
%S_sam=abs(s_sam).^2;
%S_inter=(s_sam).*conj(s_ref)+conj(s_sam).*(s_ref);

FFT_S=fft(S,N_t);
%FFT_S_target=fft(S_target,N_t);
%FFT_res=fft(S_res);
%FFT_S_inter=fft(S_inter,N_t);
%FFT_S_sam=fft(S_sam,N_t);
%FFT_f=fft(S_inter./G,N_t);
%%

subplot(4,1,1)
plot(Position_micron,real(FFT_S));
xlabel('Optical Path Difference (micron)');
ylabel('Signal (a.u.)');
xlim([0 fix(max(Position_micron)/2)]);
ylim([-1200 1200]);

%subplot(4,1,2)
%plot(Position_micron,real(FFT_S_inter));
%xlabel('Optical Path Difference (micron)');
%ylabel('Signal (a.u.)');
%xlim([0 fix(max(Position_micron)/2)]);

subplot(4,1,2)
plot(Position_micron,real(FFT_S));
xlabel('Optical Path Difference (micron)');
ylabel('Signal (a.u.)');
xlim([l_sam-l_ref-50 l_sam-l_ref+850]);
ylim([-1200 1200]);

subplot(4,1,3)
plot(Position_micron,real(FFT_S));
xlabel('Optical Path Difference (micron)');
ylabel('Signal (a.u.)');
xlim([l_sam-l_ref-50 l_sam-l_ref+50]);
ylim([-1200 1200]);

subplot(4,1,4)
plot(Position_micron,real(FFT_S));
xlabel('Optical Path Difference (micron)');
ylabel('Signal (a.u.)');
xlim([l_sam-l_ref+710 l_sam-l_ref+810]);
ylim([-1200 1200]);


plot(Position_micron,real(FFT_S));
xlabel('Optical Path Difference (micron)');
ylabel('Signal (a.u.)');
xlim([l_sam-l_ref+710 l_sam-l_ref+810]);
ylim([-1200 1200]);

%%
QQ=fftshift(real(FFT_S));
QQQ=QQ;

QQQ((fix(length(Position_micron)/2)-8000):(fix(length(Position_micron)/2)+8000))=0;

QQQQ=2*QQQ;
QQQQ(1:fix(length(Position_micron)/2))=0;
plot(Position_micron-Position_micron(fix(length(Position_micron)/2)),QQ);
xlabel('Optical Path Difference (micron)');
ylabel('Signal (a.u.)');


plot(Position_micron-Position_micron(fix(length(Position_micron)/2)),QQQ);
xlabel('Optical Path Difference (micron)');
ylabel('Signal (a.u.)');


plot(Position_micron-Position_micron(fix(length(Position_micron)/2)),QQQQ);
xlabel('Optical Path Difference (micron)');
ylabel('Signal (a.u.)');