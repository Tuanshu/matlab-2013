clear all;
close all;
clc;

N_f=4096*4;
N_t=N_f*16;

%temp=importdata('J:\100721_SDOCT with OO\100721_set2(1ms)_R+S(interfered)-3 (5times averaged).txt');
inter=importdata('J:\100816\100816_set2_align2_400mA_inter_1kHz_57deg.txt');
%ref=0;
ref=importdata('J:\100816\100816_set2_align2(ref)_400mA_inter_1kHz_57deg.txt');
sample=importdata('J:\100816\100816_set2_align2(sam)_400mA_inter_1kHz_57deg.txt');
%temp=temp.data;



pixel=1:4096;  %1~4096
[peak index_peak]=max(ref);
index_long=find(ref>0.5*max(ref),1,'last');
index_short=find(ref>0.5*max(ref),1,'first');
Q=141.31/(index_long-index_short);
lambda=((-(pixel-index_peak)*Q+745))';   %前面的負號很重要! 影響freq domain是否接近guassian (不過為什麼差一個負號dispersion的broaden效應好像也會變大? 因為carrier in lambda domain有chirp, 但在freq domain應該不太有)

%ref=0;
%sample=0;

S0=inter-ref-sample;
%S0=inter;
%plot(lambda,S0,lambda,inter,lambda,ref);
S0=S0/max(S0);


c=3E8;              %m/sec
freq=c./(lambda*1E-9);     %Hz
d_f=max(freq)/(N_f-1);
fx=0:d_f:max(freq);
S=interp1(freq,S0,fx);
S(isnan(S))=0;

zeros(1:N_f)=0;
S_padded=[S zeros];             %with minus frequency, 2*N_f
CS=real(fft(S_padded,N_t))';     %with minus time
%CS=real(fft(S_padded));     %with minus time
CS_normal=CS/max(abs(CS));
d_t=1/(d_f*N_t);
%d_t=1/(d_f*2*N_f);
time=[-0.5*(N_t-1)*d_t:d_t:0.5*N_t*d_t]'/2;%/2是因為一來一回
%time=[-(N_f-1)*d_t:d_t:N_f*d_t]/2;  %/2是因為一來一回
space=c*time;
CS_envelope=abs(hilbert(CS_normal));
plot(space,CS_normal,space,CS_envelope);

space_min=-1.782E-3;
space_max=-1.776E-3;
%space_min=-2.2E-3;
%space_max=-2.15E-3;
space_min_index=find(space>space_min, 1, 'first');
space_max_index=find(space>space_max, 1, 'first');
space=space(space_min_index:space_max_index);
time=time(space_min_index:space_max_index);
CS_normal=CS_normal(space_min_index:space_max_index);
CS_envelope=CS_envelope(space_min_index:space_max_index);
[peakvalue peakindex]=max(CS_envelope);
FWHM_max=space(find(CS_envelope>0.5*peakvalue, 1, 'last'));
FWHM_min=space(find(CS_envelope>0.5*peakvalue, 1, 'first'));
FWHM=FWHM_max-FWHM_min;
M=[space CS_normal CS_envelope];
plot(space,CS_normal,space,CS_envelope);
%dlmwrite('Interferogram.txt',M,'delimiter','\t','newline','pc');
%BW=lambda(find(S0>0.5,1,'last'))-lambda(find(S0>0.5,1,'first'));
%x_Res=x(find(CS_envelope>0.5,1,'last'))-x(find(CS_envelope>0.5,1,'first'));