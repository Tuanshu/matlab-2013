clear all

cd('D:\Users\TuanShu\131007\');
C=3E8;

LP=80;

SP=20;

Thickness_1=17.94*1E-6;
Thickness_2= 39.52*1E-6;

M=520;  %which line for calibration
N=815; %which pixel is center

Frequency_max=8E14;

N_f=8192;

fx=0:Frequency_max/(N_f-1):Frequency_max;
fx=importdata('fx.txt');


Data1_2D=importdata('3');
%OPD4=47.769897

Data1_2D=rot90(Data1_2D);

Data2_2D=importdata('4');
%OPD1=61.605919

Data2_2D=rot90(Data2_2D);


image(Data1_2D);
colormap(gray);
%%
Data1=Data1_2D(M,:);
Data2=Data2_2D(M,:);

X=1:length(Data1);

plot(X,Data1,X,Data2);
%%

FFT1=fft(Data1);


FFT1(1:SP)=0;
FFT1(LP:end)=0;
Data1_H=ifft(FFT1);
PHI1=unwrap(angle(Data1_H));

FFT2=fft(Data2);
FFT2(1:SP)=0;
FFT2(LP:end)=0;
Data2_H=ifft(FFT2);
PHI2=unwrap(angle(Data2_H));

Delta_PHI=PHI2-PHI1;

plot(X,PHI1,X,PHI2);

%%
Freq=Delta_PHI*C/(Thickness_2-Thickness_1)/4/pi;

plot(X,Delta_PHI);
plot(Freq,Delta_PHI);

Freq=Freq-Freq(N)+C/(780E-9);

plot(X,Freq);
%plot(Freq);
xlabel('CCD Pixel Number');
ylabel('Optical Frequency (Hz)');

Calibration_array=Freq';

dlmwrite('Calibration_array.txt',Calibration_array,'delimiter','\t','newline','pc');

plot(Calibration_array,Data1);

for p=1:size(Data1_2D,1)
    Data1_2D_New(p,:)=interp1(Calibration_array,Data1_2D(p,:),fx);
    Data2_2D_New(p,:)=interp1(Calibration_array,Data2_2D(p,:),fx);
    disp(p);
end
Data1_2D_New(isnan(Data1_2D_New))=0;
Data2_2D_New(isnan(Data2_2D_New))=0;

image(Data2_2D_New);
colormap(gray);
FFT_1 = abs(fft(Data1_2D_New,[],2));
FFT_2 = abs(fft(Data2_2D_New,[],2));

image(FFT_1);
colormap(gray);

c=3E8;
Time_total=1/(max(fx)/(length(fx)-1));
Time=[0:Time_total/(length(fx)-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=c*Time;
Position_micron=Position*1E6;
d_Position_micron=Position_micron(2)-Position_micron(1);