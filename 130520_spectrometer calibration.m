clear all

cd('D:\Users\TuanShu\130520\CCD\');
C=3E8;

LP=2000;

SP=40;

Data1=importdata('ccd4.txt');
%OPD4=47.769897

Data2=importdata('ccd1.txt');
%OPD1=61.605919

Data3=importdata('ccd5.txt');
%OPD5=126.923941

X=1:length(Data1);

plot(X,Data1,X,Data2,X,Data3);

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

FFT3=fft(Data3);
FFT3(1:SP)=0;
FFT3(LP:end)=0;
Data3_H=ifft(FFT3);
PHI3=unwrap(angle(Data3_H));

Delta_PHI1=PHI2-PHI1;
%13.8360
Delta_PHI2=PHI3-PHI1;
%79.1540

plot(X,PHI2,X,PHI3);
FX1=Delta_PHI1*C/(13.8360E-6)/4/pi;
FX2=Delta_PHI2*C/(79.1540E-6)/4/pi;
plot(X,Delta_PHI1,X,Delta_PHI2);
plot(FX1,Delta_PHI2);


plot(FX1,Delta_PHI2);
FX1=FX1-FX1(2200)+C/(780E-9);
FX2=FX2-FX2(2200)+C/(780E-9);
plot(X,FX1,X,FX2);
plot(FX1);
xlabel('CCD Pixel Number');
ylabel('Optical Frequency (Hz)');
A=fit(X(2000:2400)',FX1(2000:2400),'poly1')

index=[0:4095]';

Calibration_array=A.p1*index+A.p2;

dlmwrite('Calibration_array.txt',Calibration_array,'delimiter','\t','newline','pc');
