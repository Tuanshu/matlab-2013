clear all

%%% To find the phase from min voltage to max voltage ONLY

Path='D:\Users\TuanShu\130319_Yb-clad TD Signal.txt';

Data=dlmread(Path);

plot(Data(:,1),Data(:,2));

FFT=fft(Data(:,2));

plot(real(FFT));

FFT(1:8000)=FFT(8001);

FFT(60000:end)=0;

Data(:,2)=ifft(FFT);

Data(Data(:,1)>7000,2)=0;
Data(Data(:,1)<5000,2)=0;


Data(:,2)=Data(:,2)/max(abs(Data(:,2)));

FWHM=abs(Data(find(Data(:,2)>0.5,1,'first'),1)-Data(find(Data(:,2)>0.5,1,'last'),1));

[maxvalue maxindex]=max(Data(:,2));

plot(Data(:,1)-Data(maxindex,1),(Data(:,2)),Data(:,1)-Data(maxindex,1),abs(Data(:,2)));
xlabel('OPD (Micron)');
ylabel('a.u.');