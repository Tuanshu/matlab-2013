clear all

%%% To find the phase from min voltage to max voltage ONLY

Path='D:\Users\TuanShu\130419_Glass\Data44.txt';

Data=dlmread(Path);

Position=(1:length(Data))/20000*2000;

plot(Data);

FFT=fft(Data);

plot(real(FFT));

FFT(1:36000)=0;

FFT(fix(length(FFT)/2):end)=0;
%FFT(60000:end)=0;

Data=ifft(FFT);
Data((length(Data)-100):end)=0;
Data(1:100)=0;
[maxvalue maxindex]=max(Data);
Position=Position-Position(maxindex);
Data=Data/maxvalue;
plot(Position,abs(Data),Position,real(Data));
xlabel('Position (micron)');
ylabel('Signal Intensity (norm.)');


%%

Data(Data>7000,2)=0;
Data(Data<5000,2)=0;


Data(:,2)=Data(:,2)/max(abs(Data(:,2)));

FWHM=abs(Data(find(Data(:,2)>0.5,1,'first'),1)-Data(find(Data(:,2)>0.5,1,'last'),1));

[maxvalue maxindex]=max(Data(:,2));

plot(Data(:,1)-Data(maxindex,1),(Data(:,2)),Data(:,1)-Data(maxindex,1),abs(Data(:,2)));
xlabel('OPD (Micron)');
ylabel('a.u.');