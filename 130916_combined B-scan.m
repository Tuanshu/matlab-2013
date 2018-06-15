clear all

array=1:1:3999;

Max_Wavelength=700;             %nm
Min_Wavelength=400;             %nm
N_f=8192;
N_t=8192*8;


cd('D:\Users\TuanShu\130916\600+\INTER_air_2\');
%% Data Loading
Data_First=importdata(sprintf('D%i.txt',array(1)));

Data(size(Data_First,1),length(array))=0;
Data(:,1)=Data_First(:,1);
for jj=1:length(array)

%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data_Temp=importdata(sprintf('D%i.txt',array(jj)));
Data(:,jj+1)=Data_Temp(:,2)-Data_First(1,2);
end

Wavelength=Data(:,1);           %nm

%Data=importdata('inter.txt')-importdata('sam.txt')-importdata('ref.txt')+importdata('bs.txt');
%cd('D:\120222\');
%Data_R=importdata('R1 Ref.txt');
%cd('D:\120222\R1 Sam\');
%Data_S=importdata(sprintf('D%i.txt',array(jj)));

C=3E8;

Frequency=C./(Wavelength*1E-9);

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_New=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency_New=Frequency_New';

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position_micron=C*Time(1:round(length(Time)))*1E6;

Spectra_Wavelength=Data(:,2:end);%-Data_R(:,2)-Data_S(:,2);
Spectra_Frequency=(Spectra_Wavelength.*((repmat(Wavelength,1,size(Spectra_Wavelength,2))*1E-9).^2)/C);
Spectra_New=interp1(Frequency,Spectra_Frequency,Frequency_New);

Spectra_New(isnan(Spectra_New))=0;
Spectra_New(Frequency_New<Min_Frequency)=0;
plot(Frequency_New,Spectra_New);

FFT_Spectra=fft(Spectra_New,[],2);
FFT_Spectra(:,1:100)=0;
FFT_Spectra(:,700:end)=0;

Spectra_New_Complex=ifft(FFT_Spectra,[],2);
Spectra_New_Complex((N_f+1):N_t,:)=0;
Signal=(fft(Spectra_New_Complex,[],1));

%dlmwrite('Lateral_Position.txt',Lateral_Position,'delimiter','\t','newline','pc');
%%
plot(real(Spectra_New(6000,:)));
plot(real(FFT_Spectra(6000,:)));

plot(real(Spectra_New_Complex(:,40)));
plot(Position_micron,(real(Signal(:,24))));
plot((abs(Signal)));
%%
plot(Position_micron(Position_micron>1550),(abs(Signal((Position_micron>1550),145))));

%%
[max_value max_index]=max(abs(Signal));
plot(max_index);
plot((abs(max_value)));
