clear all

%% Setting
SR=200000;
Stage_Speed=2;  %mm/s

Stage_Speed_MTS=(Stage_Speed*1E-3);

C=3E8;

TH_first_peak=0.01;
TH_second_peak=0.001;
Delay_first_peak=80000; %the smaller one, pixel
Min_thickness=12000; %pixel

Max_thickness=16000; %pixel

cd('D:\Users\TuanShu\130513\');

qqq=1;
xxx=1;
yyy=1;
j=14;     %14 or 15: first peak, 11 or 12: rear peak

%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata(sprintf('Data%d.txt',j));

N_t=length(Data);
N_f=N_t;

Time_Stage=1/SR:1/SR:(1/SR)*N_t;

Time=Time_Stage*Stage_Speed_MTS/C;
Position_micron=Time*C*1E6;

dTime=(Time(2)-Time(1))*2;      %*2 becuase of round trip

Frequency_Max=1/dTime;

Frequency=Frequency_Max/N_f:Frequency_Max/N_f:Frequency_Max;

Wavelength_micron=(C./Frequency)*1E6;

Spectrum=fft(Data,N_f);


Window1=(gaussmf(Frequency,[0.2E14 2.2E14]));
Window1(Frequency>2.2E14)=1;
Window2=(gaussmf(Frequency,[0.2E14 3.6E14]));
Window2(Frequency<3.6E14)=1;
Window=(Window1.*Window2)';
Spectrum=Window.*Spectrum;
Spectrum((round(length(Spectrum)/2)+1):end)=0;

plot(Frequency,Spectrum);

Data_New=ifft(Spectrum);
Data_New=Data_New(1:N_t);

Max_Wavelength_micron=1;

Min_Wavelength_micron=0.6;

Max_Wavelength_micron_index=find(Wavelength_micron<Max_Wavelength_micron,1,'first');

Min_Wavelength_micron_index=find(Wavelength_micron<Min_Wavelength_micron,1,'first');

Data_New(1:40)=Data_New(40);
Data_New((length(Data_New)-40):end)=Data_New((length(Data_New)-40));

plot(Wavelength_micron(Max_Wavelength_micron_index:Min_Wavelength_micron_index),Spectrum(Max_Wavelength_micron_index:Min_Wavelength_micron_index));

xlabel('Wavelength (micron)');
ylabel('Interference Power (a.u.)');

[C ind]=max(abs(Data_New));
Position_micron=Position_micron-Position_micron(ind);
plot(Position_micron(abs(Position_micron-200)<800),Data_New(abs(Position_micron-200)<800)/max(abs(Data_New)),Position_micron(abs(Position_micron-200)<800),abs(Data_New(abs(Position_micron-200)<800))/max(abs(Data_New)));
xlabel('OPD (micron)');
ylabel('Interference Signal');
%max_array_even_filtered=max_array_even-(max_array_even_fitting(1).*[1:length(max_array_even)]+max_array_even_fitting(2));

%max_array_odd_filtered=max_array_odd-(max_array_odd_fitting(1).*[1:length(max_array_odd)]+max_array_odd_fitting(2));


%SD=std(abs(Data_New(1.26E5:1.36E5)));

%SNR_even=10*log10((max(abs(max_array_even))^2)/(SD^2));
%SNR_odd=10*log10((max(abs(max_array_odd))^2)/(SD^2));
%plot(Wavelength_micron,Spectrum);
%xlabel('Wavelength (micron)');
%label('Interference Spectrum');
%cd('D:\120222\');

%dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');
%plot(Position_micron,Signal_Bscan_Envelope);