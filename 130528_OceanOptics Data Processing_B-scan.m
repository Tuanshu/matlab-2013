clear all

%% Setting
Motor_Speed=0.1;%0.0025;   %mm/sec
Integration_Time=8;   %40%total integration time,ms

Pixel_Average_Axial=2;
Lateral_Spacing=Motor_Speed*Integration_Time;   %micron

Max_Wavelength=700;             %nm
Min_Wavelength=400;             %nm
N_f=8192;
N_t=8192*32;
ROI_ratio=1/16;                  %only consider the first ROI_ratio data in TD

DC_Cutoff=5;                 %micron

array=0:999;

cd('D:\Users\TuanShu\130528\GPR 1\');
%% Data Loading

Signal_Bscan_Envelope(1:round(N_t*ROI_ratio/Pixel_Average_Axial),1:length(array))=0;

for jj=1:length(array)

%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata(sprintf('D%i.txt',array(jj)));

Wavelength=Data(:,1);           %nm

%Data=importdata('inter.txt')-importdata('sam.txt')-importdata('ref.txt')+importdata('bs.txt');
%cd('D:\120222\');
%Data_R=importdata('R1 Ref.txt');
%cd('D:\120222\R1 Sam\');
%Data_S=importdata(sprintf('D%i.txt',array(jj)));

if jj==1

C=3E8;

Frequency=C./(Wavelength*1E-9);

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_New=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency_New=Frequency_New';

end

Spectrum=Data(:,2)-Data(1,2);%-Data_R(:,2)-Data_S(:,2);
Spectrum_Frequency=(Spectrum.*((Wavelength*1E-9).^2)/C)/max(Spectrum.*((Wavelength*1E-9).^2)/C);
Spectrum_New=interp1(Frequency,Spectrum_Frequency,Frequency_New);

Spectrum_New(isnan(Spectrum_New))=0;
Spectrum_New(Frequency_New<Min_Frequency)=0;
%plot(Frequency_New,Spectrum_New);

%% To time domain

Spectrum_New((N_f+1):N_t)=0;
Signal=fft(Spectrum_New);

Spectrum_New=Spectrum_New(1:N_f);
Signal=Signal(1:round(length(Signal)*ROI_ratio));
Signal_Bscan_Envelope(:,jj)=downsample(conv(abs(Signal),(ones(Pixel_Average_Axial,1))/Pixel_Average_Axial,'same'),Pixel_Average_Axial);

if jj==1

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=downsample(C*Time(1:round(length(Time)*ROI_ratio)),Pixel_Average_Axial);
Position_micron=Position*1E6;

end

Signal_Bscan_Envelope(Position_micron<DC_Cutoff,jj)=0;

disp(jj)
end
Signal_Bscan_Envelope=Signal_Bscan_Envelope/max(max(Signal_Bscan_Envelope));

plot(Position_micron,Signal_Bscan_Envelope(:,2));
xlabel('OPD (micron)');
ylabel('Interference Signal');

%axis equal
imagesc(Signal_Bscan_Envelope,'xdata',Lateral_Spacing:Lateral_Spacing:size(Signal_Bscan_Envelope,2)*Lateral_Spacing,'ydata',Position_micron);
%colormap(gray);
%axis equal
xlabel('Lateral Position (micron)');
ylabel('Optical Path Difference (micron)');



[maxvalue maxindex]=max(Signal_Bscan_Envelope(:,1));

plot(Lateral_Spacing:Lateral_Spacing:size(Signal_Bscan_Envelope,2)*Lateral_Spacing,Signal_Bscan_Envelope(maxindex,:));

%% Circshift
Range_1=35;  %micron

Range_2=45;


[maxvalue maxindex]=max(Signal_Bscan_Envelope(find(Position_micron>Range_1,1,'first'):find(Position_micron>Range_2,1,'first'),:),[],1);

needshift=maxindex(1)-maxindex;

%needshift=round(conv(needshift,(ones(20,1))/20,'same'));

for pp=1:length(needshift)
    Signal_Bscan_Envelope_Shifted(:,pp)=circshift(Signal_Bscan_Envelope(:,pp),needshift(pp));
end

imagesc(Signal_Bscan_Envelope_Shifted,'xdata',Lateral_Spacing:Lateral_Spacing:size(Signal_Bscan_Envelope,2)*Lateral_Spacing,'ydata',Position_micron);
%colormap(gray);
%axis equal
xlabel('Lateral Position (micron)');
ylabel('Optical Path Difference (micron)');

%% Optical Thickness Measurement (SD only)

Range_3=46;  %micron

Range_4=52;


[maxvalue2 maxindex2]=max(Signal_Bscan_Envelope_Shifted(find(Position_micron>Range_3,1,'first'):find(Position_micron>Range_4,1,'first'),:),[],1);

Signal_Bscan_Envelope_Shifted=Signal_Bscan_Envelope_Shifted/max(maxvalue2);
SD=std(Position_micron(maxindex2))

%% Show

Range_3=15;  %micron

Range_4=75;

imagesc((Signal_Bscan_Envelope_Shifted(find(Position_micron>Range_3,1,'first'):find(Position_micron>Range_4,1,'first'),:)),'xdata',Lateral_Spacing:Lateral_Spacing:size(Signal_Bscan_Envelope,2)*Lateral_Spacing,'ydata',Position_micron(find(Position_micron>Range_3,1,'first'):find(Position_micron>Range_4,1,'first')));
colormap(gray);
%axis equal
set(gca, 'XTick', []);
set(gca, 'YTick', []);
xlabel('Lateral Position (micron)');
ylabel('Optical Path Difference (micron)');
%dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');
%plot(Position_micron,Signal_Bscan_Envelope);