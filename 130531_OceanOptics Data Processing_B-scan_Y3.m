clear all

%% Setting
Motor_Speed=0.5;%0.0025;   %mm/sec
Integration_Time=5;   %40%total integration time,ms

Pixel_Average_Axial=1;
Pixel_Average_Lateral=10;

Lateral_Spacing=Motor_Speed*Integration_Time;   %micron

Max_Wavelength=700;             %nm
Min_Wavelength=400;             %nm
N_f=8192;
N_t=8192*8;
ROI_ratio=1/16;                  %only consider the first ROI_ratio data in TD

DC_Cutoff=5;                 %micron

array=0:1999;

cd('D:\Users\TuanShu\130531\YPR 3\');
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
Signal=downsample(conv(abs(Signal),(ones(Pixel_Average_Axial,1))/Pixel_Average_Axial,'same'),Pixel_Average_Axial);

Signal_Bscan_Envelope(:,jj)=Signal(1:size(Signal_Bscan_Envelope,1));

if jj==1

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position_micron=downsample(C*Time(1:round(length(Time)*ROI_ratio)),Pixel_Average_Axial)*1E6;
Position_micron=Position_micron(1:size(Signal_Bscan_Envelope,1));

end

Signal_Bscan_Envelope(Position_micron<DC_Cutoff,jj)=0;

disp(jj)
end
Signal_Bscan_Envelope=Signal_Bscan_Envelope/max(max(Signal_Bscan_Envelope));
for pp=1:size(Signal_Bscan_Envelope,1)
    Signal_Bscan_Envelope(pp,:)=conv(Signal_Bscan_Envelope(pp,:),ones(Pixel_Average_Lateral,1)/Pixel_Average_Lateral,'same');
end
Lateral_Position=[Lateral_Spacing:Lateral_Spacing:size(Signal_Bscan_Envelope,2)*Lateral_Spacing]';



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

plot(Lateral_Position,Signal_Bscan_Envelope(maxindex,:));

%% Circshift
Range_1=45;  %micron

Range_2=55.3;


[maxvalue maxindex]=max(Signal_Bscan_Envelope(find(Position_micron>Range_1,1,'first'):find(Position_micron>Range_2,1,'first'),:),[],1);

needshift=maxindex(1)-maxindex;

absolute_index_record=find(Position_micron>Range_1,1,'first')+maxindex(1)-1;

%needshift=round(conv(needshift,(ones(20,1))/20,'same'));

for pp=1:length(needshift)
    Signal_Bscan_Envelope_Shifted(:,pp)=circshift(Signal_Bscan_Envelope(:,pp),needshift(pp));
end

imagesc(Signal_Bscan_Envelope_Shifted,'xdata',Lateral_Position,'ydata',Position_micron);
%colormap(gray);
%axis equal
xlabel('Lateral Position (micron)');
ylabel('Optical Path Difference (micron)');

%% Optical Thickness Measurement

Range_3=55;  %micron

Range_4=62;


[maxvalue2 maxindex2]=max(Signal_Bscan_Envelope_Shifted(find(Position_micron>Range_3,1,'first'):find(Position_micron>Range_4,1,'first'),:),[],1);

absolute_index_2_record=find(Position_micron>Range_3,1,'first')+maxindex2-1;

Signal_Bscan_Envelope_Shifted=Signal_Bscan_Envelope_Shifted/max(maxvalue2);
SD=std(Position_micron(maxindex2));

Thickness_Profile=Position_micron(absolute_index_2_record)-Position_micron(absolute_index_record);
plot(Lateral_Position,Thickness_Profile);
dlmwrite('Lateral_Position.txt',Lateral_Position,'delimiter','\t','newline','pc');
dlmwrite('Thickness_Profile.txt',Thickness_Profile,'delimiter','\t','newline','pc');

%% line by line normalization
for qq=1:size(Signal_Bscan_Envelope_Shifted,2)
    Signal_Bscan_Envelope_Shifted(:,qq)=Signal_Bscan_Envelope_Shifted(:,qq)/max(Signal_Bscan_Envelope_Shifted(find(Position_micron>Range_3,1,'first'):find(Position_micron>Range_4,1,'first'),qq));
end
%% Show

Range_3=46;  %micron

Range_4=66;
set(gcf, 'Position',[200 200 500 200]);
set(gca, 'Position',[0 0 1 1]);
%imagesc((Signal_Bscan_Envelope_Shifted(find(Position_micron>Range_3,1,'first'):find(Position_micron>Range_4,1,'first'),:)),'xdata',Lateral_Position,'ydata',Position_micron(find(Position_micron>Range_3,1,'first'):find(Position_micron>Range_4,1,'first')));
imagesc((10*log10(Signal_Bscan_Envelope_Shifted(find(Position_micron>Range_3,1,'first'):find(Position_micron>Range_4,1,'first'),:))),'xdata',Lateral_Position,'ydata',Position_micron(find(Position_micron>Range_3,1,'first'):find(Position_micron>Range_4,1,'first')));
colormap(gray);
caxis([-20 0]);
%axis equal
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca,'XColor','white');
set(gca,'YColor','white');
%xlabel('Lateral Position (micron)');
%ylabel('Optical Path Difference (micron)');
%dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');
%plot(Position_micron,Signal_Bscan_Envelope);