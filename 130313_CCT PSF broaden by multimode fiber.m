clear all

%% Setting


thickness_temp=2.68E-6;
n_should=1.5;

Starting_pixel_f_considered=250;
Ending_pixel_f_considered=350;

range_specified=1;

Max_Wavelength=1100;             %nm
Min_Wavelength=900;             %nm
N_f=8192;
N_t=8192*8;
ROI_ratio=1/4;                  %only consider the first ROI_ratio data in TD


%% Cornea Setting

Spot_Size=4:20:4002;  %micron
%Spot_Size=2900;

N_Grid=100;                 %min N
R_Grid=0.1;                   %micron

RoC=7.8;    %mm

RoC_Rear=6.4;    %mm


%% Data Loading


cd('D:\Users\TuanShu\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata('130311_YbSpectrum.txt');
%cd('D:\120222\');
%Data_R=importdata('R1 Ref.txt');
%cd('D:\120222\R1 Sam\');
%Data_S=importdata(sprintf('D%i.txt',array(jj)));

    
Wavelength=Data(:,1);           %nm

C=3E8;

Frequency=C./(Wavelength*1E-9);

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_New=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency_New=Frequency_New';


Spectrum=Data(:,2)-Data(1,2);%-Data_R(:,2)-Data_S(:,2);
Spectrum_Frequency=(Spectrum.*((Wavelength*1E-9).^2)/C)/max(Spectrum.*((Wavelength*1E-9).^2)/C);
Spectrum_New=interp1(Frequency,Spectrum_Frequency,Frequency_New);

Spectrum_New(isnan(Spectrum_New))=0;
Spectrum_New(Frequency_New<Min_Frequency)=0;
plot(Frequency_New,Spectrum_New);

%% To time domain

Spectrum_New((N_f+1):N_t)=0;



Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;
Position_micron_Scale=Position_micron;
Signal=fft(Spectrum_New);
Spectrum_New=Spectrum_New(1:N_f);
Signal=fftshift(Signal)./max(abs(fftshift(Signal)));
[maxvalue maxindex]=max(abs(Signal));
Position_micron=Position_micron-Position_micron(maxindex);

%% Signal Weighting

for q=1:length(Spot_Size)
    if fix((Spot_Size(q)/2)/R_Grid)>N_Grid
        R=R_Grid:R_Grid:(Spot_Size(q)/2);
    else
        R=((Spot_Size(q)/2)/N_Grid):((Spot_Size(q)/2)/N_Grid):(Spot_Size(q)/2);
    end
    OPD_Mapping=abs(((RoC*1E3)^2-R.^2).^0.5-RoC*1E3);
    OPD_Mapping_Rear=abs(((RoC_Rear*1E3)^2-R.^2).^0.5-RoC_Rear*1E3);
    Weight_Function=(R.*2*pi.*R)./sum((R.*2*pi.*R));

    Signal_Sum=0;    
    Signal_Sum_Rear=0;
    for p=1:length(R)
        Signal_Sum=Signal_Sum+Weight_Function(p)*circshift(Signal,find(Position_micron_Scale>OPD_Mapping(p),1,'first'));
        Signal_Sum_Rear=Signal_Sum_Rear+Weight_Function(p)*circshift(Signal,find(Position_micron_Scale>OPD_Mapping_Rear(p),1,'first'));
        fprintf('N_R=%d/%d, N_Spot=%d/%d\n',p,length(R),q,length(Spot_Size));
    end
    Signal_Env=abs(Signal_Sum)/max(abs(Signal_Sum));
    Signal_Env_Rear=abs(Signal_Sum_Rear)/max(abs(Signal_Sum_Rear));
    [maxvalue maxindex]=max(Signal_Env);
    Max_Position_Shift(q)=Position_micron(maxindex);
    FWHM(q)=abs(Position_micron(find(Signal_Env>0.5,1,'first'))-Position_micron(find(Signal_Env>0.5,1,'last')));

    [maxvalue maxindex]=max(Signal_Env_Rear);
    Max_Position_Shift_Rear(q)=Position_micron(maxindex);
end

subplot(1,2,1);
plot(Spot_Size,Max_Position_Shift,Spot_Size,Max_Position_Shift_Rear);
xlabel('Focus Spot Size (micron)','fontsize',12);
ylabel('PSF Shift (micron)','fontsize',12);
xlim([0 max(Spot_Size)]);
xlim([0 2000]);
title('PSF Shift Due to Finite Spot Size','fontsize',12);

legend('Front Interface','Rear Interface');
subplot(1,2,2);
plot(Spot_Size,Max_Position_Shift_Rear-Max_Position_Shift);
xlabel('Focus Spot Size (micron)','fontsize',12);
ylabel('PSF FWHM (micron)','fontsize',12);
xlim([0 max(Spot_Size)]);
xlim([0 2000]);
title('CCT Measurement Deviation Due to Finite Spot Size','fontsize',12);

dlmwrite('Spot_Size.txt',Spot_Size,'delimiter','\t','newline','pc');

dlmwrite('Max_Position_Shift_Front.txt',Max_Position_Shift,'delimiter','\t','newline','pc');

dlmwrite('Max_Position_Shift_Rear.txt',Max_Position_Shift_Rear,'delimiter','\t','newline','pc');

dlmwrite('Max_Position_Shift_Difference.txt',Max_Position_Shift_Rear-Max_Position_Shift,'delimiter','\t','newline','pc');

dlmwrite('FWHM.txt',FWHM,'delimiter','\t','newline','pc');

plot(Spot_Size,FWHM);
xlabel('Focus Spot Size (micron)','fontsize',12);
ylabel('FWHM (micron)','fontsize',12);
xlim([0 2000]);
title('PSF Broadened Due to Finite Spot Size','fontsize',16);

plot(Position_micron(abs(Position_micron)<100),abs(Signal(abs(Position_micron)<100)),Position_micron(abs(Position_micron)<100),Signal_Env(abs(Position_micron)<100));
xlabel('OPD (micron)','fontsize',12);
ylabel('Signal Intensity (micron)','fontsize',12);
legend('Spot Size < 10 micron','Spot Size = 1 mm');

plot(Position_micron(abs(Position_micron)<500),Signal_Env(abs(Position_micron)<500));
xlabel('OPD (micron)','fontsize',12);
ylabel('Signal Intensity (micron)','fontsize',12);
