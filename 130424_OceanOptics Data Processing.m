clear all

%% Data Loading

cd('D:\Users\TuanShu\');

Data=importdata('130424_Air_incident.txt');

%Data_sam=importdata('130424_Substrate_incident_with reference cover glass_3_sam.txt');
%Data(:,2)=Data(:,2)-Data_sam(:,2);
%% Setting

thickness_temp=2.68E-6;
n_should=1.5;

Starting_pixel_f_considered=250;
Ending_pixel_f_considered=350;

range_specified=1;

Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=8192*8;
ROI_ratio=1/4;                  %only consider the first ROI_ratio data in TD

DC_cutoff=35;                 %in TD

array=1;

%% Data Loading

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
%plot(Frequency_New,Spectrum_New);

%% To time domain

Spectrum_New((N_f+1):N_t)=0;



Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2�O�]���@�Ӥ@�^
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;

Time=Time(1:round(length(Time)*ROI_ratio));
Position=Position(1:round(length(Position)*ROI_ratio));
Position_micron=Position_micron(1:round(length(Position_micron)*ROI_ratio));


Signal=fft(Spectrum_New);

Spectrum_New=Spectrum_New(1:N_f);
Signal=Signal(1:round(length(Signal)*ROI_ratio));

Window=(gaussmf(Position_micron,[2 8]));
Window(Position_micron>8)=1;
%Signal=Signal.*Window;
%Signal(1:DC_cutoff)=0;
Signal_Carrier=real(Signal);
Signal_Envelope=abs(Signal);

%% ROI


%% Shift

Signal_Envelope_Cut=Signal_Envelope(Position_micron>DC_cutoff)/max(Signal_Envelope(Position_micron>DC_cutoff));
Position_micron_cut=Position_micron(Position_micron>DC_cutoff);
FWHM=Position_micron_cut(find(Signal_Envelope_Cut>0.5,1,'last'))-Position_micron_cut(find(Signal_Envelope_Cut>0.5,1,'first'));
FWHM_should=2*Position_micron(find((Signal_Envelope./max(Signal_Envelope))<0.5,1,'first'));


plot(Wavelength,Spectrum);
xlabel('Wavelength (nm)');
ylabel('Spectrum');

plot(Position_micron,Signal_Carrier);
xlabel('OPD (micron)');
ylabel('Interference Signal');

fprintf('FWHM = %f micron .\n',FWHM);
