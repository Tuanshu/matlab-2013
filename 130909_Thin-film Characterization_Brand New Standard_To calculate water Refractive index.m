%% Options

clear all
index_calculation=1;
load_nk=0;
incident=1;     %1: sub, 2: air
N_Shift=[74]';     %Physically must be integer, 若不是integer那會算錯
Max_Number_of_Loop=100;
%For single wavelength,
%A:1.989262e-04/1.000000e-04
%B:1.913883e-06/1.000000e-04
%C:1.643246e-10/1.000000e-04
%Hopefully this errorhe
%
Ratio_Lower2Upper=0.3;     %if the rear interface has larger interference eff than front interface, Ratio_Lower2Upper>1
Ratio_Upper2Reference=0.95;%0.835; %if the front interface has larger interference eff than referential glass, Ratio_Upper2Reference>1


MSE_A_OK=0.0001;         %mean(delta_A./A)
MSE_B_OK=0.0001;         %mean(delta_B./B)
MSE_C_OK=0.002;         %mean(delta_C./C)
MSE_total_OK=0.005;

Wavelength_Center_for_Thickness_Weighting=540;
Wavelength_BW_for_Thickness_Weighting=50;

Wavelength_Thickness_Detection=0.54;    %micron

Only_Calc_The_First=0;

Lambda=0;

Sample_Path='D:\Users\TuanShu\130830\3k+\sample 3\';
Reference_Path='D:\Users\TuanShu\130830\3k+\sample 3\';
Spectroscopy_Path='D:\Users\TuanShu\1208030331_5micron 1.jws.txt';
cd('D:\Users\TuanShu\');
%cd(sprintf('%s\\',Spectroscopy_Path));
Data_Spectroscopy=importdata(Spectroscopy_Path);
array=[0];

Frequency_Downsample_Ratio=5;
Lateral_Downsample_Ratio=1;


index_used_reference_plane=100;

lateral_index_calibration_start=50;
lateral_index_calibration_end=700;

lateral_index_reference_start=50;
lateral_index_reference_end=250;


lateral_index_sample_start=650;
lateral_index_sample_end=650;


Thickness=[16]'*1E-6;      %Initial Value
if load_nk ==1
    n_should=dlmread('n.txt')+1i*dlmread('k.txt');
elseif load_nk ==0
    n_should=1.33;
end

%%%
DC_cutoff=5;   %works for all lateral position
Least_Separation=5;
Window_Size_Right_1=3.8;  %use for the left side of 1st interface and right side of second interface
Window_Size_Left_1=4.6;  %use for the left side of 1st interface and right side of second interface
Window_Size_Right_2=5;  %use for the left side of 1st interface and right side of second interface
Window_Size_Left_2=5.1;  %use for the left side of 1st interface and right side of second interface
Reference_Window_Size=4;
No_Second_Interface_Threshold=0.05;

%%%
Wavelength_Center=540;

Center_Wavelength_micron=0.54;
Wavelength_Considered_Min=510; %510         %nm
Wavelength_Considered_Max=590;  %580

Max_Wavelength=700;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=N_f*16;


%%% Global arrays generation

c=3E8;

Max_Frequency=c/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=c/(Max_Wavelength*1E-9);             %Hz

Frequency_Center=c/(Wavelength_Center*1E-9);
Frequency_Considered_Min=c/(Wavelength_Considered_Max*1E-9);             %Hz
Frequency_Considered_Max=c/(Wavelength_Considered_Min*1E-9);             %Hz


if length(array)>1
    cd(Sample_Path);
    Data=importdata('D0.txt');      % Data_2: the glass data 1
elseif exist(Sample_Path,'file')  == 7
    cd(Sample_Path);
    Data_inter=importdata('inter.txt');      % Data_2: the glass data 1
    Data_Sam=importdata('sam.txt');      % Data_2: the glass data 1
    Data_Ref=importdata('ref.txt');      % Data_2: the glass data 1
    Data(:,2)=Data_inter(:,2)-Data_Sam(:,2)-Data_Ref(:,2);
    Data(:,1)=Data_inter(:,1);
else
    Data=importdata(Sample_Path);      % Data_2: the glass data 1
end

Wavelength=Data(:,1);           %nm
Frequency_Old=c./(Wavelength*1E-9);
Frequency=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency=Frequency';
Wavelength_micron=(c./Frequency)*1E6;

Frequency_Center_Index=find(Frequency>Frequency_Center,1,'first');
Frequency_Considered_Min_Index=find(Frequency>Frequency_Considered_Min,1,'first');
Frequency_Considered_Max_Index=find(Frequency>Frequency_Considered_Max,1,'first');

Frequency_Center_for_Thickness_Weighting=c./(Wavelength_Center_for_Thickness_Weighting*1E-9);
Frequency_BW_for_Thickness_Weighting=c.*(Wavelength_BW_for_Thickness_Weighting*1E-9)./(Wavelength_Center_for_Thickness_Weighting*1E-9).^2;

% Time-domain

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=c*Time;
Position_micron=Position*1E6;

%%% T 



Spectrum_Spectroscopy_Old=Data_Spectroscopy(end:-1:1,2)/100;
Wavelength_Spectroscopy=Data_Spectroscopy(end:-1:1,1);     
Frequency_Spectroscopy=c./(Wavelength_Spectroscopy*1E-9);

T=interp1(Frequency_Spectroscopy,Spectrum_Spectroscopy_Old,Frequency,'spline'); 

%%% Theory - Sample Model (n1 - n - n2)


% n1 = AN100
C1 = 1.03961212; 
C2 = 0.00600069867; 
C3 = 0.231792344; 
C4 = 0.0200179144; 
C5 = 1.01046945; 
C6 = 103.560653;


n_AN100=(C1*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C2)+C3*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C4)+C5*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C6)+1).^0.5;

n_AN100=abs(n_AN100);

n_AN100(isnan(n_AN100))=0;

if incident==1
    t_AN100=(2*(1)./(n_AN100+1));
    n1=n_AN100;
    n2=n_AN100;
elseif incident==2
    n1=n_AN100;
    n2=n_AN100;
    t_AN100=(2*(n_AN100)./(n_AN100+1));
end


%%% To generate the original spectrum

r_AN100=((1-n_AN100)./(n_AN100+1));

%%% Data Loading
Spectrum_Sample_Frequency(1:length(Wavelength),1:fix(length(array)/Lateral_Downsample_Ratio))=0;

if length(array)>1
    cd(Sample_Path);
    for array_number=1:fix(length(array)/Lateral_Downsample_Ratio)
        Data_Temp=importdata(sprintf('D%i.txt',array(1+(array_number-1)*Lateral_Downsample_Ratio)));      % Data_2: the glass data 1
        Spectrum_Sample_Frequency(:,array_number)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
        fprintf('Loading B-Scan Data...%d/%d\n',array_number,length(array));
    end
elseif exist(Sample_Path,'file') == 7
    cd(Sample_Path);
    Data_inter=importdata('inter.txt');      % Data_2: the glass data 1
    Data_Sam=importdata('sam.txt');      % Data_2: the glass data 1
    Data_Ref=importdata('ref.txt');      % Data_2: the glass data 1
    Data_Temp(:,2)=Data_inter(:,2)-Data_Sam(:,2)-Data_Ref(:,2);
    Data_Temp(:,1)=Data_inter(:,1);
    Spectrum_Sample_Frequency(:,1)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
else

    Data_Temp=importdata(Sample_Path);      % Data_2: the glass data 1
    Spectrum_Sample_Frequency(:,1)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
end

Spectrum_Reference_Frequency(1:length(Wavelength),1:fix(length(array)/Lateral_Downsample_Ratio))=0;
if length(array)>1
    cd(Reference_Path);
    for array_number=1:fix(length(array)/Lateral_Downsample_Ratio)
        Data_Temp=importdata(sprintf('D%i.txt',array(1+(array_number-1)*Lateral_Downsample_Ratio)));      % Data_2: the glass data 1
        Spectrum_Reference_Frequency(:,array_number)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
        fprintf('Loading B-Scan Data...%d/%d\n',array_number,length(array));
    end
elseif exist(Sample_Path,'file') == 7
    cd(Reference_Path);
    Data_inter=importdata('inter.txt');      % Data_2: the glass data 1
    Data_Sam=importdata('sam.txt');      % Data_2: the glass data 1
    Data_Ref=importdata('ref.txt');      % Data_2: the glass data 1
    Data_Temp(:,2)=Data_inter(:,2)-Data_Sam(:,2)-Data_Ref(:,2);
    Data_Temp(:,1)=Data_inter(:,1);
    Spectrum_Reference_Frequency(:,1)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
else
    Data_Temp=importdata(Reference_Path);
    Spectrum_Reference_Frequency(:,1)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
end


%%%

Spectrum_Sample=interp1(Frequency_Old,Spectrum_Sample_Frequency,Frequency); 
Spectrum_Sample(isnan(Spectrum_Sample))=0;
Spectrum_Sample(Frequency<Min_Frequency,:)=0;


plot(Wavelength_micron(Wavelength_micron<1),Spectrum_Sample(Wavelength_micron<1));
%xlabel('Wavelength (micron)');
%ylabel('Spectral Power (a.u.)');
Spectrum_Sample(N_f+1:N_t,:)=0;
Signal_Sample=fft(Spectrum_Sample,[],1);  
Signal_Sample(round(size(Signal_Sample,1)/2+1):end,:)=0;
Spectrum_Sample=Spectrum_Sample(1:N_f,:);


Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency); 
Spectrum_Reference(isnan(Spectrum_Reference))=0;
Spectrum_Reference(Frequency<Min_Frequency,:)=0;
Spectrum_Reference(N_f+1:N_t,:)=0;
Signal_Reference=fft(Spectrum_Reference,[],1);  
Signal_Reference(round(size(Signal_Reference,1)/2+1):end,:)=0;
Spectrum_Reference=Spectrum_Reference(1:N_f,:);


clear Spectrum_Sample_Frequency Spectrum_Reference_Frequency

plot(Wavelength_micron(Wavelength_micron<1),Spectrum_Sample(Wavelength_micron<1));
xlabel('Wavelength (micron)');
ylabel('Spectral Intesity (a.u.)');
xlim([0.45 0.7]);

plot(Position_micron,Signal_Sample);%,Position_micron,Signal_Reference);
xlabel('Optical Path Difference (micron)');
ylabel('Amplitude (a.u.)');
xlim([0 50]);

%%% Spectrum generation
DC_cutoff_index=find(Position_micron>DC_cutoff,1,'first');
Least_Separation_index=find(Position_micron>Least_Separation,1,'first');
Window_Size_Right_1_index=find(Position_micron>Window_Size_Right_1,1,'first');
Window_Size_Left_1_index=find(Position_micron>Window_Size_Left_1,1,'first');
Window_Size_Right_2_index=find(Position_micron>Window_Size_Right_2,1,'first');
Window_Size_Left_2_index=find(Position_micron>Window_Size_Left_2,1,'first');
Reference_Window_Size_index=find(Position_micron>Reference_Window_Size,1,'first');
Signal_Sample(1:DC_cutoff_index,:)=0;
Signal_Reference(1:DC_cutoff_index,:)=0;
[Max_Value Max_index]=max(abs(Signal_Sample));
[Max_Value_Ref Max_index_Ref]=max(abs(Signal_Reference));
Signal_Sample_First=Signal_Sample;
Signal_Sample_Second=Signal_Sample;

for p=1:length(Max_index)
    Signal_Sample_Second((Max_index(p)-Least_Separation_index):(Max_index(p)+Least_Separation_index),p)=0;
end

[Max_Value_Second Max_index_Second]=max(abs(Signal_Sample_Second));
If_Second_Exist=Max_Value_Second>Max_Value*No_Second_Interface_Threshold;
If_Second_Lower=Max_index_Second>Max_index;                                    %1 if second is lower, 0 if second is upper

if If_Second_Lower ==1
    [minvalue minindex]=min(abs(Signal_Sample_Second((Max_index+Least_Separation_index+1):Max_index_Second)));
    Separation_index=minindex+Max_index+Least_Separation_index;
elseif If_Second_Lower ==0
    [minvalue minindex]=min(abs(Signal_Sample_Second(Max_index_Second:(Max_index-Least_Separation_index-1))));
    Separation_index=minindex+Max_index_Second-1;
end
    
round((Max_index+Max_index_Second)/2);

First_No_Second_Index=find(If_Second_Exist==0,1,'first');
Last_With_Second_Index=find(If_Second_Exist==1,1,'last');

Safe_Reference_Index=round((First_No_Second_Index+length(If_Second_Exist))/2);
Safe_Sample_Index=round((Last_With_Second_Index+1)/2);

if (Last_With_Second_Index>First_No_Second_Index)
    disp('Second Peak Condition Working Badly!');
else
    fprintf('Film caoted from index: 1 to %d.\n',Last_With_Second_Index);
    if xor(abs(Max_index(Safe_Reference_Index)-Max_index(Safe_Sample_Index))>abs(Max_index(Safe_Reference_Index)-Max_index_Second(Safe_Sample_Index)),mean(If_Second_Lower(1:Safe_Sample_Index))>0.5)
        disp('Film facing down.');
    else
        disp('Film facing up.');
    end
        
end
if Only_Calc_The_First == 1
    Active_Array_Size=1;
else
    Active_Array_Size=length(array(1:Last_With_Second_Index));
end

Signal_Sample_1(1:N_t,1:length(array))=0;
Signal_Sample_2(1:N_t,1:length(array))=0;

for p=1:length(Max_index)
    if If_Second_Exist(p)==0 %Second non-exist
        Signal_Sample_First(1:(Max_index(p)-Window_Size_Left_1_index),p)=0;
        Signal_Sample_First((Max_index(p)+Window_Size_Right_1_index):end,p)=0;
        Signal_Sample_Second(:,p)=0;
        Max_index_Second(p)=Max_index(p); %to avoid wrong calc of thickness (OPD)
        Signal_Sample_1(:,p)=Signal_Sample_First(:,p);
        
    elseif If_Second_Lower(p)==0 %Second exist and second is upper 
        Signal_Sample_First(1:max(Separation_index(p),Max_index(p)-Window_Size_Left_2_index),p)=0;
        %Signal_Sample_First(1:Separation_index(p),p)=0;
        Signal_Sample_First((Max_index(p)+Window_Size_Right_2_index):end,p)=0;
        Signal_Sample_Second(1:(Max_index_Second(p)-Window_Size_Left_1_index),p)=0;
        Signal_Sample_Second(min(Separation_index(p),Max_index_Second(p)+Window_Size_Right_1_index):end,p)=0;
        %Signal_Sample_Second(Separation_index(p):end,p)=0;
        Signal_Sample_1(:,p)=Signal_Sample_Second(:,p);
        Signal_Sample_2(:,p)=Signal_Sample_First(:,p);

    else %Second exist and second is lower        
        Signal_Sample_First(min(Separation_index(p),Max_index(p)+Window_Size_Right_1_index):end,p)=0;
        %Signal_Sample_First(Separation_index(p):end,p)=0;
        Signal_Sample_First(1:(Max_index(p)-Window_Size_Left_1_index),p)=0;
        Signal_Sample_Second((Max_index_Second(p)+Window_Size_Right_2_index):end,p)=0;
        Signal_Sample_Second(1:max(Separation_index(p),(Max_index_Second(p)-Window_Size_Left_2_index)),p)=0;
        %Signal_Sample_Second(1:Separation_index(p),p)=0;
        
        Signal_Sample_1(:,p)=Signal_Sample_First(:,p);
        Signal_Sample_2(:,p)=Signal_Sample_Second(:,p);
        
    end
    Signal_Reference(1:(Max_index_Ref(p)-Reference_Window_Size_index),p)=0;
    Signal_Reference((Max_index_Ref(p)+Reference_Window_Size_index):end,p)=0;
end

clear Signal_Sample
%clear Signal_Sample_First Signal_Sample_Second


Spectrum_Sample_1=ifft(Signal_Sample_1,[],1);
Spectrum_Sample_1=2*Spectrum_Sample_1(1:N_f,:);                

Spectrum_Sample_2=ifft(Signal_Sample_2,[],1);
Spectrum_Sample_2=2*Spectrum_Sample_2(1:N_f,:);


Spectrum_Reference=ifft(Signal_Reference,[],1);
Spectrum_Reference=2*Spectrum_Reference(1:N_f,:);                



plot(real(Signal_Sample_1));
 

Spectrum_Sample_1=Spectrum_Sample_1(:,1:Active_Array_Size);
Spectrum_Sample_2=Spectrum_Sample_2(:,1:Active_Array_Size);   %!!!!!!!!!! 因為接下來就要算n k 了, 把沒有兩個介面的放進來情況會很怪

Spectrum_Devided=Spectrum_Sample_1./repmat(Spectrum_Reference,1,Active_Array_Size);
Spectrum_Devided_2=Spectrum_Sample_2./repmat(Spectrum_Reference,1,Active_Array_Size);
Spectrum_Devided_3=Spectrum_Sample_2./Spectrum_Sample_1;

subplot(1,3,1);

plot(Position_micron(Position_micron<100),Signal_Sample_1(Position_micron<100),Position_micron(Position_micron<100),Signal_Sample_2(Position_micron<100));%,Position_micron(Position_micron<100),Signal_Reference(Position_micron<100));
xlabel('Optical Path Difference (micron)');
ylabel('Amplitude (a.u.)');
xlim([0 40]);

%plot(Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Devided(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Devided_2(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),(abs(Spectrum_Devided_2(Wavelength_micron<1)./abs(Spectrum_Devided(Wavelength_micron<1)))));
%xlabel('Wavelength (micron)');
%ylabel('Interference Power spectral density (a.u.)');
subplot(1,3,2);
plot(Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Sample_1(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Sample_2(Wavelength_micron<1)));%,Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Reference(Wavelength_micron<1)));
xlabel('Wavelength (micron)');
ylabel('Interference Power spectral density (a.u.)');
fprintf('Peak height ratio:%f\n',Max_Value/Max_Value_Second);
xlim([0.45 0.7]);

subplot(1,3,3);
plot(Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Sample_2(Wavelength_micron<1))./abs(Spectrum_Sample_1(Wavelength_micron<1)));%,Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Sample_2(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Reference(Wavelength_micron<1)));
xlabel('Wavelength (micron)');
ylabel('Interference Power spectral density (a.u.)');
fprintf('Peak height ratio:%f\n',Max_Value/Max_Value_Second);
xlim([0.5 0.65]);
Thickness_Profile=abs((Max_index_Second-Max_index)*(Position_micron(2)-Position_micron(1)));
