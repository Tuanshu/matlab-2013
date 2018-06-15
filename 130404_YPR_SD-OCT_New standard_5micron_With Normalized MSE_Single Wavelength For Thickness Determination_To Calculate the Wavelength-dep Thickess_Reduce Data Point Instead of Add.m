%% Options

clear all
index_calculation=1;


Max_Number_of_Loop=3000;

Merit_A_OK=0.0001;         %mean(delta_A./A)
Merit_B_OK=0.0001;         %mean(delta_B./B)
Merit_C_OK=0.001;         %mean(delta_C./C)

Wavelength_Center_for_Thickness_Weighting=540;
Wavelength_BW_for_Thickness_Weighting=0.01;

Only_Calc_The_First=0;

Thickness_Determination_Wavelength=0.540; %micron

Sample_Path='D:\Users\TuanShu\130220_YPR\5micron\';
Reference_Path='D:\Users\TuanShu\130220_YPR\glass\';
Spectroscopy_Path='D:\Users\TuanShu\';

cd(sprintf('%s\\',Spectroscopy_Path));
Data_Spectroscopy=importdata('1208030331_5micron 1.jws.txt');
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


Ratio_Upper2Reference=1.56;
Thickness=(4.6).*1E-6;      %Initial Value
n_should=1.6;
Ratio_Lower2Upper=1.06;

%%
DC_cutoff=13;   %works for all lateral position
Least_Separation=3;
Window_Size=3.6;  %use for the left side of 1st interface and right side of second interface
No_Second_Interface_Threshold=0.1;

%%
Wavelength_Center=540;

Center_Wavelength_micron=0.56;
Wavelength_Considered_Min=539.2;          %nm
Wavelength_Considered_Max=540;

Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=N_f*8;


%% Global arrays generation

c=3E8;

Max_Frequency=c/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=c/(Max_Wavelength*1E-9);             %Hz

Frequency_Center=c/(Wavelength_Center*1E-9);
Frequency_Considered_Min=c/(Wavelength_Considered_Max*1E-9);             %Hz
Frequency_Considered_Max=c/(Wavelength_Considered_Min*1E-9);             %Hz


cd(Sample_Path);
Data=importdata('D0.txt');      % Data_2: the glass data 1

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

%% T 



Spectrum_Spectroscopy_Old=Data_Spectroscopy(end:-1:1,2)/100;
Wavelength_Spectroscopy=Data_Spectroscopy(end:-1:1,1);     
Frequency_Spectroscopy=c./(Wavelength_Spectroscopy*1E-9);

T=interp1(Frequency_Spectroscopy,Spectrum_Spectroscopy_Old,Frequency,'spline'); 

%% Theory - Sample Model (n1 - n - n2)

n2=1;

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
n1=n_AN100;

%% To generate the original spectrum

r_AN100=((1-n_AN100)./(n_AN100+1));
t_AN100=(2*(n_AN100)./(n_AN100+1));

%% Data Loading
cd(Sample_Path);
Spectrum_Sample_Frequency(1:length(Wavelength),1:fix(length(array)/Lateral_Downsample_Ratio))=0;
for array_number=1:fix(length(array)/Lateral_Downsample_Ratio)
    Data_Temp=importdata(sprintf('D%i.txt',array(1+(array_number-1)*Lateral_Downsample_Ratio)));      % Data_2: the glass data 1
    Spectrum_Sample_Frequency(:,array_number)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
    fprintf('Loading B-Scan Data...%d/%d\n',array_number,length(array));
end

cd(Reference_Path);
Spectrum_Reference_Frequency(1:length(Wavelength),1:fix(length(array)/Lateral_Downsample_Ratio))=0;
for array_number=1:fix(length(array)/Lateral_Downsample_Ratio)
    Data_Temp=importdata(sprintf('D%i.txt',array(1+(array_number-1)*Lateral_Downsample_Ratio)));      % Data_2: the glass data 1
    Spectrum_Reference_Frequency(:,array_number)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
    fprintf('Loading B-Scan Data...%d/%d\n',array_number,length(array));
end

%plot(Wavelength,Spectrum_Sample_Wavelength);
%xlabel('Wavelength (micron)');
%ylabel('Spectral Power (a.u.)');

%%

Spectrum_Sample=interp1(Frequency_Old,Spectrum_Sample_Frequency,Frequency); 
Spectrum_Sample(isnan(Spectrum_Sample))=0;
Spectrum_Sample(Frequency<Min_Frequency,:)=0;
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

plot(Position_micron,Signal_Sample,Position_micron,Signal_Reference);


%% Spectrum generation
DC_cutoff_index=find(Position_micron>DC_cutoff,1,'first');
Least_Separation_index=find(Position_micron>Least_Separation,1,'first');
Window_Size_index=find(Position_micron>Window_Size,1,'first');
Signal_Sample(1:DC_cutoff_index,:)=0;
Signal_Reference(1:DC_cutoff_index,:)=0;
[Max_Value Max_index]=max(abs(Signal_Sample));
[Max_Value_Ref Max_index_Ref]=max(abs(Signal_Reference));
Signal_Sample_First=Signal_Sample;
Signal_Sample_Second=Signal_Sample;
clear Signal_Sample
for p=1:length(Max_index)
    Signal_Sample_Second((Max_index(p)-Least_Separation_index):(Max_index(p)+Least_Separation_index),p)=0;
end

[Max_Value_Second Max_index_Second]=max(abs(Signal_Sample_Second));
If_Second_Exist=Max_Value_Second>Max_Value*No_Second_Interface_Threshold;
If_Second_Lower=Max_index_Second>Max_index;                                    %1 if second is lower, 0 if second is upper

Separation_index=round((Max_index+Max_index_Second)/2);

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
        Signal_Sample_First(1:(Max_index(p)-Window_Size_index),p)=0;
        Signal_Sample_First((Max_index(p)+Window_Size_index):end,p)=0;
        Signal_Sample_Second(:,p)=0;
        Max_index_Second(p)=Max_index(p); %to avoid wrong calc of thickness (OPD)
        Signal_Sample_1(:,p)=Signal_Sample_First(:,p);
        
    elseif If_Second_Lower(p)==0 %Second exist and second is upper 
        Signal_Sample_First(1:max(Separation_index(p),Max_index(p)-Window_Size_index),p)=0;
        Signal_Sample_First((Max_index(p)+Window_Size_index):end,p)=0;
        Signal_Sample_Second(1:(Max_index_Second(p)-Window_Size_index),p)=0;
        Signal_Sample_Second(min(Separation_index(p),Max_index_Second(p)+Window_Size_index):end,p)=0;
        Signal_Sample_1(:,p)=Signal_Sample_Second(:,p);
        Signal_Sample_2(:,p)=Signal_Sample_First(:,p);

    else %Second exist and second is lower        
        Signal_Sample_First(min(Separation_index(p),Max_index(p)+Window_Size_index):end,p)=0;
        Signal_Sample_First(1:(Max_index(p)-Window_Size_index),p)=0;
        Signal_Sample_Second((Max_index_Second(p)+Window_Size_index):end,p)=0;
        Signal_Sample_Second(1:max(Separation_index(p),(Max_index_Second(p)-Window_Size_index)),p)=0;
        
        Signal_Sample_1(:,p)=Signal_Sample_First(:,p);
        Signal_Sample_2(:,p)=Signal_Sample_Second(:,p);
        
    end
    Signal_Reference(1:(Max_index_Ref(p)-Window_Size_index),p)=0;
    Signal_Reference((Max_index_Ref(p)+Window_Size_index):end,p)=0;
end

plot(Position_micron(Position_micron<50),Signal_Sample_1(Position_micron<50),Position_micron(Position_micron<50),Signal_Sample_2(Position_micron<50),Position_micron(Position_micron<50),Signal_Reference(Position_micron<50));
xlabel('Optical Path Difference (micron)');
ylabel('Amplitude (a.u.)');

clear Signal_Sample_First Signal_Sample_Second


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

plot(Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Devided(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Devided_2(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),(abs(Spectrum_Devided_2(Wavelength_micron<1)./abs(Spectrum_Devided(Wavelength_micron<1)))));
xlabel('Wavelength (micron)');
ylabel('Interference Power spectral density (a.u.)');

plot(Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Reference(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Sample_1(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Sample_2(Wavelength_micron<1)));
xlabel('Wavelength (micron)');
ylabel('Interference Power spectral density (a.u.)');

Thickness_Profile=abs((Max_index_Second-Max_index)*(Position_micron(2)-Position_micron(1)));

plot((1:size(Thickness_Profile,2))*5,Thickness_Profile(1,:));
xlabel('Lateral Position (micron)');
ylabel('Thickness (OPD, micron)');

%%


Frequency_Considered=repmat(Frequency(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
Frequency_Considered=downsample(Frequency_Considered,Frequency_Downsample_Ratio);
Aexp=abs(Spectrum_Devided(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,:));  
Aexp=downsample(Aexp,Frequency_Downsample_Ratio);
Bexp=abs(Spectrum_Devided_2(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,:));  
Bexp=downsample(Bexp,Frequency_Downsample_Ratio);
Cexp=angle(Spectrum_Devided_3(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,:)); 
Cexp=downsample(Cexp,Frequency_Downsample_Ratio);

clear Spectrum_Devided Spectrum_Devided_2 Spectrum_Devided_3

T_Considered=repmat(T(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
T_Considered=downsample(T_Considered,Frequency_Downsample_Ratio);
%Weight_Function=(abs(Spectrum_Reference(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index)).*T_Considered)./max((abs(Spectrum_Reference(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index).*T_Considered)));
Wavelength_micron_Considered=Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
Wavelength_micron_Considered=downsample(Wavelength_micron_Considered,Frequency_Downsample_Ratio);
Center_Wavelength_micron_Index=find(Wavelength_micron_Considered<Center_Wavelength_micron,1,'first');  
n1_Considered=repmat(n1(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
n1_Considered=downsample(n1_Considered,Frequency_Downsample_Ratio);

n2_Considered=repmat(n2,size(n1_Considered,1),Active_Array_Size);


r_AN100_Considered=repmat(r_AN100(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
r_AN100_Considered=downsample(r_AN100_Considered,Frequency_Downsample_Ratio);
t_AN100_Considered=repmat(t_AN100(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
t_AN100_Considered=downsample(t_AN100_Considered,Frequency_Downsample_Ratio);

Cexp=unwrap(Cexp,[],1);

Cexp=diff(Cexp,[],1);

if size(Cexp,1)>1
    Cexp(size(Cexp,1)+1,:)=2*Cexp(size(Cexp,1),:)-Cexp(size(Cexp,1)-1,:);
else
    Cexp(size(Cexp,1)+1,:)=Cexp(size(Cexp,1),:);
end

Weighting_Function_for_Thickness=gaussmf(Frequency_Considered,[Frequency_BW_for_Thickness_Weighting Frequency_Center_for_Thickness_Weighting]);
Weighting_Function_for_Thickness=Weighting_Function_for_Thickness/sum(Weighting_Function_for_Thickness);

plot(Wavelength_micron_Considered,Aexp,Wavelength_micron_Considered,Bexp);
xlabel('Wavelength (micron)');
ylabel('Normalized Inteference Spectrum');
%for qq=1:Active_Array_Size            no more phase shifting
%    Phase_diff=Cexp(Center_Wavelength_micron_Index,qq);
%    NN=fix(Phase_diff/(2*pi));
%    Cexp(:,qq)=Cexp(:,qq)-NN*2*pi;
    %while Phase_diff > pi
    %    if Phase_diff>pi
    %        Cexp(:,qq)=Cexp(:,qq)-2*pi;
    %    elseif Phase_diff<pi
    %        Cexp(:,qq)=Cexp(:,qq)+2*pi;
    %    end
    %    Phase_diff=Cexp(Center_Wavelength_micron_Index,qq);
    %end
%end

plot(Wavelength_micron_Considered,Cexp);
xlabel('Wavelength (micron)');
ylabel('Continuous Phase Spectrum');


Thickness_Profile=repmat(Thickness_Profile(1:Active_Array_Size),length(Frequency_Considered),1);
%%

t1= @ (n) (2*n1_Considered)./(n1_Considered+n); 
t1_r= @ (n) (2*n)./(n2_Considered+n);
t2= @ (n) (2*n)./(n+n2_Considered); 
r1= @ (n) (n1_Considered-n)./(n1_Considered+n);
r1_r= @ (n) (n-n1_Considered)./(n+n1_Considered);
r2= @ (n) (n-n2_Considered)./(n+n2_Considered);

% n得比A, B, C要長1 point
% 但dn會跟A, B, C一樣長

A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) abs(Ratio_Upper2Reference_Now./r_AN100_Considered.*(r1(n)));
B = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) abs(Ratio_Upper2Reference_Now./r_AN100_Considered.*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(1i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c))));
C = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) unwrap(angle(1./r1(n).*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(1i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c)))),[],1); %雖然沒有絕對的phase, 但是是不是應該只shift pi的整數?
                   
det_array = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9) Q1.*Q5.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 -Q3.*Q5.*Q7 -Q2.*Q4.*Q9 -Q1.*Q6.*Q8;
 
T_abs = @ (n,Thickness_Now) abs(t_AN100_Considered.*t1(n).*t2(n).*exp(1i*2*pi*Frequency_Considered.*n.*Thickness_Now/c)./(1-r1_r(n).*r2(n).*exp(1i*4*pi*Frequency_Considered.*n.*Thickness_Now/c))).^2;          %Note, 120827 remove 2 from 分子4*pi > 2*pi
                    
dn=1E-10;
dt=1E-16;

n(1:length(Frequency_Considered),1:Active_Array_Size)=n_should;
Thickness_Now(1:length(Frequency_Considered),1:Active_Array_Size)=Thickness*Thickness_Profile/max(max(Thickness_Profile));
%plot(Wavelength_micron_Considered,Aexp,Wavelength_micron_Considered,A(n,Thickness,Ratio_Lower2Upper(1),Ratio_Upper2Reference(1)));
%plot(Wavelength_micron_Considered,Bexp,Wavelength_micron_Considered,B(n,Th
%ickness,Ratio_Lower2Upper(1),Ratio_Upper2Reference(1)));
%plot(Wavelength_micron_Considered,Cexp,Wavelength_micron_Considered,C(n,35E-6,Ratio_Lower2Upper(1),Ratio_Upper2Reference(1)));


dn_real_record=0;
dn_imag_record=0;
dthickness_record=0;
%% Start the n k fitting    
Merit_A=99999999;
Merit_B=99999999;
Merit_C=99999999;
Current_Progress=0;
if index_calculation==1   
    Merit_Best=999999999999999999999;
    Ratio_Lower2Upper_Best=1;
    Ratio_Upper2Reference_Best=1;
    for p=1:length(Thickness)
        for q=1:length(Ratio_Lower2Upper)
            for w=1:length(Ratio_Upper2Reference)
                    

                    %Thickness_Now=Thickness(p);
                    Ratio_Lower2Upper_Now=Ratio_Lower2Upper(q);
                    Ratio_Upper2Reference_Now=Ratio_Upper2Reference(w);
                    n(1:length(Frequency_Considered),1:Active_Array_Size)=n_should;
                    Thickness_Now(1:length(Frequency_Considered),1:Active_Array_Size)=Thickness*Thickness_Profile/max(max(Thickness_Profile));
                    Current_Loop=1;
                    qqqq=1;
                    while (Current_Loop<Max_Number_of_Loop) && ((Merit_A>Merit_A_OK)||(Merit_B>Merit_B_OK)||(Merit_C>Merit_C_OK))

                        C_Temp=C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now);
                        C_Temp=diff(C_Temp,[],1);
                        %if size(C_Temp,1)>1
                        %    C_Temp(size(C_Temp,1)+1,:)=2*C_Temp(size(C_Temp,1),:)-C_Temp(size(C_Temp,1)-1,:);
                        %else
                        %    C_Temp(size(C_Temp,1)+1,:)=C_Temp(size(C_Temp,1),:);
                        %end
                        
                        dA_dn_real = (A(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dB_dn_real = (B(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dC_dn_real = (C(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                    
                        dA_dn_imag = (A(n+1i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dB_dn_imag = (B(n+1i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dC_dn_imag = (C(n+1i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        
                        
                        dA_dthickness = (A(n,Thickness_Now+dt,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dt;
                        dB_dthickness = (B(n,Thickness_Now+dt,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dt;
                        dC_dthickness = (C(n,Thickness_Now+dt,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dt;
                        
                        delta_A=Aexp-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now);
                        delta_B=Bexp-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now);
                        delta_C=Cexp-C_Temp;                        
                        MSE_A=(sum((delta_A./Aexp).^2)/length(delta_A)).^0.5;
                        MSE_B=(sum((delta_B./Bexp).^2)/length(delta_B)).^0.5;
                        MSE_C=(sum((delta_C./Cexp).^2)/length(delta_C)).^0.5;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        dn_real = det_array(delta_A,dA_dn_imag,dA_dthickness,delta_B,dB_dn_imag,dB_dthickness,delta_C,dC_dn_imag,dC_dthickness)./det_array(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness);
                        dn_imag = det_array(dA_dn_real,delta_A,dA_dthickness,dB_dn_real,delta_B,dB_dthickness,dC_dn_real,delta_C,dC_dthickness)./det_array(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness);
                        dthickness = det_array(dA_dn_real,dA_dn_imag,delta_A,dB_dn_real,dB_dn_imag,delta_B,dC_dn_real,dC_dn_imag,delta_C)./det_array(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness);
                        %dn_real_record(length(dn_real_record)+1)=sum(dn_real.^2).^0.5;
                        %dn_imag_record(length(dn_real_record)+1)=sum(dn_imag.^2).^0.5;
                        %dthickness_record(length(dn_real_record)+1)=sum(dthickness.^2).^0.5;
                        n=n+0.1*dn_real+0.6*1i*dn_imag;
                        % (The thickness is different for different wavelength) Thickness_Now=Thickness_Now+0.5*dthickness;   %How if I make the thickness a array during the iteration?
                        %Thickness_Now=Thickness_Now+0.2*sum(Weighting_Function_for_Thickness.*dthickness);   %How if I make the thickness a array during the iteration?
                        %Thickness_Now=Thickness_Now+0.2*dthickness(find(Wavelength_micron_Considered<Thickness_Determination_Wavelength,1,'first'));
                        Thickness_Now=Thickness_Now+0.2*dthickness;
                        %Thickness_Now=Thickness_Now+0.2*dthickness;%(find(Wavelength_micron_Considered<Thickness_Determination_Wavelength,1,'first'));
                        Current_Loop=Current_Loop+1;
                        fprintf('%d/%d     %f%%\n',Current_Loop,Max_Number_of_Loop,Current_Progress);
                        fprintf('%d/%d\n',MSE_A,Merit_A_OK);
                        fprintf('%d/%d\n',MSE_B,Merit_B_OK);
                        fprintf('%d/%d\n',MSE_C,Merit_C_OK);

                        n(isnan(n))=n_should;
                        n(isinf(n))=n_should;
                        if (max(isnan(Thickness_Now)>0))
                            Thickness_Now(isnan(Thickness_Now))=Thickness*Thickness_Profile(isnan(Thickness_Now))/max(max(Thickness_Profile));
                        end
                        if (max(isinf(Thickness_Now)>0))
                            Thickness_Now(isinf(Thickness_Now))=Thickness*Thickness_Profile(isnan(Thickness_Now))/max(max(Thickness_Profile));
                        end

                        
                    end
                        Current_Progress=100*((p-1)+((q-1)+(w-1)/length(Ratio_Upper2Reference))/length(Ratio_Lower2Upper))/length(Thickness);
                        %fprintf('%f%%\n',Current_Progress);
                        Merit=sum((T_Considered-T_abs(n,Thickness_Now)).^2,1);
                        If_Merit_min=Merit<Merit_Best;
                        if max(If_Merit_min)
                            Merit_Best=Merit.*If_Merit_min+Merit_Best.*(1-If_Merit_min);
                        %dthickness_Best=dthickness;
                            Ratio_Lower2Upper_Best=Ratio_Lower2Upper_Now.*If_Merit_min+Ratio_Lower2Upper_Best.*(1-If_Merit_min);
                            Ratio_Upper2Reference_Best=Ratio_Upper2Reference_Now.*If_Merit_min+Ratio_Upper2Reference_Best.*(1-If_Merit_min);
                            T_Best=T_abs(n,Thickness_Now);
                            %C_Best=C_Temp;
                            n_Best=n;
                            
                            Thickness_Best=mean(Thickness_Now,1);
                            %Thickness_Best=Thickness_Now(Center_Wavelength_micron_Index,:);
                        end
                        %Thickness_Best_Record=(Thickness_Now);
                            

            end
        end
    end
end



subplot(3,1,1)
plot(Frequency_Considered,Aexp,Frequency_Considered,A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now));



subplot(3,1,2)
plot(Frequency_Considered,Bexp,Frequency_Considered,B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now));


subplot(3,1,3)
plot(Frequency_Considered,Cexp,Frequency_Considered,C_Temp);




plot(Wavelength_micron_Considered,real(n_Best),Wavelength_micron_Considered,imag(n_Best));

plot(Wavelength_micron_Considered,T_Best,Wavelength_micron_Considered,T_Considered);
%% Plot
subplot(3,1,1)
plot(Wavelength_micron_Considered,real(n_Best));%,Wavelength_micron_Considered,n1_Considered(:,1));
xlabel('Wavelength (micron)');
ylabel('n');
ylim([1.55 1.651]);
xlim([0.51 0.6]);


%legend('Calculated Film Refractive Index','Refractive Index of Substrate','Location','Best')

subplot(3,1,2)
plot(Wavelength_micron_Considered,imag(n_Best));
xlabel('Wavelength (micron)');
ylabel('k');
ylim([0.008 0.012]);
xlim([0.51 0.6]);

%%

subplot(3,1,3)
plot(Wavelength_micron_Considered,T_Best,Wavelength_micron_Considered,T_Considered(:,1));
xlabel('Wavelength (micron)');
ylabel('T');
ylim([0.5 1]);
xlim([0.51 0.6]);

%legend('Calculated Transmission','Measured Transmission','Location','Best')
%ylim([0.4 1]);
%set(gca,'YTick',0.6:0.1:1)

%subplot(6,1,4)
%plot(Wavelength_micron_Considered,Aexp,Wavelength_micron_Considered,A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now));
%xlabel('Wavelength (micron)');
%ylabel('A');

%subplot(6,1,5)
%plot(Wavelength_micron_Considered,Bexp,Wavelength_micron_Considered,B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now));
%xlabel('Wavelength (micron)');
%ylabel('B');

%subplot(6,1,6)
%plot(Wavelength_micron_Considered,Cexp,Wavelength_micron_Considered,C_Temp);
%xlabel('Wavelength (micron)');
%ylabel('C');

%subplot(6,1,4)
%plot(2:(length(dn_real_record)-1),dn_real_record(3:end));
%xlabel('Number');
%ylabel('dn_real_record');

%subplot(6,1,5)
%plot(2:(length(dn_imag_record)-1),dn_imag_record(3:end));
%xlabel('Number');
%ylabel('dn_imag_record');

%subplot(6,1,6)
%plot(2:(length(dthickness_record)-1),dthickness_record(3:end));
%xlabel('Number');
%ylabel('dthickness_record');

%plot((1:length(Thickness_Best))*5/1000,Thickness_Best*1E6,(1:length(Thickness_Best))*5/1000,Thickness*1E6*Thickness_Profile(1,:)/max(max(Thickness_Profile)));
%xlabel('Lateral Position (mm)');
%ylabel('Thickness (micron)');
%legend('Calculated Thickness','Predicted Thickness (Lieanr to OPD)')

%%

%plot(Wavelength_micron_Considered,mean(T_Best,2),Wavelength_micron_Considered,T_Considered(:,1));
%xlabel('Wavelength (micron)');
%ylabel('T');

%dlmwrite('Lateral Position.txt',((1:length(Thickness_Best))*5/1000)','delimiter','\t','newline','pc');
%dlmwrite('Thickness.txt',Thickness_Best'*1E6,'delimiter','\t','newline','pc');
%dlmwrite('OPD.txt',Thickness_Profile(1,:)','delimiter','\t','newline','pc');

%dlmwrite('T_Best.txt',mean(T_Best(:,1:100),2),'delimiter','\t','newline','pc');

%dlmwrite('T_Considered.txt',T_Considered(:,1),'delimiter','\t','newline','pc');

%dlmwrite('Wavelength_micron_Considered.txt',Wavelength_micron_Considered,'delimiter','\t','newline','pc');

%dlmwrite('n.txt',mean(real(n_Best(:,1:100)),2),'delimiter','\t','newline','pc');

%dlmwrite('k.txt',mean(imag(n_Best(:,1:100)),2),'delimiter','\t','newline','pc');

%plot(Wavelength_micron_Considered,A(n_Best,Thickness_Best,Ratio_Lower2Upper_Best,Ratio_Upper2Reference_Best),Wavelength_micron_Considered,Aexp);
%plot(Wavelength_micron_Considered,B(n_Best,Thickness_Best,Ratio_Lower2Upper_Best,Ratio_Upper2Reference_Best),Wavelength_micron_Considered,Bexp);
%plot(Wavelength_micron_Consi
mean(Thickness_Now)
MSE_A_Record=(sum((delta_A./Aexp).^2)/length(delta_A)).^0.5;
MSE_B_Record=(sum((delta_B./Bexp).^2)/length(delta_B)).^0.5;
MSE_C_Record=(sum((delta_C./Cexp).^2)/length(delta_C)).^0.5;
fprintf('MSE_A=%e, MSE_B=%e, MSE_C=%f\n',MSE_A_Record,MSE_B_Record,MSE_C_Record);

cd('D:\Users\TuanShu\');


dlmwrite('Wavelength_micron_Considered_5micron_NewAlgorithm.txt',Wavelength_micron_Considered,'delimiter','\t','newline','pc');
dlmwrite('n_5micron_NewAlgorithm.txt',real(n_Best),'delimiter','\t','newline','pc');
dlmwrite('k_5micron_NewAlgorithm.txt',imag(n_Best),'delimiter','\t','newline','pc');
dlmwrite('T_5micron_NewAlgorithm.txt',T_Best,'delimiter','\t','newline','pc');
dlmwrite('T_exp_5micron_NewAlgorithm.txt',T_Considered,'delimiter','\t','newline','pc');