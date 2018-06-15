%% Options


clear all
index_calculation=1;
Number_of_Loop=75;
Sampling_Rate=1035.197/10;
Frequency_Downsample_Ratio=1;
PZT_Speed=0.895972828129950;         %micron/sec


Select_Position=[31 40];

% Parameters
Ratio_Upper2Reference=1.92;%+[-0.05:0.005:0.05];
Thickness=(4.7).*1E-6;% ([-0.1:0.1:0.1]+2.55).*1E-6; In this file, it is the initial guess of thickness
n_should=1.55;
Ratio_Lower2Upper=0.75;%+[-0.05:0.005:0.05];%+[-0.2:0.02:0.2];  %([-1:0.02:1]+10.82)*1E-6;%[0.019:0.001:0.024]*1E-6;%[0.02:0.0003:0.023]*1E-6;



TD_OCT_point_Spacing=PZT_Speed/Sampling_Rate;      %micron �]���ڨ����I�ƬO������|��, �G�n�����ӬO169.8/4800/4=0.00884375

Wavelength_micron_Filter_min=0.40;

Wavelength_micron_Filter_max=0.80;

Only_Calc_The_First=0;


Sample_Path='D:\Users\TuanShu\130117\130117_Sam\Bscan\';
Reference_Path='D:\Users\TuanShu\130117\130117_Ref\Bscan\';
Spectroscopy_Path='D:\Users\TuanShu\';

delta_n=0.1;
Wavelength_Center=540;

Center_Wavelength_micron=0.54;
Wavelength_Considered_Min=485;          %nm
Wavelength_Considered_Max=600;

Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm



% Data Loading

Considered_Lateral_Index_Range=1:60;
Considered_Axial_Index_Range=1:3000;
Number_of_Axial_point_After_Exapnding=100000;

%
Time=(1:(Number_of_Axial_point_After_Exapnding))/Sampling_Rate;

Distance_of_Largest_Signal_of_Sample_from_Zero=20;           %micron, a set value (for shifting), ��Sample��X�һ�shift���I�ƫ�, sample�Mreference shift�ۦP���q (for YPR, �����U����)

cd(Sample_Path);
Frame_Temp=dlmread(sprintf('Bscan_Sam_%i.txt',Select_Position(2))); 
Signal_Sam=Frame_Temp(:,:);
cd(Reference_Path);
Frame_Temp=dlmread(sprintf('Bscan_Ref_%i.txt',Select_Position(2))); 
Signal_Ref=Frame_Temp(:,:);


%Signal_Sam=[Signal_Sam Signal_Ref];                                         %%%%%%%%%%%%%%%%%%%%%%%%% NOTE!!!!!!! �Ȯɳo�˼g, ��ɭ�FF-OCT�i��|���Ӥ@��

Signal_Sam=Signal_Sam(Considered_Axial_Index_Range,Considered_Lateral_Index_Range);
Signal_Ref=Signal_Ref(Considered_Axial_Index_Range,Considered_Lateral_Index_Range);     %�Y��shift����I����, �|�y���b�W�ФW���@�ܤj��phase

[maxvalue_sam maxindex_sam]=max(Signal_Sam,[],1);

Shift_Amount_sam=1500-maxindex_sam(1);  %���ܪ����βĤ@���shift, circshift�Y1d array�n���|�ǩ�

Signal_Sam=circshift(Signal_Sam,Shift_Amount_sam);      %�L�u�����C��shift���P���q��? �n�����S��, ���M�Ӳ{�b���g�k, reference���Ӥ]�|�Qset��20 (�i�O�S��)

[maxvalue_ref maxindex_ref]=max(Signal_Sam,[],1);

Shift_Amount_ref=1500-maxindex_ref(1);  %���ܪ����βĤ@���shift, circshift�Y1d array�n���|�ǩ�

Signal_Ref=circshift(Signal_Ref,Shift_Amount_ref);      %�L�u�����C��shift���P���q��? �n�����S��, ���M�Ӳ{�b���g�k, reference���Ӥ]�|�Qset��20 (�i�O�S��)

Time=Time-Time(1000);

%Position_micron=(1:(Number_of_Axial_point_After_Exapnding))*TD_OCT_point_Spacing;
%Position_micron=Position_micron-max(Position_micron)/200;

Signal_Sam((length(Considered_Axial_Index_Range)+1):Number_of_Axial_point_After_Exapnding,:)=0;
Signal_Ref((length(Considered_Axial_Index_Range)+1):Number_of_Axial_point_After_Exapnding,:)=0;

Position_micron=(1:(Number_of_Axial_point_After_Exapnding))*TD_OCT_point_Spacing;

Spectrum_Sam=ifft(Signal_Sam,[],1);
Spectrum_Ref=ifft(Signal_Ref,[],1);

c=3E8;

df=c/((max(Position_micron))*(1E-6)*2); %the OPD is roundtrp
Frequency=df:df:df*(length(Position_micron));  
Frequency=Frequency';
Wavelength_micron=(c./Frequency)*1E6;

plot(Wavelength_micron(Wavelength_micron<1),real(Spectrum_Sam(Wavelength_micron<1,:)),Wavelength_micron(Wavelength_micron<1),real(Spectrum_Ref(Wavelength_micron<1,:)));           %�Y�O�W�иѪR�פӮt, �O���O�]�|�]��Signal

Frequency_Filter_max=Frequency(find(Wavelength_micron<Wavelength_micron_Filter_min,1,'first'));

Frequency_Filter_min=Frequency(find(Wavelength_micron<Wavelength_micron_Filter_max,1,'first'));

Spectrum_Sam(Frequency<Frequency_Filter_min,:)=0;

Spectrum_Sam(Frequency>Frequency_Filter_max,:)=0;

Spectrum_Ref(Frequency<Frequency_Filter_min,:)=0;

Spectrum_Ref(Frequency>Frequency_Filter_max,:)=0;

Signal_Sam_Filtered=fft(Spectrum_Sam,[],1);
Signal_Ref_Filtered=fft(Spectrum_Ref,[],1);


%plot(Wavelength_micron(Wavelength_micron<1),real(Spectrum_Sam(Wavelength_micron<1,:))); 

%
%plot(Position_micron(Position_micron<50),Signal_Sam_Filtered(Position_micron<50,1),Position_micron(Position_micron<50),Signal_Sam_Filtered(Position_micron<50,2));%,Position_micron,Signal_Ref_Filtered(:,:));



% Signal Separation
Least_Separation=3.5;
Window_Size=5;
No_Second_Interface_Threshold=0.05;
Least_Separation_index=find(Position_micron>Least_Separation,1,'first');
Window_Size_index=find(Position_micron>Window_Size,1,'first');
[Max_Value Max_index]=max(abs(Signal_Sam_Filtered));
[Max_Value_Ref Max_index_Ref]=max(abs(Signal_Ref_Filtered));
Signal_Sample_First=Signal_Sam_Filtered;
Signal_Sample_Second=Signal_Sam_Filtered;
Signal_Reference=Signal_Ref_Filtered;

for p=1:length(Max_index)
    Signal_Sample_Second((Max_index(p)-Least_Separation_index):(Max_index(p)+Least_Separation_index),p)=0;
end

[Max_Value_Second Max_index_Second]=max(abs(Signal_Sample_Second));
If_Second_Exist=Max_Value_Second>Max_Value*No_Second_Interface_Threshold;
If_Second_Lower=Max_index_Second>Max_index;                                    %1 if second is lower, 0 if second is upper

Separation_index=round((Max_index+Max_index_Second)/2);
%Separation_index(1)=find(Position_micron>10.7,1,'first');


First_With_Second_Index=find(If_Second_Exist==1,1,'first');
Last_No_Second_Index=find(If_Second_Exist==0,1,'last');

Safe_Reference_Index=1;
Safe_Sample_Index=round((First_With_Second_Index+length(Considered_Lateral_Index_Range))/2);

if (Last_No_Second_Index>First_With_Second_Index)
    disp('Second Peak Condition Working Badly!');
else
    fprintf('Film caoted from index: %d to End.\n',First_With_Second_Index);
    if xor(abs(Max_index(Safe_Reference_Index)-Max_index(Safe_Sample_Index))>abs(Max_index(Safe_Reference_Index)-Max_index_Second(Safe_Sample_Index)),mean(If_Second_Lower(1:Safe_Sample_Index))>0.5)
        disp('Film facing down.');
    else
        disp('Film facing up.');
    end
        
end
if Only_Calc_The_First == 1
    Active_Array_Size=1;
else
    %Active_Array_Size=length(Considered_Lateral_Index_Range(1:Last_With_Second_Index));
    Active_Array_Size=length(Considered_Lateral_Index_Range);
end

Signal_Sample_1(1:length(Position_micron),1:length(Considered_Lateral_Index_Range))=0;
Signal_Sample_2(1:length(Position_micron),1:length(Considered_Lateral_Index_Range))=0;

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
    Signal_Reference((Max_index_Ref(p)+Window_Size_index):end,p)=0;
    Signal_Reference(1:(Max_index_Ref(p)-Window_Size_index),p)=0;
end


plot(Position_micron(Position_micron<50),Signal_Sample_1(Position_micron<50,30),Position_micron(Position_micron<50),Signal_Sample_2(Position_micron<50,30),Position_micron(Position_micron<50),Signal_Reference(Position_micron<50,30));
xlabel('Optical Path Difference (micron)');
ylabel('Amplitude (a.u.)');
%
%plot(Position_micron(Position_micron<50),Signal_Sample_1(Position_micron<50,1),Position_micron(Position_micron<50),Signal_Sample_1(Position_micron<50,2));

Spectrum_Sample_1=ifft(Signal_Sample_1,[],1);
%Spectrum_Sample_1=2*Spectrum_Sample_1(1:N_f,:);

Spectrum_Sample_2=ifft(Signal_Sample_2,[],1);
%Spectrum_Sample_2=2*Spectrum_Sample_2(1:N_f,:);

Spectrum_Reference=ifft(Signal_Reference,[],1);

%plot(Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Sample_1(Wavelength_micron<1,1)),Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Sample_2(Wavelength_micron<1,1)),Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Sample_1(Wavelength_micron<1,2)));           %�Y�O�W�иѪR�פӮt, �O���O�]�|�]��Signal 
%xlabel('Wavelength (micron)');
%ylabel('Normalized Spectrum');
%legend('Front Interface','Rear Interface','Reference Interface');

%
Spectrum_Sample_1=Spectrum_Sample_1(:,1:Active_Array_Size);
Spectrum_Sample_2=Spectrum_Sample_2(:,1:Active_Array_Size);   %!!!!!!!!!! �]�����U�ӴN�n��n k �F, ��S����Ӥ�������i�ӱ��p�|�ܩ�

Spectrum_Devided=Spectrum_Sample_1./Spectrum_Reference;
Spectrum_Devided_2=Spectrum_Sample_2./Spectrum_Reference;
Spectrum_Devided_3=Spectrum_Sample_2./Spectrum_Sample_1;

%%%%% Till Here 130205



Thickness_Profile=abs((Max_index_Second-Max_index)*(Position_micron(2)-Position_micron(1)));

clear Signal_Sample_First Signal_Sample_Second


% Global arrays generation


Max_Frequency=c/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=c/(Max_Wavelength*1E-9);             %Hz

Frequency_Center=c/(Wavelength_Center*1E-9);
Frequency_Considered_Min=c/(Wavelength_Considered_Max*1E-9);             %Hz
Frequency_Considered_Max=c/(Wavelength_Considered_Min*1E-9);             %Hz


cd(Sample_Path);

Frequency_Center_Index=find(Frequency>Frequency_Center,1,'first');
Frequency_Considered_Min_Index=find(Frequency>Frequency_Considered_Min,1,'first');
Frequency_Considered_Max_Index=find(Frequency>Frequency_Considered_Max,1,'first');

% T 


cd(sprintf('%s\\',Spectroscopy_Path));
Data_Spectroscopy=importdata('1208030331_5micron 1.jws.txt');
SDOCT_n=importdata('SDOCT_n.txt');
SDOCT_k=importdata('SDOCT_k.txt');
SDOCT_T=importdata('SDOCT_T_Best.txt');
SDOCT_Wavelength_micron=importdata('SDOCT_Wavelength_micron_Considered.txt');

Spectrum_Spectroscopy_Old=Data_Spectroscopy(end:-1:1,2)/100;
Wavelength_Spectroscopy=Data_Spectroscopy(end:-1:1,1);     
Frequency_Spectroscopy=c./(Wavelength_Spectroscopy*1E-9);

T=interp1(Frequency_Spectroscopy,Spectrum_Spectroscopy_Old,Frequency,'spline'); 

% Theory - Sample Model (n1 - n - n2)

n2=1;

C1 = 1.03961212; 
C2 = 0.00600069867; 
C3 = 0.231792344; 
C4 = 0.0200179144; 
C5 = 1.01046945; 
C6 = 103.560653;


n_bk7=(C1*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C2)+C3*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C4)+C5*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C6)+1).^0.5;

n_bk7=abs(n_bk7);

n_bk7(isnan(n_bk7))=0;
n1=n_bk7;

% To generate the original spectrum, assuming reference is BK7, too

r_BK7=((1-n_bk7)./(n_bk7+1));
t_AN100=(2*(n_bk7)./(n_bk7+1));

%

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


r_BK7_Considered=repmat(r_BK7(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
r_BK7_Considered=downsample(r_BK7_Considered,Frequency_Downsample_Ratio);
t_AN100_Considered=repmat(t_AN100(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
t_AN100_Considered=downsample(t_AN100_Considered,Frequency_Downsample_Ratio);

Cexp=unwrap(Cexp,[],1);

%plot(Wavelength_micron_Considered,Aexp,Wavelength_micron_Considered,Bexp);
%xlabel('Wavelength (micron)');
%ylabel('Normalized Inteference Spectrum');
for qq=1:Active_Array_Size
    Phase_diff=Cexp(Center_Wavelength_micron_Index,qq);
    NN=fix(Phase_diff/(2*pi));
    Cexp(:,qq)=Cexp(:,qq)-NN*2*pi;
    %while Phase_diff > pi
    %    if Phase_diff>pi
    %        Cexp(:,qq)=Cexp(:,qq)-2*pi;
    %    elseif Phase_diff<pi
    %        Cexp(:,qq)=Cexp(:,qq)+2*pi;
    %    end
    %    Phase_diff=Cexp(Center_Wavelength_micron_Index,qq);
    %end
end

%plot(Wavelength_micron_Considered,Cexp);
%xlabel('Wavelength (micron)');
%ylabel('Continuous Phase Spectrum');


Thickness_Profile=repmat(Thickness_Profile(1:Active_Array_Size),length(Frequency_Considered),1);
%

t1= @ (n) (2*n1_Considered)./(n1_Considered+n); 
t1_r= @ (n) (2*n)./(n2_Considered+n);
t2= @ (n) (2*n)./(n+n2_Considered); 
r1= @ (n) (n1_Considered-n)./(n1_Considered+n);
r1_r= @ (n) (n-n1_Considered)./(n+n1_Considered);
r2= @ (n) (n-n2_Considered)./(n+n2_Considered);
                    %A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now) exp(i*4*pi.*Frequency_Considered.*Position_Upper2Reference_Now/c).*Ratio_Upper2Reference_Now./r_BK7_Considered.*(r1(n)+Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(i*2.*Frequency_Considered.*n.*Thickness_Now/c)./(1-r1_r(n).*r2(n).*exp(i*2*Frequency_Considered.*n.*Thickness_Now/c))));
                    %A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now) abs(exp(i*4*pi.*Frequency_Considered.*Position_Upper2Reference_Now/c).*Ratio_Upper2Reference_Now./r_BK7_Considered.*(r1(n)));
                    %B = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now) abs(exp(i*4*pi.*Frequency_Considered.*Position_Upper2Reference_Now/c).*Ratio_Upper2Reference_Now./r_BK7_Considered.*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(i*2.*Frequency_Considered.*n.*Thickness_Now/c))));
A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) abs(Ratio_Upper2Reference_Now./r_BK7_Considered.*(r1(n)));
B = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) abs(Ratio_Upper2Reference_Now./r_BK7_Considered.*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(1i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c))));
C = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) unwrap(angle(1./r1(n).*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(1i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c)))),[],1); %���M�S�����諸phase, ���O�O���O���ӥushift pi�����?
                   
det_array = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9) Q1.*Q5.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 -Q3.*Q5.*Q7 -Q2.*Q4.*Q9 -Q1.*Q6.*Q8;
                    %A = @ (n,Ratio_Upper2Reference_Now,Position_Upper2Reference) exp(i*4*pi.*Frequency_Considered.*Position_Upper2Reference_Now/c).*Ratio_Upper2Reference_Now./r_BK7_Considered.*r1(n);
                    
                    %A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now) Ratio_Lower2Upper_Now./r1(n).*(t1(n).*t1_r(n).*r2(n).*exp(i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c)./(1-r1_r(n).*r2(n).*exp(i*4*pi*Frequency_Considered.*n.*Thickness_Now/c)));                    
                    %T_abs= @ (n,Thickness_Now) abs((2.*n2_Considered./(n2_Considered+n1_Considered)).*(2*n1_Considered./(n+n1_Considered)).*(2*n./(n+n2_Considered)).*exp(i*2*pi.*Frequency_Considered/c.*n.*Thickness_Now)./(1-((n-n1_Considered)./(n+n1_Considered)).*((n-n2_Considered)./(n+n2_Considered)).*exp(i*4*pi.*Frequency_Considered/c.*n.*Thickness_Now))).^2;
T_abs = @ (n,Thickness_Now) abs(t_AN100_Considered.*t1(n).*t2(n).*exp(1i*2*pi*Frequency_Considered.*n.*Thickness_Now/c)./(1-r1_r(n).*r2(n).*exp(1i*4*pi*Frequency_Considered.*n.*Thickness_Now/c))).^2;          %Note, 120827 remove 2 from ���l4*pi > 2*pi
                    
dn=1E-10;
dt=1E-16;

n(1:length(Frequency_Considered),1:Active_Array_Size)=n_should;
Thickness_Now(1:length(Frequency_Considered),1:Active_Array_Size)=Thickness*Thickness_Profile/max(max(Thickness_Profile));

% Start the n k fitting    
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
                    Number_of_Loop_Checking=25;   
                    qqqq=1;
                    while (Current_Loop<Number_of_Loop)
                        if mod(Current_Loop,Number_of_Loop_Checking)==0
                            CONT_left=1;
                            CONT_right=1;
                            while(CONT_left || CONT_right)
                                index_needshift_left=find(abs(diff(n(1:Center_Wavelength_micron_Index)))>0.008,1,'last');
                                index_needshift_right=find(abs(diff(n((Center_Wavelength_micron_Index+1):end)))>0.008,1,'first')+Center_Wavelength_micron_Index;
                                CONT_left=0;
                                CONT_right=0;
                                if (index_needshift_left>5)
                                    CONT_left=1;
                                    n(1:index_needshift_left)=n(1:index_needshift_left)-n(index_needshift_left)+n(index_needshift_left+1);
                                elseif (index_needshift_right<(length(n)-5)) 
                                    CONT_right=1;
                                    n((index_needshift_right+1):end)=n((index_needshift_right+1):end)-n(index_needshift_right+1)+n(index_needshift_right);
                                end

                            end

                        end
                        dA_dn_real = (A(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dB_dn_real = (B(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dC_dn_real = (C(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                    
                        dA_dn_imag = (A(n+1i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dB_dn_imag = (B(n+1i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dC_dn_imag = (C(n+1i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        
                        
                        dA_dthickness = (A(n,Thickness_Now+dt,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dt;
                        dB_dthickness = (B(n,Thickness_Now+dt,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dt;
                        dC_dthickness = (C(n,Thickness_Now+dt,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dt;
                        
                        C_Temp=C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now);
                        for qq=1:Active_Array_Size
                            Phase_diff=C_Temp(Center_Wavelength_micron_Index,qq);
                            NN=fix(Phase_diff/(2*pi));
                            C_Temp(:,qq)=C_Temp(:,qq)-NN*2*pi;
                        end
                        
                        delta_A=Aexp-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now);
                        delta_B=Bexp-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now);
                        delta_C=Cexp-C_Temp;
                        
                        %dn_real=(delta_A.*dB_dn_imag-dA_dn_imag.*delta_B)./(dA_dn_real.*dB_dn_imag-dA_dn_imag.*dB_dn_real);
                        %dn_imag=(dA_dn_real.*delta_B-delta_A.*dB_dn_real)./(dA_dn_real.*dB_dn_imag-dA_dn_imag.*dB_dn_real);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        dn_real = det_array(delta_A,dA_dn_imag,dA_dthickness,delta_B,dB_dn_imag,dB_dthickness,delta_C,dC_dn_imag,dC_dthickness)./det_array(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness);
                        dn_imag = det_array(dA_dn_real,delta_A,dA_dthickness,dB_dn_real,delta_B,dB_dthickness,dC_dn_real,delta_C,dC_dthickness)./det_array(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness);
                        dthickness = det_array(dA_dn_real,dA_dn_imag,delta_A,dB_dn_real,dB_dn_imag,delta_B,dC_dn_real,dC_dn_imag,delta_C)./det_array(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness);
                        
                        n=n+0.5*(dn_real+1i*dn_imag);
                        Thickness_Now=Thickness_Now+0.5*dthickness;   %How if I make the thickness a array during the iteration?
                        Current_Loop=Current_Loop+1;
                        %fprintf('%d/%d\n',Current_Loop,Number_of_Loop);
                        n(isnan(n))=n_should;
                        n(isinf(n))=n_should;
                        if isnan(Thickness_Now)>0
                            Thickness_Now(isnan(Thickness_Now))=Thickness*Thickness_Profile(isnan(Thickness_Now))/max(max(Thickness_Profile));
                            Thickness_Now(isinf(Thickness_Now))=Thickness*Thickness_Profile(isnan(Thickness_Now))/max(max(Thickness_Profile));
                        end
                        %if Current_Loop>70
                        %    n_Record(:,qqqq)=n;
                        %    Thickness_Record(:,qqqq)=Thickness_Now;
                        %    qqqq=qqqq+1;
                        %end
                        
                    end
                        Current_Progress=100*((p-1)+((q-1)+(w-1)/length(Ratio_Upper2Reference))/length(Ratio_Lower2Upper))/length(Thickness);
                        fprintf('%f%%\n',Current_Progress);
                        %Merit=sum(Weight_Function.*(T_Considered-T_abs(n,Thickness_Now)).^2);
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

%% Time Frequency Analysis

%plot(Wavelength_micron_Considered,real(n_Best(:,1)),Wavelength_micron_Considered,n2_Considered(:,1));
xlabel('Wavelength (micron)');
ylabel('n');

%legend('Calculated Film Refractive Index','Refractive Index of Substrate')

%plot(Wavelength_micron_Considered,imag(n_Best(:,1)));
xlabel('Wavelength (micron)');
ylabel('k');

%%

%plot(Wavelength_micron_Considered,T_Best(:,1),Wavelength_micron_Considered,T_Considered(:,1));
xlabel('Wavelength (micron)');
ylabel('T');
Data_Reduction_Factor=25;
clear FFOCT_n FFOCT_k FFOCT_T FFOCT_Wavelength_micron_Considered
FFOCT_Start_Wavelength=0.51;
FFOCT_End_Wavelength=0.56;
FFOCT_Wavelength_Spacing=0.005;
FFOCT_Wavelength_micron_Considered=FFOCT_Start_Wavelength:FFOCT_Wavelength_Spacing:FFOCT_End_Wavelength;

FFOCT_n=interp1(Wavelength_micron_Considered,real(n_Best),FFOCT_Wavelength_micron_Considered);
FFOCT_k=interp1(Wavelength_micron_Considered,imag(n_Best),FFOCT_Wavelength_micron_Considered);
FFOCT_T=interp1(Wavelength_micron_Considered,T_Best,FFOCT_Wavelength_micron_Considered);


subplot(3,1,1)
plot(FFOCT_Wavelength_micron_Considered(abs(FFOCT_Wavelength_micron_Considered-0.54)<0.05),FFOCT_n(abs(FFOCT_Wavelength_micron_Considered-0.54)<0.05),SDOCT_Wavelength_micron(abs(SDOCT_Wavelength_micron-0.54)<0.05),SDOCT_n(abs(SDOCT_Wavelength_micron-0.54)<0.05));

set(gca,'FontSize',10);
xlabel('Wavelength (micron)');
ylabel('n');

subplot(3,1,2)
plot(FFOCT_Wavelength_micron_Considered(abs(FFOCT_Wavelength_micron_Considered-0.54)<0.05),FFOCT_k(abs(FFOCT_Wavelength_micron_Considered-0.54)<0.05),SDOCT_Wavelength_micron(abs(SDOCT_Wavelength_micron-0.54)<0.05),SDOCT_k(abs(SDOCT_Wavelength_micron-0.54)<0.05));

set(gca,'FontSize',10);
xlabel('Wavelength (micron)');
ylabel('k');

subplot(3,1,3)
plot(FFOCT_Wavelength_micron_Considered(abs(FFOCT_Wavelength_micron_Considered-0.54)<0.05),FFOCT_T(abs(FFOCT_Wavelength_micron_Considered-0.54)<0.05),SDOCT_Wavelength_micron(abs(SDOCT_Wavelength_micron-0.54)<0.05),SDOCT_T(abs(SDOCT_Wavelength_micron-0.54)<0.05));

set(gca,'FontSize',10);
xlabel('Wavelength (micron)');
ylabel('T');
legend('FF-OCT','SD-OCT');



dlmwrite('FFOCT_n.txt',FFOCT_n','delimiter','\t','newline','pc');

dlmwrite('FFOCT_k.txt',FFOCT_k','delimiter','\t','newline','pc');

dlmwrite('FFOCT_T.txt',FFOCT_T','delimiter','\t','newline','pc');
dlmwrite('FFOCT_Wavelength_micron_Considered.txt',FFOCT_Wavelength_micron_Considered','delimiter','\t','newline','pc');
