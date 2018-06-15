%% Options

clear all


Sample_Path='D:\Users\TuanShu\130512\Y2\sam1\ref\D0.txt';
cd('D:\Users\TuanShu\');

Frequency_Downsample_Ratio=5;

    n=dlmread('n.txt');
    k=dlmread('k.txt');
    G=dlmread('G.txt');
G=abs(G)./max(abs(G));
Ratio_Lower2Upper=1;     %if the rear interface has larger interference eff than front interface, Ratio_Lower2Upper>1

%%%
%%%
Wavelength_Center=540;

Center_Wavelength_micron=0.54;
Wavelength_Considered_Min=510; %510         %nm
Wavelength_Considered_Max=590;  %580

Max_Wavelength=800;             %nm
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


Data=importdata(Sample_Path);      % Data_2: the glass data 1
Wavelength=Data(:,1);           %nm
Frequency_Old=c./(Wavelength*1E-9);
Frequency=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency=Frequency';
Wavelength_micron=(c./Frequency)*1E6;

Frequency_Center_Index=find(Frequency>Frequency_Center,1,'first');
Frequency_Considered_Min_Index=find(Frequency>Frequency_Considered_Min,1,'first');
Frequency_Considered_Max_Index=find(Frequency>Frequency_Considered_Max,1,'first');

Wavelength_micron_Considered=Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
Wavelength_micron_Considered=downsample(Wavelength_micron_Considered,Frequency_Downsample_Ratio);
Wavelength_nm_Considered=Wavelength_micron_Considered*1000;
% Time-domain

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=c*Time;
Position_micron=Position*1E6;

%%% T 




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

    t_AN100=(2*(1)./(n_AN100+1));
    n1=n_AN100;
    n2=1;


%%% To generate the original spectrum

r_AN100=((1-n_AN100)./(n_AN100+1));

%%% Data Loading

    n1_Considered=n1(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
    n1_Considered=downsample(n1_Considered,Frequency_Downsample_Ratio);
    n2_Considered=repmat(n2,size(n1_Considered,1),1);

t1= @ (n) (2*n1_Considered)./(n1_Considered+n); 
t1_r= @ (n) (2*n)./(n1_Considered+n);
t2= @ (n) (2*n)./(n+n2_Considered); 
r1= @ (n) (n1_Considered-n)./(n1_Considered+n);
r1_r= @ (n) (n-n1_Considered)./(n+n1_Considered);
r2= @ (n) (n-n2_Considered)./(n+n2_Considered);


A= @ (n,k) abs(G.*r1(n+i*k));
E= @ (n,k) abs(G.*(t1_r(n+i*k).*t1(n+i*k)./r1_r(n+i*k)));

dn=0.000000001;
plot(Wavelength_nm_Considered,abs(G));
xlabel('Wavelength (nm)');
ylabel('G');

dA_dn = (A(n+dn,k)-A(n,k))/dn;
dE_dn = (E(n+dn,k)-E(n,k))/dn;
dA_dk = (A(n,k+dn)-A(n,k))/dn;
dE_dk = (E(n,k+dn)-E(n,k))/dn;
%%

subplot(2,1,1);
plot(Wavelength_nm_Considered,A(n,k));
xlabel('Wavelength (nm)');
ylabel('A');

subplot(2,1,2);
plot(Wavelength_nm_Considered,E(n,k));
xlabel('Wavelength (nm)');
ylabel('E');
%%

subplot(2,2,1);
plot(Wavelength_nm_Considered,dA_dn);
xlabel('Wavelength (nm)');
ylabel('dA/dn');
set(gca,'XTick',500:20:600);

subplot(2,2,2);
plot(Wavelength_nm_Considered,dA_dk);
xlabel('Wavelength (nm)');
ylabel('dA/dk');
set(gca,'XTick',500:20:600);

subplot(2,2,3);
plot(Wavelength_nm_Considered,dE_dn);
xlabel('Wavelength (nm)');

ylabel('dE/dn');
set(gca,'XTick',500:20:600);

subplot(2,2,4);
plot(Wavelength_nm_Considered,dE_dk);
xlabel('Wavelength (nm)');

ylabel('dE/dk');
set(gca,'XTick',500:20:600);
