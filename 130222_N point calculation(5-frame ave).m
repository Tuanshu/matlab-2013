clear all

cd('D:\Users\TuanShu\');
Data=dlmread('130222_Tsai N point data.txt'); 

Time=Data(:,1);

Signal=Data(:,2);

Number_of_Point_per_Carrier_Before_Averaging=140;

N_Array=[3 4 5 7 10 20 35 140];
 
Averaging_2=4;
Window_array=ones(Averaging_2,1);

for qq=1:length(N_Array)

N=N_Array(qq);
if N==3
    Averaging=47;
elseif N==4
    Averaging=35;
elseif N==7
    Averaging=20;    
elseif N==10
    Averaging=14;  
elseif N==20
    Averaging=7;  
elseif N==35
    Averaging=4;  
elseif N==140
    Averaging=1;
end


%N=Number_of_Point_per_Carrier_Before_Averaging/Averaging;

%Window_array=ones(Averaging,1);

%Signal_ave=conv(Signal,Window_array,'same')/Averaging;

%C= @ (n,m) factorial(n)/factorial(m)/factorial(n-m);

clear Signal_New Time_New

Signal_New(1:fix(length(Signal)/Averaging))=0;
Time_New(1:fix(length(Signal)/Averaging))=0;
Signal_New=Signal_New';
Time_New=Time_New';


for p=1:fix(length(Signal)/Averaging)
    Signal_New(p)=mean(Signal((1+(p-1)*Averaging):(p*Averaging)));
    Time_New(p)=Time(1+(p-1)*Averaging);
end

Better_Array_Length=fix(length(Signal_New)/N)*N;
Signal_New=Signal_New(1:Better_Array_Length);
Time_New=Time_New(1:Better_Array_Length);

clear Signal_New_Env Time_New_Env
Signal_New_Env(1:fix(length(Signal_New)/N))=0;
Time_New_Env(1:fix(length(Time_New)/N))=0;

for p=1:fix(length(Signal_New)/N)
    Temp=0;
    for m=1:N
       for n=1:N
           Temp=Temp+(Signal_New((p-1)*N+m)-Signal_New((p-1)*N+n))^2;
       end
    end
    Signal_New_Env(p)=(Temp^0.5)/N;
    Time_New_Env(p)=Time_New(1+(p-1)*N);
end

Signal_New_Env=conv(Signal_New_Env,Window_array,'same')/Averaging_2;

if qq==1
    Signal_New_Env_1=Signal_New_Env;
    Time_New_Env_1=Time_New_Env;
elseif qq==2
    Signal_New_Env_2=Signal_New_Env;
    Time_New_Env_2=Time_New_Env;
elseif qq==3
    Signal_New_Env_3=Signal_New_Env;
    Time_New_Env_3=Time_New_Env;
elseif qq==4
    Signal_New_Env_4=Signal_New_Env;
    Time_New_Env_4=Time_New_Env;
elseif qq==5
    Signal_New_Env_5=Signal_New_Env;
    Time_New_Env_5=Time_New_Env;   
elseif qq==6
    Signal_New_Env_6=Signal_New_Env;
    Time_New_Env_6=Time_New_Env;
elseif qq==7
    Signal_New_Env_7=Signal_New_Env;
    Time_New_Env_7=Time_New_Env;
end
end
%plot(Time_New,Signal_New,Time_New_Env_1,(Signal_New_Env_1));
plot(Time_New_Env_1,(Signal_New_Env_1),Time_New_Env_2,(Signal_New_Env_2),Time_New_Env_3,(Signal_New_Env_3),Time_New_Env_4,(Signal_New_Env_4),Time_New_Env_5,(Signal_New_Env_5),Time_New_Env_6,(Signal_New_Env_6),Time_New_Env_7,(Signal_New_Env_7));

plot(Time_New_Env_1,log10(Signal_New_Env_1)*10,Time_New_Env_2,log10(Signal_New_Env_2)*10,Time_New_Env_3,log10(Signal_New_Env_3)*10,Time_New_Env_4,log10(Signal_New_Env_4)*10,Time_New_Env_5,log10(Signal_New_Env_5)*10,Time_New_Env_6,log10(Signal_New_Env_6)*10,Time_New_Env_7,log10(Signal_New_Env_7)*10);

%plot(Time_New_Env_1,log10(Signal_New_Env_1-Signal_New_Env_1(end))*10,Time_New_Env_2,log10(Signal_New_Env_2-Signal_New_Env_2(end))*10,Time_New_Env_3,log10(Signal_New_Env_3-Signal_New_Env_3(end))*10);
plot(Time_New_Env_1,log10(Signal_New_Env_1)*10,Time_New_Env_2,log10(Signal_New_Env_2)*10,Time_New_Env_3,log10(Signal_New_Env_3)*10,Time_New_Env_4,log10(Signal_New_Env_4)*10,Time_New_Env_5,log10(Signal_New_Env_5)*10,Time_New_Env_6,log10(Signal_New_Env_6)*10,Time_New_Env_7,log10(Signal_New_Env_7)*10);
xlabel('Position (micron)');
ylabel('Amplitude (log)');
legend('3-point','4-point','7-point','10-point','20-point','35-point','140-point');


