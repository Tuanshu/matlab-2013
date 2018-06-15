clear all

Mode=2;

%Mode=1 the number of pixels is equal to the original data
%Mode=2 the number of pixels is number of period

Ndataset=2;

if Ndataset==1
    cd('D:\Users\TuanShu\130228_FF-OCT data\');
    Data=dlmread('filtered.txt'); 
    Data=Data(1:35000,:);
    
elseif Ndataset==2
    cd('D:\Users\TuanShu\130429_for CCTsai_Compare N=3 and more\');
    Data(:,1)=dlmread('Position.txt'); 
    Data(:,2)=dlmread('130429_Original_Data.txt'); 

elseif Ndataset==3   
    cd('D:\Users\TuanShu\130429_for CCTsai_Compare N=3 and more\');
    Data(:,1)=dlmread('Position.txt'); 
    Data(:,2)=dlmread('130429_Differentiated_Data.txt'); 
    
elseif Ndataset==4          %generated data  
    cd('D:\Users\TuanShu\130429_for CCTsai_Compare N=3 and more\');
    Data(:,1)=dlmread('Position.txt'); 
    Temp_Signal=dlmread('130429_Differentiated_Data.txt'); 
    %[maxvalue maxindex]=max(Temp_Signal);
    %[minvalue minindex]=min(Temp_Signal(maxindex:end));
    %minindex=minindex+maxindex-1;
    %[secondmax secondmaxindex]=max(Temp_Signal(minindex:end));
    period=150;%abs(secondmax-maxindex);
    index_array=1:length(Data(:,1));
    Generated_Data=cos(index_array/period*2*pi);
    
elseif Ndataset==5          %generated data  
    cd('D:\Users\TuanShu\130429_for CCTsai_Compare N=3 and more\');
    Data(:,1)=dlmread('Position.txt'); 
    Temp_Signal=dlmread('130429_Differentiated_Data.txt'); 
    %[maxvalue maxindex]=max(Temp_Signal);
    %[minvalue minindex]=min(Temp_Signal(maxindex:end));
    %minindex=minindex+maxindex-1;
    %[secondmax secondmaxindex]=max(Temp_Signal(minindex:end));
    index_array=1:length(Data(:,1));
    period=152;%abs(secondmax-maxindex);  %(index)
    center_near=25000;                         %(index)
    center=find(cos(index_array(center_near:end)/period*2*pi)>0.999,1,'first')+center_near-1;
    fwhm=830;
    index_array=1:length(Data(:,1));
    Generated_Data=gaussmf(index_array,[fwhm/((2*log(2))^0.5)/2 center]).*cos(index_array/period*2*pi);
    Data(:,2)=Generated_Data;
    
end

Position_UpperLimit=110;

Position=Data(:,1);

Signal=Data(:,2);

Signal=Signal(Position<Position_UpperLimit);
Position=Position(Position<Position_UpperLimit);

Number_of_Point_per_Carrier_Before_Averaging=140;

N_Array=[3 4 38 152];

Number_of_Points_per_Period=152;

 
Averaging_2=1;
Window_array_2=ones(Averaging_2,1);
NN=1;

for qq=1:length(N_Array)

    N=N_Array(qq);
    
    %Number_of_Points_per_Period=N_Array(qq);
    %if N==3
    %    Averaging=48;
    %elseif N==4
    %    Averaging=36;
    %elseif N==6
    %    Averaging=24;    
    %elseif N==9
    %    Averaging=16;  
    %elseif N==18
    %    Averaging=8;  
    %elseif N==36
    %    Averaging=4;  
    %elseif N==144
    %    Averaging=1;
    %end
    Averaging=fix(Number_of_Points_per_Period/N);

    Window_array_1=ones(Averaging,1);

    %N=Number_of_Point_per_Carrier_Before_Averaging/Averaging;

    %Window_array=ones(Averaging,1);

    %Signal_ave=conv(Signal,Window_array,'same')/Averaging;

    %C= @ (n,m) factorial(n)/factorial(m)/factorial(n-m);

    clear Signal_New Position_New

    %Signal_New(1:fix(length(Signal)/Averaging))=0;
    %Position_New(1:fix(length(Signal)/Averaging))=0;
    %Signal_New=Signal_New';
    %Position_New=Position_New';

    %for p=1:fix(length(Signal)/Averaging)

    %    Signal_New(p)=mean(Signal((1+(p-1)*Averaging):(p*Averaging)));

    %    Position_New(p)=Position(1+(p-1)*Averaging);
    %    disp(p);
    %end
    
    Signal_New=conv(Signal,Window_array_1,'same')/Averaging;
    Position_New=Position;
    
    %Signal_New(1:fix(length(Signal)))=0;
    %Position_New(1:fix(length(Signal)))=0;
    %Signal_New=Signal_New';
    %Position_New=Position_New';

    %for p=1:fix(length(Signal))
    %    if p+Averaging < length(Signal)
    %        Signal_New(p)=mean(Signal(p:(p+Averaging)));
    %    else
    %        Signal_New(p)=mean(Signal(p:end));
    %    end        
    %    Position_New(p)=Position(p);
    %    disp(p);
    %end

    %for p=1:fix(length(Signal))
    %    if p+Averaging < length(Signal)
    %        Signal_New(p)=mean(Signal(p:(p+Averaging)));
    %    else
    %        Signal_New(p)=mean(Signal(p:end));
    %    end
    %    Position_New(p)=Position(p);
    %    disp(p);
    %end
    
    
    %Better_Array_Length=fix(length(Signal_New)/N)*N;
    %Signal_New=Signal_New(1:Better_Array_Length);
    %Position_New=Position_New(1:Better_Array_Length);
    %Position_New=Position;
    
    clear Signal_New_Env Position_New_Env

    %for p=1:fix(length(Signal_New)/N)
    %    Temp=0;
    %    for m=1:N
    %       for n=1:N
    %           Temp=Temp+(Signal_New((p-1)*N+m)-Signal_New((p-1)*N+n))^2;
    %       end
    %    end
    %    Signal_New_Env(p)=(Temp^0.5)/N;
    %    Position_New_Env(p)=Position_New(1+(p-1)*N);
    %end
    %for p=1:(fix((length(Signal_New)-N)/NN))
    %    Temp=0;
    %    for m=1:N
    %       for n=1:N
    %           Temp=Temp+(Signal_New(1+(p-1)*NN+m)-Signal_New(1+(p-1)*NN+n))^2;
    %       end
    %    end
    %    Signal_New_Env(p)=(Temp^0.5)/N;
    %    Position_New_Env(p)=Position_New(1+(p-1)*NN);
    %    disp(p);
    %end
    if Mode == 1
        
        Signal_New_Env(1:fix(length(Signal_New)))=0;
        Position_New_Env(1:fix(length(Position_New)))=0;
        for p=1:(fix((length(Signal_New))))
            Temp=0;
            for m=1:N
               for n=1:N
                   Temp=Temp+(Signal_New(min(p+(m-1)*Averaging,length(Signal_New)))-Signal_New(min(p+(n-1)*Averaging,length(Signal_New))))^2;
               end
            end
            Signal_New_Env(p)=(Temp^0.5)/N;
            Position_New_Env(p)=Position_New(p)+Position(Number_of_Points_per_Period)/2;
            disp(p);
        end
    elseif Mode ==2
        
        Signal_New_Env(1:fix(length(Signal_New)/Number_of_Points_per_Period))=0;
        Position_New_Env(1:fix(length(Position_New)/Number_of_Points_per_Period))=0;
        for p=1:(fix((length(Signal_New))/Number_of_Points_per_Period))
            Temp=0;
            for m=1:N
               for n=1:N
                   Temp=Temp+(Signal_New(min(1+(p-1)*Number_of_Points_per_Period+(m-1)*Averaging,length(Signal_New)))-Signal_New(min(1+(p-1)*Number_of_Points_per_Period+(n-1)*Averaging,length(Signal_New))))^2;
               end
            end
            Signal_New_Env(p)=(Temp^0.5)/N;
            Position_New_Env(p)=Position_New(1+(p-1)*Number_of_Points_per_Period)+Position(Number_of_Points_per_Period)/2;
            disp(p);
        end
    end
    Signal_New_Env=conv(Signal_New_Env,Window_array_2,'same')/Averaging_2;

    if qq==1
        Signal_New_Env_1=Signal_New_Env;
        Position_New_Env_1=Position_New_Env;
    elseif qq==2
        Signal_New_Env_2=Signal_New_Env;
        Position_New_Env_2=Position_New_Env;
    elseif qq==3
        Signal_New_Env_3=Signal_New_Env;
        Position_New_Env_3=Position_New_Env;
    elseif qq==4
        Signal_New_Env_4=Signal_New_Env;
        Position_New_Env_4=Position_New_Env;
    elseif qq==5
        Signal_New_Env_5=Signal_New_Env;
        Position_New_Env_5=Position_New_Env;   
    elseif qq==6
        Signal_New_Env_6=Signal_New_Env;
        Position_New_Env_6=Position_New_Env;
    elseif qq==7
        Signal_New_Env_7=Signal_New_Env;
        Position_New_Env_7=Position_New_Env;
    end
end

FFT=fft(Data(:,2));

FFT(1:250)=0;

FFT(500:end)=0;
Data_filtered=2*real(ifft(FFT));
Data_hilbert=2*abs(ifft(FFT));

plot(Position,Data(:,2),Position,Data_filtered,Position,Data_hilbert,Position_New_Env_1,Signal_New_Env_1,Position_New_Env_2,Signal_New_Env_2,Position_New_Env_3,Signal_New_Env_3,Position_New_Env_4,Signal_New_Env_4);
xlabel('Position (micron)');
ylabel('Amplitude (a.u.)');



Output=[Position_New_Env_1' Signal_New_Env_1' Signal_New_Env_2' Signal_New_Env_3' Signal_New_Env_4'];

Original=[Position Data(:,2) Data_filtered Data_hilbert];

dlmwrite('Output.txt',Output,'delimiter','\t','newline','pc');

dlmwrite('Original.txt',Original,'delimiter','\t','newline','pc');
%plot(Position_New_Env_1,Signal_New_Env_1,Position,Data(:,2));
%xlabel('Position (micron)');
%ylabel('Amplitude (a.u.)');
%plot(Position_New,Signal_New,Position_New_Env_1,(Signal_New_Env_1));
%plot(Position_New_Env_1,(Signal_New_Env_1),Position_New_Env_2,(Signal_New_Env_2),Position_New_Env_3,(Signal_New_Env_3),Position_New_Env_4,(Signal_New_Env_4),Position_New_Env_5,(Signal_New_Env_5),Position_New_Env_6,(Signal_New_Env_6),Position_New_Env_7,(Signal_New_Env_7));

%plot(Position_New_Env_1,log10(Signal_New_Env_1)*10,Position_New_Env_2,log10(Signal_New_Env_2)*10,Position_New_Env_3,log10(Signal_New_Env_3)*10,Position_New_Env_4,log10(Signal_New_Env_4)*10,Position_New_Env_5,log10(Signal_New_Env_5)*10,Position_New_Env_6,log10(Signal_New_Env_6)*10,Position_New_Env_7,log10(Signal_New_Env_7)*10);

%plot(Position_New_Env_1,log10(Signal_New_Env_1-Signal_New_Env_1(end))*10,Position_New_Env_2,log10(Signal_New_Env_2-Signal_New_Env_2(end))*10,Position_New_Env_3,log10(Signal_New_Env_3-Signal_New_Env_3(end))*10);
%plot(Position_New_Env_1,log10(Signal_New_Env_1)*10,Position_New_Env_2,log10(Signal_New_Env_2)*10,Position_New_Env_3,log10(Signal_New_Env_3)*10,Position_New_Env_4,log10(Signal_New_Env_4)*10,Position_New_Env_5,log10(Signal_New_Env_5)*10,Position_New_Env_6,log10(Signal_New_Env_6)*10,Position_New_Env_7,log10(Signal_New_Env_7)*10);
%xlabel('Position (micron)');
%ylabel('Amplitude (dB)');
%legend('3-point','4-point','6-point','9-point','18-point','36-point','144-point');

%plot(Position_New_Env_1,Signal_New_Env_1,Position_New_Env_2,Signal_New_Env_2,Position_New_Env_3,Signal_New_Env_3,Position_New_Env_4,Signal_New_Env_4,Position_New_Env_5,Signal_New_Env_5,Position_New_Env_6,Signal_New_Env_6,Position_New_Env_7,Signal_New_Env_7);
%xlabel('Position (micron)');
%ylabel('Amplitude (a.u.)');
%legend('3-point','4-point','6-point','9-point','18-point','36-point','144-point');

%legend('148-point','150-point','152-point');

%dlmwrite('Signal_filtered.txt',Signal_New_Env_1','delimiter','\t','newline','pc');

%dlmwrite('Position.txt',Position_New_Env_1','delimiter','\t','newline','pc');