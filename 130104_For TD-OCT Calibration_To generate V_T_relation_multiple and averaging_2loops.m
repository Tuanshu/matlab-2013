clear all

%%% To find the phase from min voltage to max voltage ONLY
%%%��M�o�{�ڸӬ������n�����ӬO�����q�� �Ӥ��O�^�Ǫ��q��

Path='D:\Users\TuanShu\130103\';

Parameter='Frequency';

Frequency=[0.004];
%Frequency=[0.005];
Sampling_rate=500; %Hz
smooth_window=10000;
Range=10;        %V
min_voltage=(10-Range)/2;
max_voltage=min_voltage+Range;
dvoltage_min=0.001;
dvoltage_max=0.1;
Wavelength_Laser=0.6328;     %micron

Lower_Band=20;
Upper_Band=600;%fix(length(Spectrum)/2);

dt=1/Sampling_rate;
cd(Path);
clf
%gcf: handle of current fiqure
%gca: handle of current axes
%%%rect = [left, bottom, width, height]
set(gcf,'Position',get(0,'ScreenSize')+[200 100 -400 -200])             %0�i��O�ù�, gca�O�ثe�ϧΪ�handle (�p�G�S���Ϫ��ܷ|�ۤv�}�s���@��handle)
set(gca,'FontSize',16);

Time_length_reduction_factor=1;
N_Waveform=100000;  %Full period

required_number_of_loop=100;

Data=dlmread(sprintf('130103_f%g_V%d.txt',Frequency,Range));                   %Asuume: line 1: time, line 2: interference signal, line 3: voltage 121219_f0.001_V2to8.txt�ǩ�

Voltage_all=Data(:,1);
Signal_all=Data(:,2);

Estimated_Number_of_Pixel_for_One_Trip=1/Frequency/2*Sampling_rate;
Voltage_UniformBase=min_voltage:(max_voltage-min_voltage)/(Estimated_Number_of_Pixel_for_One_Trip-1):max_voltage;

p=1;
Voltage_previous_max_index=1;
Voltage_min_index=1;
Voltage_next_max_index=1;

while Voltage_previous_max_index<size(Signal_all,1)-2*(Voltage_next_max_index-Voltage_min_index)

    Voltage_min_index_temp=Voltage_previous_max_index+find(Voltage_all((Voltage_previous_max_index+1):end)<(min_voltage+dvoltage_min),1,'first');
    Voltage_next_max_index=Voltage_min_index_temp+find(Voltage_all((Voltage_min_index_temp+1):end)>(max_voltage-dvoltage_max),1,'first');
    [value Voltage_min_index]=min(Voltage_all((Voltage_min_index_temp+1):Voltage_next_max_index));
    Voltage_min_index=Voltage_min_index_temp+Voltage_min_index;
    Voltage_previous_max_index=Voltage_next_max_index;
    
    Voltage=Voltage_all(Voltage_min_index:Voltage_next_max_index);
    Signal=Signal_all(Voltage_min_index:Voltage_next_max_index);

    Time_index=1:length(Signal);

    Spectrum=ifft(Signal);
    Spectrum(1:Lower_Band)=0;
    Spectrum(Upper_Band:end)=0;

    Signal_New=fft(Spectrum,[],1);
    Phase_original=angle(Signal_New);
    Phase=unwrap(Phase_original);

    Position=-1*Wavelength_Laser*Phase/(2*pi)/2;
    %Position=Position-Position(find(Voltage>5,1,'first')); �]���ϥ��n������,�b�����ᰵ�o��Ƥ���n
    %Velocity=diff(Position)
    Position=smooth(Position,smooth_window);
    Voltage=smooth(Voltage,smooth_window);
   
    plot(Voltage,Position);
    
    Position_UniformBase_temp=interp1(Voltage,Position,Voltage_UniformBase);
    Position_UniformBase(:,p)=Position_UniformBase_temp-Position_UniformBase_temp(find(Voltage_UniformBase>=5,1,'first'));
    disp(p);
    p=p+1;
    
end
plot(Voltage_UniformBase,Position_UniformBase);
Acquired_Number_of_Loop=p;


    Position_max=max(Position);
    Position_min=min(Position);
    Position_range=Position_max-Position_min;
    
    Predict_Speed=Position_range/(1/Frequency/2);
    
    Time_old=Position/Predict_Speed;
    
    Time_waveform_1=Time_old(1)*Time_length_reduction_factor:(Time_old(end)*Time_length_reduction_factor-Time_old(1)*Time_length_reduction_factor)/(N_Waveform/2-1):Time_old(end)*Time_length_reduction_factor;
    
    Voltage_waveform_1=interp1(Time_old,Voltage,Time_waveform_1);
    
    Time_waveform_2=Time_waveform_1-Time_waveform_1(1); %start from 0
    Time_waveform(1:length(Time_waveform_2),p)=Time_waveform_2';
    Time_waveform((length(Time_waveform_2)+1):(2*length(Time_waveform_2)),p)=Time_waveform_2(end)+Time_waveform_2';
    Voltage_waveform(:,p)=[Voltage_waveform_1 max(Voltage_waveform_1)-min(Voltage_waveform_1)-Voltage_waveform_1]';
    %Position_relative=Position-(Voltage-5)*10;
    %plot(Voltage,Position_relative);
    %plot(Time_waveform,Voltage_waveform);
 
    
%% To generate time array

%xlabel('Voltage (V)','FontSize',18);
%ylabel('PZT Position (micron)','FontSize',18);
%hleg=legend('Frequency = 0.001 Hz','Frequency = 0.002 Hz','Frequency = 0.005 Hz','Frequency = 0.01 Hz','Frequency = 0.02 Hz','Frequency = 0.05 Hz','Positon','Best');

%plot(Time_waveform,Voltage_waveform);
xlabel('Time waveform (sec)','FontSize',18);
ylabel('Voltage waveform (V)','FontSize',18);
dlmwrite('Time_waveform.txt',Time_waveform,'delimiter','\t','newline','pc');
dlmwrite('Voltage_waveform.txt',Voltage_waveform,'delimiter','\t','newline','pc');