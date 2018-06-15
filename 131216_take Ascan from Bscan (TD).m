clear all

cd('D:\Users\TuanShu\131216\SNR\');

Calibration_Array=importdata('Calibration_Array.txt');
fx=importdata('fx_new.txt');

%Array_Number=119:273;
Array_Number=626:826;

M=620;
N=1000;
LP=10;
%%
Data1_1D_array(1:2040,1:length(Array_Number))=0;
for p=1:length(Array_Number)
    Data1_2D_Temp=dlmread(sprintf('sig_%i.txt',Array_Number(p)));
    Data1_1D_array(:,p)=Data1_2D_Temp(:,M);
    %Data1_2D_New=interp1(Calibration_Array,Data1_2D_Temp,fx);
    %Data1_2D_New(isnan(Data1_2D_New))=0;
    %FFT = abs(fft(Data1_2D_New,[],1));
    %FFT=FFT(1:round(size(FFT,1)/2),:);
    %FFT(1:LP,:)=0;
    %Image=Image+FFT;
    disp(p);
end
    Data1_1D_array_New=interp1(Calibration_Array,Data1_1D_array,fx);
    Data1_1D_array(isnan(Data1_1D_array))=0;
    FFT = abs(fft(Data1_1D_array,[],1));
    FFT=FFT(1:round(size(FFT,1)/2),:);
    FFT(1:LP,:)=0;
    
    SIG=mean(max(FFT));
    STD=std(max(FFT));
    
plot(1:length(Data1_1D_array(N,:)),Data1_1D_array(N,:));
xlim([0 length(Data1_1D_array(N,:))]);
xlabel('Number of Scan');
ylabel('Signal (at the center wavelength)');

STD=std(Data1_1D_array(N,:));

SIG=mean(Data1_1D_array(N,:));

SNR=10*log10((SIG/STD)^2);
