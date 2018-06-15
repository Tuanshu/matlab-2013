clear all

cd('D:\Users\TuanShu\131224\');

Calibration_Array=importdata('Calibration_Array.txt');
fx=importdata('fx_new.txt');

NofPixil=2040;
XX=1:2040;
%%
Weight_function=gaussmf(XX,[2040/6 NofPixil/2]);         %to suppress the signal near the boundary of CMOS
plot(Weight_function);
ylim([0 1]);
%%
%Array_Number=119:273;
Array_Number=875:925;

LP=100;
%%
Image(1:2048,1:1088)=0;
%Weight_function_2D=repmat()
for p=1:length(Array_Number)
    Data1_2D_Temp=importdata(sprintf('131224_FOREARM_1_%i',Array_Number(p)));
    Ground=(Data1_2D_Temp(1,:)+Data1_2D_Temp(end,:))/2;
    Ground_2D=repmat(Ground,[2040 1]);
    %Data1_2D_Temp=importdata(sprintf('131111_onion_inter_%i',Array_Number(p)));
    Data1_2D_New=interp1(Calibration_Array,Data1_2D_Temp-Ground_2D,fx);
    Data1_2D_New(isnan(Data1_2D_New))=0;

    %imagesc(Data1_2D_New);
    %axis equal
    %colormap(gray);
    FFT = abs(fft(Data1_2D_New,[],1));
    FFT=FFT(1:round(size(FFT,1)/2),:);
    FFT(1:LP,:)=0;
    Image=Image+FFT;
    disp(p);
end

Image=Image-min(min(Image));
Image_dB=20*log10(Image);
Image_dB=Image_dB-max(max(Image_dB));
imagesc(Image_dB,'xdata',1*[1:size(Image,2)],'ydata',0.2*[1:size(Image,1)]);
xlabel('X Position (micron)');
ylabel('X Position (micron)');
%axis equal
colormap(gray);
%% SNR
Sig_max=max(Image(:,700));
noise=std(Image(1000:end,700));
SNR=20*log10(Sig_max/noise);

%%
plot(Data1_2D_New(:,450));

%%

Image_single=FFT;
Image_single_dB=20*log10(Image_single);
Image_single_dB=Image_single_dB-max(max(Image_single_dB));
imagesc(Image_single_dB,'xdata',1*[1:size(Image,2)],'ydata',0.2*[1:size(Image,1)]);
%axis equal

%plot(real(FFT_1(:,500)))
c=3E8;
Time_total=1/(max(fx)/(length(fx)-1));
Time=[0:Time_total/(length(fx)-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=c*Time;
Position_micron=Position*1E6;
d_Position_micron=Position_micron(2)-Position_micron(1);