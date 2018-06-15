clear all

%Too bad axial resolution
fitting=0;
new_bit_depth=16;

Sample_Path='D:\Users\TuanShu\130314_LCDCrack_01\';
Sample_Path_Bscan='D:\Users\TuanShu\130314_LCDCrack_01\Enface\';

set_peak_index=250;
Averaging_Factor=1;
Sampling_Rate=381.388/Averaging_Factor;
PZT_Speed=0.895972828129950;         %micron/sec

TD_OCT_point_Spacing=PZT_Speed/Sampling_Rate;   
Binning_Factor=1;

Objective_Focal_Length=4/0.85;   %cm
Porjection_Lens_Focal_Length=50;    %cm
Pixel_Size=0.72;       %micron

%Lateral_Spacing=Pixel_Size/Porjection_Lens_Focal_Length*Objective_Focal_Length;
Lateral_Spacing=600/(538-330);

image_index=0:500;

ROI=[1 760;1 1000];      %up, down, left, right

Select_Position=[380 500];

%ROI=[43 342;363 662];      %up, down, left, right
%ROI=[1 384;257 768];      %up, down, left, right
%%ROI=[151 234;457 568];      %up, down, left, right
Image_Stack(1:ceil(length(image_index)/Averaging_Factor),1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;
Low_Pass=30;   %pixel
High_Pass=130;
Image_Temp_2(1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;

        cd(Sample_Path);

    for p=1:length(image_index)
        Image_Temp=imread(sprintf('%i.png',image_index(p)));   

        [m,n]=size(Image_Temp); %M is the original matrix

        Image_Temp=sum(reshape(Image_Temp,Binning_Factor,[]),1);
        Image_Temp=reshape(Image_Temp,m/Binning_Factor,[]).';

        Image_Temp=sum(reshape(Image_Temp,Binning_Factor,[]) ,1);
        Image_Temp=reshape(Image_Temp,n/Binning_Factor,[]).';
        if rem(p,Averaging_Factor)==0
            Image_Temp_2=Image_Temp(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2))+Image_Temp_2;
            Image_Stack(ceil(p/Averaging_Factor),:,:)=Image_Temp_2;
            Image_Temp_2=0;
        else
            Image_Temp_2=Image_Temp(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2))+Image_Temp_2;
        end
        %Image_Stack(ceil(p/Averaging_Factor),:,:)=Image_Temp;
        %Image=imread(sprintf('1.jpg',image_index(p)));    
        disp(p);
    end
    Image_Stack(size(Image_Stack,1),:,:)=Image_Stack(size(Image_Stack,1)-1,:,:);
    Position_micron=(1:size(Image_Stack,1))*TD_OCT_point_Spacing;

    Start_DC=Image_Stack(1,:,:);
    End_DC=Image_Stack(end,:,:);
    DC_volume(1:size(Image_Stack,1),1:size(Image_Stack,2),1:size(Image_Stack,3))=0;
    for q=1:size(Image_Stack,1)
        DC_volume(q,:,:)=Start_DC+(End_DC-Start_DC)*q/size(Image_Stack,1);
    end
    %Image_Stack=Image_Stack-DC_volume;
    plot(real(Image_Stack(:,500,500)));
    clear DC_volume
    [Max_value Max_index]=max(Image_Stack(:,120,160));
    Image_Stack=circshift(Image_Stack,set_peak_index-Max_index);

    FFT_Stack=fft(Image_Stack,[],1);
    FFT_Stack(1:Low_Pass,:,:)=0;
    FFT_Stack(round(size(FFT_Stack,1)/2):end,:,:)=0;
    FFT_Stack(High_Pass:end,:,:)=0;

    Image_Stack_New=(ifft(FFT_Stack,[],1));
    %Max_Value_of_The_Volume=max(abs(Image_Stack_New(:)));
    %Image_Stack_New_Normalized=Image_Stack_New./max(abs(Image_Stack_New(:)));
    %Image_Stack_New_Env=log10(abs(Image_Stack_New)/max(abs(Image_Stack_New(:)))*10);
    Image_Stack_New_Env=(abs(Image_Stack_New)/max(abs(Image_Stack_New(:))));
    [maxvalue maxindex]=max(Image_Stack_New_Env);
    Height_Map(:,:)=(maxindex-maxindex(1,1,1))*TD_OCT_point_Spacing;
    %Image_Stack_New_Env=Image_Stack_New_Env/max(abs(Image_Stack_New_Env(:)));
        cd(Sample_Path_Bscan);
        %dlmwrite('Max_Value_Sam.txt',Max_Value_of_The_Volume,'delimiter','\t','newline','pc');
        for rr=1:size(Image_Stack_New_Env,1)
            Temp_Bscan(:,:)=(Image_Stack_New_Env(rr,:,:));
            imwrite(Temp_Bscan,sprintf('Enface_Sam_%i.bmp',rr));%,'Bitdepth',new_bit_depth);
            %dlmwrite(sprintf('Bscan_Sam_%i.txt',rr),Temp_Bscan,'delimiter','\t','newline','pc');
            disp(rr);
        end
imagesc(Height_Map,'xdata',Pixel_Size:Pixel_Size:Pixel_Size*320,'ydata',Pixel_Size:Pixel_Size:Pixel_Size*240);
xlabel('(micron)');
ylabel('(micron)');
axis equal

plot(Pixel_Size:Pixel_Size:Pixel_Size*240,Height_Map(:,round(100/Pixel_Size)));

%%
%plot(Position_micron,Image_Stack_New_Env(:,120,160));
Frame(:,:)=abs(Image_Stack_New(170,:,:));
imagesc(Frame);

