clear all

%Too bad axial resolution
fitting=0;
new_bit_depth=16;

Sample_Path='D:\Users\TuanShu\130129\BinnedSam\';
Sample_Path_Bscan='D:\Users\TuanShu\130129\BinnedSam\Bscan\';
Reference_Path='D:\Users\TuanShu\130129\BinnedRef\';
Reference_Path_Bscan='D:\Users\TuanShu\130129\BinnedRef\\Bscan\';

set_peak_index=1500;
Averaging_Factor=10;
Sampling_Rate=1035.197/10;
PZT_Speed=0.895972828129950;         %micron/sec

TD_OCT_point_Spacing=PZT_Speed/Sampling_Rate;   

Binning_Factor=1;

Stage_speed=1.5;  %micron/sec
Sampling_rate=28;   %Hz

Frame_Axial_Spacing=Stage_speed/Sampling_rate*5/3;  %micron

Objective_Focal_Length=4/0.85;   %cm
Porjection_Lens_Focal_Length=50;    %cm
Pixel_Size=7.4;       %micron

%Lateral_Spacing=Pixel_Size/Porjection_Lens_Focal_Length*Objective_Focal_Length;
Lateral_Spacing=600/(538-330);

image_index=0:79999;
ROI=[1 61;1 80];      %up, down, left, right

Select_Position=[30 40];

%ROI=[43 342;363 662];      %up, down, left, right
%ROI=[1 384;257 768];      %up, down, left, right
%%ROI=[151 234;457 568];      %up, down, left, right
Image_Stack(1:ceil(length(image_index)/Averaging_Factor),1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;
Low_Pass=50;   %pixel
High_Pass=400;
Image_Temp_2(1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;

for qqq=1:2
    if qqq==1
        cd(Sample_Path);
    else
        cd(Reference_Path);
    end
    for p=1:length(image_index)
        Image_Temp=imread(sprintf('%06i.tif',image_index(p)));   

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
    Image_Stack=Image_Stack-DC_volume;
    clear DC_volume
    [Max_value Max_index]=max(Image_Stack(:,30,40));
    Image_Stack=circshift(Image_Stack,set_peak_index-Max_index);

    FFT_Stack=fft(Image_Stack,[],1);
    FFT_Stack(1:Low_Pass,:,:)=0;
    FFT_Stack(round(size(FFT_Stack,1)/2):end,:,:)=0;
    FFT_Stack(High_Pass:end,:,:)=0;

    Image_Stack_New=(ifft(FFT_Stack,[],1));
    %Max_Value_of_The_Volume=max(abs(Image_Stack_New(:)));
    %Image_Stack_New_Normalized=Image_Stack_New./max(abs(Image_Stack_New(:)));
    Image_Stack_New_Real=real(Image_Stack_New);
    if qqq==1
        cd(Sample_Path_Bscan);
        %dlmwrite('Max_Value_Sam.txt',Max_Value_of_The_Volume,'delimiter','\t','newline','pc');
        for rr=1:size(Image_Stack_New_Real,3)
            Temp_Bscan(:,:)=Image_Stack_New_Real(:,:,rr);
            %imwrite(Temp_Bscan,sprintf('Bscan_Sam_%i.png',rr),'Bitdepth',new_bit_depth);
            dlmwrite(sprintf('Bscan_Sam_%i.txt',rr),Temp_Bscan,'delimiter','\t','newline','pc');
            disp(rr);
        end
    else
        cd(Reference_Path_Bscan);
        %dlmwrite('Max_Value_Ref.txt',Max_Value_of_The_Volume,'delimiter','\t','newline','pc');
        for rr=1:size(Image_Stack_New_Real,3)
            Temp_Bscan(:,:)=Image_Stack_New_Real(:,:,rr);
            %imwrite(Temp_Bscan,sprintf('Bscan_Ref_%i.png',rr),'Bitdepth',new_bit_depth);
            dlmwrite(sprintf('Bscan_Ref_%i.txt',rr),Temp_Bscan,'delimiter','\t','newline','pc');
            disp(rr);
        end
    end

end

%%
plot(Image_Stack_New(:,30,40));
Frame=Image_Stack_New(1000,:,:);
imagesc(Frame);

