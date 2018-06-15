clear all

%Too bad axial resolution
fitting=0;
new_bit_depth=8;
cd('D:\Users\TuanShu\130114\Yellow_Sam_0.004_1035.197_Upsidedown\');
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
%ROI=[43 342;363 662];      %up, down, left, right
%ROI=[1 384;257 768];      %up, down, left, right
%%ROI=[151 234;457 568];      %up, down, left, right
Image_Stack(1:ceil(length(image_index)/Averaging_Factor),1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;
Low_Pass=50;   %pixel
High_Pass=400;
Image_Temp_2(1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;
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


[Max_value Max_index]=max(Image_Stack(:,30,40));

Image_Stack=circshift(Image_Stack,set_peak_index-Max_index);



FFT_Stack=fft(Image_Stack,[],1);

FFT_Stack(1:Low_Pass,:,:)=0;
FFT_Stack(round(size(FFT_Stack,1)/2):end,:,:)=0;
FFT_Stack(High_Pass:end,:,:)=0;

%clear Image_Stack Image Image_Temp
%FFT_Stack(round(size(FFT_Stack,1)/2):end,:,:)=0;
Image_Stack_New=(ifft(FFT_Stack,[],1));

posindex=[10 20 30 40 50 60 70 80];

hold all
for qq=1:length(posindex)
    plot(Position_micron(Position_micron<20),real(Image_Stack_New(Position_micron<20,30,posindex(qq))));
end
hold off
xlabel('Optical Path Difference (micron)');
ylabel('Amplitude (a.u.)');
%[Max_value Max_index]=max(Image_Stack_New(:,30,40));

%Image_Stack_New=circshift(Image_Stack_New,set_peak_index-Max_index);

%plot(Position_micron,real(Image_Stack(:,30,40)));

New_Max=max(Image_Stack_New(:));

%Image_Stack_New_Normalized=Image_Stack_New./max(Image_Stack_New(:));



%

Sam=Image_Stack_New(:,30,40);
Sam_2=Image_Stack_New(:,10,40);
Sam_real=real(Sam);
Sam_2_real=real(Sam_2);
%plot(1:length(Sam),Sam,1:length(Ref),Ref);
%dlmwrite('Sam.txt',Sam,'delimiter','\t','newline','pc');

%dlmwrite('Sam_2.txt',Sam_2,'delimiter','\t','newline','pc');

cd('D:\Users\TuanShu\130114\');

%dlmwrite('Sam_real.txt',Sam_real,'delimiter','\t','newline','pc');

%dlmwrite('Sam_2_real.txt',Sam_2_real,'delimiter','\t','newline','pc');


dlmwrite('Ref_real.txt',Sam_real,'delimiter','\t','newline','pc');

dlmwrite('Ref_2_real.txt',Sam_2_real,'delimiter','\t','newline','pc');
%for rr=1:size(Image_Stack_New_Normalized,3)
%    Temp_Bscan(:,:)=Image_Stack_New_Normalized(:,:,rr);
%    imwrite(Temp_Bscan,sprintf('Bscan_%i.png',rr),'Bitdepth',new_bit_depth);
%    disp(rr);
%end

%for q=1:size(Image_Stack,1)
%    Temp_Slides(:,:)=Image_Stack_New_Normalized(q,:,:);
%    imwrite(Temp_Slides,sprintf('Enface_Binned_by_%i_%i.png',Binning_Factor,image_index(q)),'Bitdepth',new_bit_depth);
%    disp(q);
%end
%Temp_Sides_Read=imread(sprintf('%i_New.png',image_index(q)));








