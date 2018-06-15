clear all

%Too bad axial resolution
fitting=1;
new_bit_depth=8;
cd('D:\Users\TuanShu\121205_Eye by 830-SLD2\');

Binning_Factor=4;

Stage_speed=1.5;  %micron/sec
Sampling_rate=28;   %Hz

Frame_Axial_Spacing=Stage_speed/Sampling_rate;  %micron

Objective_Focal_Length=4/0.85;   %cm
Porjection_Lens_Focal_Length=50;    %cm
Pixel_Size=7.4;       %micron

Lateral_Spacing=1.8*Binning_Factor;

image_index=0:16500;
ROI=[1 384;1 1024];      %up, down, left, right
%ROI=[43 342;363 662];      %up, down, left, right
%ROI=[1 384;257 768];      %up, down, left, right
%%ROI=[151 234;457 568];      %up, down, left, right
Image_Stack(length(image_index),1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;
Low_Pass=1500;   %pixel
High_Pass=2500;
for p=1:length(image_index)
    Image_Temp=imread(sprintf('%i.png',image_index(p)));   
    
    [m,n]=size(Image_Temp); %M is the original matrix

    Image_Temp=sum( reshape(Image_Temp,Binning_Factor,[]) ,1 );
    Image_Temp=reshape(Image_Temp,m/Binning_Factor,[]).';

    Image_Temp=sum( reshape(Image_Temp,Binning_Factor,[]) ,1);
    Image_Temp=reshape(Image_Temp,n/Binning_Factor,[]).';
    
    %Image=imread(sprintf('1.jpg',image_index(p)));    
    Image_Stack(p,:,:)=Image_Temp;%(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),1);
    disp(p);
end


FFT_Stack=fft(Image_Stack,[],1);

FFT_Stack(1:Low_Pass,:,:)=0;
FFT_Stack(High_Pass:end,:,:)=0;

clear Image_Stack Image Image_Temp
%FFT_Stack(round(size(FFT_Stack,1)/2):end,:,:)=0;
Image_Stack_New=abs(ifft(FFT_Stack,[],1));
%plot(real(FFT_Stack(:,500,400)));
New_Max=max(Image_Stack_New(:));

Image_Stack_New_Normalized=Image_Stack_New./max(Image_Stack_New(:));

%Temp_Sides_Read=imread(sprintf('%i_New.png',image_index(q)));


Lower_Limit_of_Front=2500;

Upper_Limit_of_Rear=10000;

Image_Stack_Front=Image_Stack_New;
Image_Stack_Rear=Image_Stack_New;
Image_Stack_Front(Lower_Limit_of_Front:end,:,:)=0;
Image_Stack_Rear(1:Upper_Limit_of_Rear,:,:)=0;

%plot(real(FFT_Stack(:,300,400)));

Frame_Axial_Spacing=Stage_speed/Sampling_rate;  %micron
Lateral_Spacing=1.8*Binning_Factor;

[max_value_front max_index_front]=max(Image_Stack_Front,[],1);
MaxValue_Front(:,:)=max_value_front(1,:,:);
Height_Front(:,:)=Frame_Axial_Spacing*max_index_front(1,:,:);

[max_value_rear max_index_rear]=max(Image_Stack_Rear,[],1);
MaxValue_Rear(:,:)=max_value_rear(1,:,:);
Height_Rear(:,:)=Frame_Axial_Spacing*max_index_rear(1,:,:);

Thickness=Height_Rear-Height_Front;
Height_Front=max(max(Height_Front))-Height_Front;
Height_Rear=max(max(Height_Rear))-Height_Rear;

imagesc(Height_Front,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Height_Front,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Height_Front,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
colorbar
caxis([0 60])
axis equal
axis([0 Lateral_Spacing*(size(Height_Front,2)-1) 0 Lateral_Spacing*(size(Height_Front,1)-1)])

imagesc(Height_Rear,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Height_Rear,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Height_Rear,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
colorbar
caxis([0 60])
axis equal
axis([0 Lateral_Spacing*(size(Height_Rear,2)-1) 0 Lateral_Spacing*(size(Height_Rear,1)-1)])

%% Height generation
clear Xgrid Ygrid
N=2000;
Xgrid(1:N,1:N)=0;
Ygrid(1:N,1:N)=0;

K_given=1/10000;
conic=1;
R_given=7800;

for p=1:size(Xgrid,2)
    Xgrid(:,p)=(p-1)*Lateral_Spacing;
end
for q=1:size(Ygrid,1)
    Ygrid(q,:)=(q-1)*Lateral_Spacing;
end


X_Center=N/2*Lateral_Spacing;
Y_Center=N/2*Lateral_Spacing;


R=((Xgrid-X_Center).^2+(Ygrid-Y_Center).^2).^0.5;

%Height_Sim=(K_given.*R.^2)./(1+(1-(1+conic).*(K_given.^2).*(R.^2)).^0.5);

Height_Sim=R_given-(R_given^2-R.^2).^0.5;

Height_Sim(abs(imag(Height_Sim))>0)=0;
imagesc(Height_Sim)
%%

Binning_Factor_Height=1;
Smooth_Factor_Height=1;

Pixel_Size_New=Lateral_Spacing*Binning_Factor_Height;

Height_Temp=Height_Sim;

Height_Temp=conv2(Height_Temp,ones(Smooth_Factor_Height,Smooth_Factor_Height),'same')/Smooth_Factor_Height/Smooth_Factor_Height;

[m,n]=size(Height_Temp);

Height_Temp=sum(reshape(Height_Temp,Binning_Factor_Height,[]),1);
Height_Temp=reshape(Height_Temp,m/Binning_Factor_Height,[]).';

Height_Temp=sum( reshape(Height_Temp,Binning_Factor_Height,[]) ,1);
Height_Temp=reshape(Height_Temp,n/Binning_Factor_Height,[]).';


%Height_Temp=Height_Temp/Binning_Factor_Height^2;

Height_Binned=Height_Temp;

Height_Binned=Height_Sim;
imagesc(Height_Binned);
%%

dz_dx=diff(Height_Binned,1,1)/Pixel_Size_New;
dz_dx(size(dz_dx,1)+1,:)=dz_dx(size(dz_dx,1),:);
d2z_dx2=diff(Height_Binned,2,1)/Pixel_Size_New^2;
d2z_dx2=diff(diff(Height_Binned,1,1),1,1)/Pixel_Size_New^2;

d2z_dx2(size(d2z_dx2,1)+1,:)=d2z_dx2(size(d2z_dx2,1),:);
d2z_dx2(size(d2z_dx2,1)+1,:)=d2z_dx2(size(d2z_dx2,1),:);

dz_dy=diff(Height_Binned,1,2)/Pixel_Size_New;
dz_dy(:,size(dz_dy,2)+1)=dz_dy(:,size(dz_dy,2));

d2z_dy2=diff(Height_Binned,2,2)/Pixel_Size_New^2;
d2z_dy2(:,size(d2z_dy2,2)+1)=d2z_dy2(:,size(d2z_dy2,2));
d2z_dy2(:,size(d2z_dy2,2)+1)=d2z_dy2(:,size(d2z_dy2,2));

%%
Kx=d2z_dx2./(1+(dz_dx).^2).^(3/2);
Ky=d2z_dy2./(1+(dz_dy).^2).^(3/2);
Kx_Cut=Kx(3:size(Kx,1)-2,2:size(Kx,2)-2);
Ky_Cut=Ky(3:size(Ky,1)-2,2:size(Ky,2)-2);
%%

Rx=1./Kx_Cut/1000;  %mm
Ry=1./Ky_Cut/1000;  %mm
R_total=1./(1./Rx+1./Ry);
imagesc(Rx,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Rx,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Rx,1)-1));
caxis([0 10]);
axis equal