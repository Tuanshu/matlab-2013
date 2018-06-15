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


for q=1:length(image_index)
    Temp_Slides(:,:)=Image_Stack_New_Normalized(q,:,:);
    imwrite(Temp_Slides,sprintf('Enface_Binned_by_%i_%i.png',Binning_Factor,image_index(q)),'Bitdepth',new_bit_depth);
    disp(q);
end


for rr=1:size(Image_Stack_New_Normalized,3)
    Temp_Bscan(:,:)=Image_Stack_New_Normalized(:,:,rr);
    imwrite(Temp_Bscan,sprintf('Bscan_%i.png',rr),'Bitdepth',new_bit_depth);
    disp(rr);
end

%Temp_Sides_Read=imread(sprintf('%i_New.png',image_index(q)));





Lower_Limit_of_Front=2500;

Upper_Limit_of_Rear=10000;

Image_Stack_Front=Image_Stack_New;
Image_Stack_Rear=Image_Stack_New;
Image_Stack_Front(Lower_Limit_of_Front:end,:,:)=0;
Image_Stack_Rear(1:Upper_Limit_of_Rear,:,:)=0;

if fitting ==1
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

[maxvalue1D maxindex1D]=max(Height_Front(:));
Xgrid(1:size(Height_Front,1),1:size(Height_Front,2))=0;
Ygrid(1:size(Height_Front,1),1:size(Height_Front,2))=0;

for p=1:size(Height_Front,2)
    Xgrid(:,p)=(p-1)*Lateral_Spacing;
end
for q=1:size(Height_Front,1)
    Ygrid(q,:)=(q-1)*Lateral_Spacing;
end

FitSurface= fittype( @(c,r,a, b, x, y) c+(r^2-(x-a).^2-(y-b).^2).^0.5, 'independent', {'x', 'y'},'dependent', 'z');    %(x-A)^2+(y-B)^2+(z-C)^2=R^2
                                                                                                                        %z=C+(R^2-(x-A)^.2-(y-B)^.2).^0.5
                                                                                                                        
Starting_X_Front=128*Lateral_Spacing;                                                                                                                                                                                                                              
Starting_Y_Front=48*Lateral_Spacing;
R_Weight=40*Lateral_Spacing;
R_Weight_Inner=0;
Weight_Func_Front=MaxValue/max(max(MaxValue));
Weight_Func_Front(:,:)=1;
Weight_Func_Front((Xgrid-Starting_X_Front).^2+(Ygrid-Starting_Y_Front).^2>R_Weight^2)=0;
Weight_Func_Front((Xgrid-Starting_X_Front).^2+(Ygrid-Starting_Y_Front).^2<R_Weight_Inner^2)=0;
FitPara_Front=fit([Xgrid(:),Ygrid(:)],Height_Front(:),FitSurface,'Weight',Weight_Func_Front(:),'StartPoint',[-15000,16000,Starting_X_Front,Starting_Y_Front]);
CheckCurve_Front=FitPara_Front.c+((FitPara_Front.r)^2-(Xgrid-FitPara_Front.a).^2-(Ygrid-FitPara_Front.b).^2).^0.5;

      
Starting_X_Rear=128*Lateral_Spacing;                                                                                                                                                                                                                              
Starting_Y_Rear=120;%48*Lateral_Spacing;
R_Weight=40*Lateral_Spacing;
R_Weight_Inner=0;
Weight_Func_Rear=MaxValue/max(max(MaxValue));
Weight_Func_Rear(:,:)=1;
Weight_Func_Rear((Xgrid-Starting_X_Rear).^2+(Ygrid-Starting_Y_Rear).^2>R_Weight^2)=0;
Weight_Func_Rear((Xgrid-Starting_X_Rear).^2+(Ygrid-Starting_Y_Rear).^2<R_Weight_Inner^2)=0;
FitPara_Rear=fit([Xgrid(:),Ygrid(:)],Height_Rear(:),FitSurface,'Weight',Weight_Func_Rear(:),'StartPoint',[-15000,16000,Starting_X_Rear,Starting_Y_Rear]);
CheckCurve_Rear=FitPara_Rear.c+((FitPara_Rear.r)^2-(Xgrid-FitPara_Rear.a).^2-(Ygrid-FitPara_Rear.b).^2).^0.5;



set(gcf,'Position',[100 100 700 700*Lateral_Spacing*(size(Height_Front,1)-1)/(Lateral_Spacing*(size(Height_Front,2)-1))])             %0可能是螢幕, gca是目前圖形的handle (如果沒有圖的話會自己開新的一個handle)

imagesc(MaxValue_Front/max(max(MaxValue_Front)),'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(MaxValue_Front,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(MaxValue_Front,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
colorbar
axis equal
axis([0 Lateral_Spacing*(size(Height_Front,2)-1) 0 Lateral_Spacing*(size(Height_Front,1)-1)])

imagesc(MaxValue_Rear/max(max(MaxValue_Rear)),'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(MaxValue_Rear,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(MaxValue_Rear,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
colorbar
axis equal
axis([0 Lateral_Spacing*(size(Height_Front,2)-1) 0 Lateral_Spacing*(size(Height_Front,1)-1)])


imagesc(Weight_Func_Front,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Weight_Func_Front,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Weight_Func_Front,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
axis equal
axis([0 Lateral_Spacing*(size(Weight_Func_Front,2)-1) 0 Lateral_Spacing*(size(Weight_Func_Front,1)-1)])



imagesc(Weight_Func_Rear,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Weight_Func_Rear,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Weight_Func_Rear,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
axis equal
axis([0 Lateral_Spacing*(size(Weight_Func_Rear,2)-1) 0 Lateral_Spacing*(size(Weight_Func_Rear,1)-1)])

imagesc(CheckCurve_Front,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(CheckCurve_Front,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(CheckCurve_Front,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
colorbar
caxis([0 60])
axis equal
axis([0 Lateral_Spacing*(size(CheckCurve_Front,2)-1) 0 Lateral_Spacing*(size(CheckCurve_Front,1)-1)])


imagesc(CheckCurve_Rear,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(CheckCurve_Rear,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(CheckCurve_Rear,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
colorbar
caxis([0 60])
axis equal
axis([0 Lateral_Spacing*(size(CheckCurve_Rear,2)-1) 0 Lateral_Spacing*(size(CheckCurve_Rear,1)-1)])

%plot((1:(size(Image_Stack_New,1)))*Frame_Axial_Spacing,real(Image_Stack_New(:,300,400)));
%xlabel('Axial Position (Micron)');
%ylabel('Amplitude (a.u.)');

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


imagesc(Thickness,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Thickness,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Thickness,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
colorbar
caxis([0 1000])
axis equal
axis([0 Lateral_Spacing*(size(Thickness,2)-1) 0 Lateral_Spacing*(size(Thickness,1)-1)])

plot(0:Lateral_Spacing:Lateral_Spacing*(size(Height_Front,2)-1),Height_Front(48,:),0:Lateral_Spacing:Lateral_Spacing*(size(Height_Front,2)-1),CheckCurve_Front(48,:));
xlabel('(Micron)');
ylabel('(Micron)');
legend('Measured height','Fitting curve')

Fitted_Curvature_Front=FitPara_Front.r %(micron)

Fitted_Curvature_Rear=FitPara_Rear.r %(micron)
end