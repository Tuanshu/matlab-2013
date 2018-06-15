clear all
loaddata=0;
if loaddata==1
    cd('D:\Users\TuanShu\121106_Lens2\');

    Stage_speed=1;  %micron/sec
    Sampling_rate=33;   %Hz

    Frame_Axial_Spacing=Stage_speed/Sampling_rate*5/3;  %micron

    Objective_Focal_Length=4/0.85;   %cm
    Porjection_Lens_Focal_Length=50;    %cm
    Pixel_Size=7.4;       %micron

    %Lateral_Spacing=Pixel_Size/Porjection_Lens_Focal_Length*Objective_Focal_Length;
    Lateral_Spacing=600/(538-330);

    image_index=1:900;
    ROI=[1 768;1 1024];      %up, down, left, right
    Image_Stack(length(image_index),1:(ROI(1,2)-ROI(1,1)+1),1:(ROI(2,2)-ROI(2,1)+1))=0;
    Low_Pass=20;   %pixel
    High_Pass=400;
    for p=1:length(image_index)
        Image=imread(sprintf('%i.png',image_index(p)));    
        %Image=imread(sprintf('1.jpg',image_index(p)));    
        Image_Stack(p,:,:)=Image(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),1);
        disp(p);
    end

    FFT_Stack=fft(Image_Stack,[],1);
    FFT_Stack(1:Low_Pass,:,:)=0;
    FFT_Stack(round(length(FFT_Stack)/2):end,:,:)=0;
    Image_Stack_New=real(ifft(FFT_Stack,[],1));
    %plot(real(FFT_Stack(:,300,400)));

    [max_value max_index]=max(Image_Stack_New,[],1);

    MaxValue(:,:)=max_value(1,:,:);

    %plot(real(FFT_Stack(:,300,400)));

    Height(:,:)=Frame_Axial_Spacing*max_index(1,:,:);
    Height=max(max(Height))-Height;

    imagesc(Height,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Height,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Height,1)-1));
    xlabel('(Micron)');
    ylabel('(Micron)');
    colorbar
    caxis([0 45])
    axis equal
    axis([0 Lateral_Spacing*(size(Height,2)-1) 0 Lateral_Spacing*(size(Height,1)-1)])
end

%% Height generation
cd('D:\Users\TuanShu\');
Lateral_Spacing_Sim=10;
N=2000;
Xgrid_Sim(1:N,1:N)=0;
Ygrid_Sim(1:N,1:N)=0;

K_given=1/7800;
conic=0.7;
R_given=7800;
R_given_2=10000;

for p=1:size(Xgrid_Sim,2)
    Xgrid_Sim(:,p)=(p-1)*Lateral_Spacing_Sim;
end
for q=1:size(Ygrid_Sim,1)
    Ygrid_Sim(q,:)=(q-1)*Lateral_Spacing_Sim;
end


X_Center=N/2*Lateral_Spacing_Sim;
Y_Center=N/2*Lateral_Spacing_Sim;


R=((Xgrid_Sim-X_Center).^2+(Ygrid_Sim-Y_Center).^2).^0.5;

%Paraaxial: Height_Sim=1./K_given-(K_given.*R.^2)./(1+(1-(1+conic).*(K_given.^2).*(R.^2)).^0.5);

%Height_Sim=(2*R_given*(Xgrid_Sim-X_Center+R_given)-(1+conic)*(Xgrid_Sim-X_Center+R_given).^2+2*R_given*(Ygrid_Sim-Y_Center+R_given)-(1+conic)*(Ygrid_Sim-Y_Center+R_given).^2-R_given^2).^0.5;
%Height_Sim=max(max(Height_Sim))-Height_Sim;
%Height_Sim=(R_given^2-R.^2).^0.5;
Height_Sim=R_given.*(1-((Xgrid_Sim-X_Center)./R_given).^2-((Ygrid_Sim-Y_Center)./R_given_2).^2).^0.5;

Height_Sim(abs(imag(Height_Sim))>0)=0;
imagesc(Height_Sim)
axis equal
%%

Binning_Factor_Height=1;
Smooth_Factor_Height=1;

%Pixel_Size_New=Lateral_Spacing*Binning_Factor_Height;

Height_Temp=Height_Sim;

Filter2D=zeros(Smooth_Factor_Height,Smooth_Factor_Height);

X_Filter2D=zeros(Smooth_Factor_Height,Smooth_Factor_Height);
Y_Filter2D=zeros(Smooth_Factor_Height,Smooth_Factor_Height);

XC_Filter2D=size(Filter2D,1)/2;

YC_Filter2D=size(Filter2D,1)/2;


for p=1:size(Filter2D,2)
    X_Filter2D(:,p)=p;
end
for q=1:size(Filter2D,1)
    Y_Filter2D(q,:)=q;
end

Filter2D(((X_Filter2D-XC_Filter2D).^2+(Y_Filter2D-YC_Filter2D).^2)<XC_Filter2D^2)=1;

%Height_Temp=conv2(Height_Temp,Filter2D,'valid')/Smooth_Factor_Height/Smooth_Factor_Height;

[m,n]=size(Height_Temp);

Height_Temp=sum(reshape(Height_Temp,Binning_Factor_Height,[]),1);
Height_Temp=reshape(Height_Temp,m/Binning_Factor_Height,[]).';

Height_Temp=sum( reshape(Height_Temp,Binning_Factor_Height,[]) ,1);
Height_Temp=reshape(Height_Temp,n/Binning_Factor_Height,[]).';


%Height_Temp=Height_Temp/Binning_Factor_Height^2;

Height_Binned=Height_Temp;

%Height_Binned=Height_Front;

imagesc(Height_Binned,'xdata',0:Lateral_Spacing_Sim:Lateral_Spacing_Sim*(size(Height_Binned,2)-1),'ydata',0:Lateral_Spacing_Sim:Lateral_Spacing_Sim*(size(Height_Binned,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
colorbar
axis equal
axis([0 Lateral_Spacing_Sim*(size(Height_Binned,2)-1) 0 Lateral_Spacing_Sim*(size(Height_Binned,1)-1)])
%caxis([7000 8000])



%%
clear Xgrid Ygrid
Xgrid(1:size(Height_Binned,1),1:size(Height_Binned,2))=0;
Ygrid(1:size(Height_Binned,1),1:size(Height_Binned,2))=0;


for p=1:size(Xgrid,2)
    Xgrid(:,p)=(p-1)*Lateral_Spacing_Sim;
end
for q=1:size(Ygrid,1)
    Ygrid(q,:)=(q-1)*Lateral_Spacing_Sim;
end


%%

[K,H,P1,P2] = surfature(Xgrid,Ygrid,Height_Binned); 
X=Xgrid;Y=Y=grid;Z=Height_Binned;

R1=1./real(P1)/1000;
R2=1./real(P2)/1000;
RK=1./real(K)/1000/1000;
RH=1./real(H)/1000;



imagesc(-R1,'xdata',0:Lateral_Spacing_Sim:Lateral_Spacing_Sim*(size(R1,2)-1),'ydata',0:Lateral_Spacing_Sim:Lateral_Spacing_Sim*(size(R1,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
colorbar
axis equal
axis([0 Lateral_Spacing_Sim*(size(R1,2)-1) 0 Lateral_Spacing_Sim*(size(R1,1)-1)])
caxis([5 15])
%caxis([0 20])

%% trans to polar


Center=[size(R1,1)/2+200 size(R1,2)/2+200];

        
       
[r c] = size(R1);
[X Y] = meshgrid(1:c,1:r);
        
NNN=1000;

r_wish=repmat([1:1000/NNN:1000],NNN,1);
theta_wish=repmat([2*pi/NNN:2*pi/NNN:2*pi]',1,NNN);
        
X_wish=r_wish.*cos(theta_wish)+Center(2);
Y_wish=r_wish.*sin(theta_wish)+Center(1);

R1_Rtheta=interp2(X,Y,double(R1),X_wish,Y_wish);
        

imagesc(-R1_Rtheta,'xdata',0:Lateral_Spacing_Sim:Lateral_Spacing_Sim*(size(R1,2)-1),'ydata',0:Lateral_Spacing_Sim:Lateral_Spacing_Sim*(size(R1,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
colorbar
axis equal
axis([0 Lateral_Spacing_Sim*(size(R1,2)-1) 0 Lateral_Spacing_Sim*(size(R1,1)-1)])
caxis([5 15])


%% again back to XY


Center_2=[size(R1,1)/2-200 size(R1,2)/2-200];

[r_rt c_rt] = size(R1_Rtheta);
[X_rt Y_rt] = meshgrid(1:c_rt,1:r_rt);
        
r_XY=((X-Center_2(2)).^2+(Y-Center_2(1)).^2).^0.5;
theta_XY=atan2((Y-Center_2(1)),(X-Center_2(2)))+pi;
        
R1_XY_New=interp2(r_wish,theta_wish,R1_Rtheta,r_XY,theta_XY);
R1_XY_New(isnan(R1_XY_New))=0;

imagesc(-R1_XY_New,'xdata',0:Lateral_Spacing_Sim:Lateral_Spacing_Sim*(size(R1,2)-1),'ydata',0:Lateral_Spacing_Sim:Lateral_Spacing_Sim*(size(R1,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
colorbar
axis equal
axis([0 Lateral_Spacing_Sim*(size(R1,2)-1) 0 Lateral_Spacing_Sim*(size(R1,1)-1)])
caxis([5 15])
%%

Smooth_Factor_Curvature=1;

Filter2D_Curvature=zeros(Smooth_Factor_Curvature,Smooth_Factor_Curvature);

X_Filter2D_Curvature=zeros(Smooth_Factor_Curvature,Smooth_Factor_Curvature);
Y_Filter2D_Curvature=zeros(Smooth_Factor_Curvature,Smooth_Factor_Curvature);

XC_Filter2D_Curvature=size(Filter2D_Curvature,1)/2;

YC_Filter2D_Curvature=size(Filter2D_Curvature,1)/2;


for p=1:size(Filter2D_Curvature,2)
    X_Filter2D_Curvature(:,p)=p;
end
for q=1:size(Filter2D_Curvature,1)
    Y_Filter2D_Curvature(q,:)=q;
end

Filter2D_Curvature(((X_Filter2D_Curvature-XC_Filter2D_Curvature).^2+(Y_Filter2D_Curvature-YC_Filter2D_Curvature).^2)<XC_Filter2D_Curvature^2)=1;

P1_Smoothed=conv2(P1,ones(Smooth_Factor_Curvature,Smooth_Factor_Curvature),'valid')/Smooth_Factor_Curvature/Smooth_Factor_Curvature;
R1_Smoothed=1./P1_Smoothed;
imagesc(R1_Smoothed)
caxis([0 20000]);


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