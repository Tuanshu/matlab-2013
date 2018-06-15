clear all

%Too bad axial resolution
fitting=1;
new_bit_depth=8;
cd('D:\Users\TuanShu\');

Binning_Factor=4;

Stage_speed=1.5;  %micron/sec
Sampling_rate=28;   %Hz

Frame_Axial_Spacing=Stage_speed/Sampling_rate;  %micron

Objective_Focal_Length=4/0.85;   %cm
Porjection_Lens_Focal_Length=50;    %cm
Pixel_Size=7.4;       %micron

Lateral_Spacing=1.8*Binning_Factor;


%% Height generation
clear Xgrid Ygrid
N=2000;
Xgrid(1:N,1:N)=0;
Ygrid(1:N,1:N)=0;

K_given=1/7800;
conic=0.5;
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

Height_Sim=(K_given.*R.^2)./(1+(1-(1+conic).*(K_given.^2).*(R.^2)).^0.5);

%Height_Sim=R_given-(R_given^2-R.^2).^0.5;

Height_Sim(abs(imag(Height_Sim))>0)=0;
imagesc(Height_Sim)

%%

[K,H,P1,P2] = surfature(Xgrid,Ygrid,Height_Sim); 

R1=1./real(P1);
R2=1./real(P2);
RK=1./real(K);
RH=1./real(H);


imagesc(R1)
caxis([7000 8000]);

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