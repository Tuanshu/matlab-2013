clear all

%% Height generation
cd('D:\Users\TuanShu\');
Lateral_Spacing=50;
NN=250;
X(1:NN,1:NN)=0;
Y(1:NN,1:NN)=0;

Ratio=0.175;

K_given=1/7800;
conic=0.7;
R_given=8500;
R_given_2=7800;

for p=1:size(X,2)
    X(:,p)=(p-1)*Lateral_Spacing;
end
for q=1:size(Y,1)
    Y(q,:)=(q-1)*Lateral_Spacing;
end


X_Center=NN/2*Lateral_Spacing;
Y_Center=NN/2*Lateral_Spacing;

Z=((R_given*R_given_2)^0.5).*(1-((X-X_Center)./R_given).^2-((Y-Y_Center)./R_given_2).^2).^0.5;
Z(abs(imag(Z))>0)=0;

Z=Z-max(max(Z))+Ratio*max(max(Z));

Z(Z<0)=0;

imagesc(Z)
axis equal
%R=((X-X_Center).^2+(Y-Y_Center).^2).^0.5;

%Paraaxial: Z=1./K_given-(K_given.*R.^2)./(1+(1-(1+conic).*(K_given.^2).*(R.^2)).^0.5);

%Z=(2*R_given*(X-X_Center+R_given)-(1+conic)*(X-X_Center+R_given).^2+2*R_given*(Y-Y_Center+R_given)-(1+conic)*(Y-Y_Center+R_given).^2-R_given^2).^0.5;
%Z=max(max(Z))-Z;
%Z=(R_given^2-R.^2).^0.5;

%%

%[K,H,P1,P2] = surfature(X,Y,Z); 
%X=X;Y=Y;Z=Z;

% First Derivatives
[Xu,Xv] = gradient(X);
[Yu,Yv] = gradient(Y);
[Zu,Zv] = gradient(Z);

% Second Derivatives
[Xuu,Xuv] = gradient(Xu);
[Yuu,Yuv] = gradient(Yu);
[Zuu,Zuv] = gradient(Zu);

[Xuv,Xvv] = gradient(Xv);
[Yuv,Yvv] = gradient(Yv);
[Zuv,Zvv] = gradient(Zv);

% Reshape 2D Arrays into Vectors
Xu = Xu(:);   Yu = Yu(:);   Zu = Zu(:); 
Xv = Xv(:);   Yv = Yv(:);   Zv = Zv(:); 
Xuu = Xuu(:); Yuu = Yuu(:); Zuu = Zuu(:); 
Xuv = Xuv(:); Yuv = Yuv(:); Zuv = Zuv(:); 
Xvv = Xvv(:); Yvv = Yvv(:); Zvv = Zvv(:); 

Xu          =   [Xu Yu Zu];
Xv          =   [Xv Yv Zv];
Xuu         =   [Xuu Yuu Zuu];
Xuv         =   [Xuv Yuv Zuv];
Xvv         =   [Xvv Yvv Zvv];

% First fundamental Coeffecients of the surface (E,F,G)
E           =   dot(Xu,Xu,2);
F           =   dot(Xu,Xv,2);
G           =   dot(Xv,Xv,2);

m           =   cross(Xu,Xv,2);
p           =   sqrt(dot(m,m,2));
n           =   m./[p p p]; 

% Second fundamental Coeffecients of the surface (L,M,N)
L           =   dot(Xuu,n,2);
M           =   dot(Xuv,n,2);
N           =   dot(Xvv,n,2);

[s,t] = size(Z);

% Gaussian Curvature
K = (L.*N - M.^2)./(E.*G - F.^2);
%K = reshape(K,s,t);

% Mean Curvature
H = (E.*N + G.*L - 2.*F.*M)./(2*(E.*G - F.^2));
%H = reshape(H,s,t);

% Principal Curvatures
Pmax = H + sqrt(H.^2 - K);
Pmin = H - sqrt(H.^2 - K);

% eigenvalue: Pmax and Pmin, Matrix

%t1= @ (n) (2*n1_Considered)./(n1_Considered+n); 

%
X_P1(1:length(Xu))=0;
Y_P1(1:length(Xu))=0;
Z_P1(1:length(Xu))=0;

X_P2(1:length(Xu))=0;
Y_P2(1:length(Xu))=0;
Z_P2(1:length(Xu))=0;

Curvature_P1(1:length(Xu))=0;
Curvature_P2(1:length(Xu))=0;

for p=1:length(Xu)
    
    I=[E(p) F(p); F(p) G(p)];
    II=[L(p) M(p); M(p) N(p)];
    
    MATRIX=I\II;
    [P_dir P_v]=eig(MATRIX);

    %P_dir_1=P_dir(:,1); %the true principal direction is P_dir_1(1)*Xu+P_dir_1(2)*Xv
    %P_dir_2=P_dir(:,2);
    
    Vector_1=[P_dir(1,1)*Xu(p,1)+P_dir(2,1)*Xv(p,1) P_dir(1,1)*Xu(p,2)+P_dir(2,1)*Xv(p,2) P_dir(1,1)*Xu(p,3)+P_dir(2,1)*Xv(p,3)];
    Vector_2=[P_dir(1,2)*Xu(p)+P_dir(2,2)*Xv(p,1) P_dir(1,2)*Xu(p,2)+P_dir(2,2)*Xv(p,2) P_dir(1,2)*Xu(p,3)+P_dir(2,2)*Xv(p,3)];
    Vector_1=Vector_1/norm(Vector_1);
    Vector_2=Vector_2/norm(Vector_2);
    
    X_P1(p)=Vector_1(1);
    Y_P1(p)=Vector_1(2);
    Z_P1(p)=Vector_1(3);
    
    X_P2(p)=Vector_2(1);
    Y_P2(p)=Vector_2(2);
    Z_P2(p)=Vector_2(3);
    %X_P1(p)=P_dir(1,1)*Xu(p)+P_dir(2,1)*Xv(p);
    %Y_P1(p)=P_dir(1,1)*Yu(p)+P_dir(2,1)*Yv(p);
    %Z_P1(p)=P_dir(1,1)*Zu(p)+P_dir(2,1)*Zv(p);

    %X_P2(p)=P_dir(1,2)*Xu(p)+P_dir(2,2)*Xv(p);
    %Y_P2(p)=P_dir(1,2)*Yu(p)+P_dir(2,2)*Yv(p);
    %Z_P2(p)=P_dir(1,2)*Zu(p)+P_dir(2,2)*Zv(p);
    
    Curvature_P1(p)=-1*P_v(1,1);
    Curvature_P2(p)=-1*P_v(2,2);
    
    disp(p);
    
end

X_P1 = reshape(X_P1,s,t);
Y_P1 = reshape(Y_P1,s,t);
Z_P1 = reshape(Z_P1,s,t);

X_P2 = reshape(X_P2,s,t);
Y_P2 = reshape(Y_P2,s,t);
Z_P2 = reshape(Z_P2,s,t);

Curvature_P1 = reshape(Curvature_P1,s,t);
Curvature_P2 = reshape(Curvature_P2,s,t);
Curvature_P1(isnan(Curvature_P1))=0;
Curvature_P2(isnan(Curvature_P2))=0;

Downsample_Factor=20;

X_resam=X(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Y_resam=Y(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Z_resam=Z(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);

X_P1_resam=X_P1(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Y_P1_resam=Y_P1(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Z_P1_resam=Z_P1(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);


X_P2_resam=X_P2(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Y_P2_resam=Y_P2(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Z_P2_resam=Z_P2(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);

test=X_P1_resam.*X_P2_resam+Y_P1_resam.*Y_P2_resam+Z_P1_resam.*Z_P2_resam;

quiver3(X_resam,Y_resam,Z_resam,X_P1_resam,Y_P1_resam,Z_P1_resam);
AA=500;
h=quiver3(X_resam,Y_resam,Z_resam,X_P2_resam,Y_P2_resam,Z_P2_resam,'linewidth',2,'color',[0 0 0]);
%adjust_quiver_arrowhead_size(h,2);
surface(X,Y,Z,'EdgeColor','none');
shading interp




%R1=1./real(P1)/1000;
%R2=1./real(P2)/1000;
%RK=1./real(K)/1000/1000;
%RH=1./real(H)/1000;



%imagesc(-R1,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(R1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(R1,1)-1));
%xlabel('(Micron)');
%ylabel('(Micron)');
%colorbar
%axis equal
%axis([0 Lateral_Spacing*(size(R1,2)-1) 0 Lateral_Spacing*(size(R1,1)-1)])
%caxis([5 15])
%caxis([0 20])


%% Next I need the radial direction (pre-experimental knowledge)
% Radial direction: in fact I want to find (polar coordinate) the tangent
% vector along R
NNN=2000;
R_consider=Lateral_Spacing*NN;
Center=[size(X,1)/2 size(X,2)/2]*Lateral_Spacing;


R=repmat([R_consider/NNN:R_consider/NNN:R_consider],NNN,1);
THETA=repmat([2*pi/NNN:2*pi/NNN:2*pi]',1,NNN);

[Rr Rtheta]=gradient(R);
[THETAr THETAtheta]=gradient(THETA);

X_wish=R.*cos(THETA)+Center(2);
Y_wish=R.*sin(THETA)+Center(1);

X_Polar=interp2(X,Y,X,X_wish,Y_wish);
Y_Polar=interp2(X,Y,Y,X_wish,Y_wish);
Z_Polar=interp2(X,Y,Z,X_wish,Y_wish);

[Xr Xtheta]=gradient(X_Polar);
[Yr Ytheta]=gradient(Y_Polar);
[Zr Ztheta]=gradient(Z_Polar);

% 現在來求
% back to XY coordinate

Center_2=[size(X,1)/2 size(X,2)/2]*Lateral_Spacing;


R_XY=((X-Center_2(2)).^2+(Y-Center_2(1)).^2).^0.5;
THETA_XY=atan2((Y-Center_2(1)),(X-Center_2(2)))+pi;

Xr_XY=interp2(R,THETA,Xr,R_XY,THETA_XY);
Yr_XY=interp2(R,THETA,Yr,R_XY,THETA_XY);
Zr_XY=interp2(R,THETA,Zr,R_XY,THETA_XY);

NORMr_XY=(Xr_XY.^2+Yr_XY.^2+Zr_XY.^2).^0.5;
Xr_XY=Xr_XY./NORMr_XY;
Yr_XY=Yr_XY./NORMr_XY;
Zr_XY=Zr_XY./NORMr_XY;      %這就是我要的radial方向的tangent vector

%現在來求tangent vector和P1 vector的夾角phi

PHI=acos(X_P1.*Xr_XY+Y_P1.*Yr_XY+Z_P1.*Zr_XY);

Curvature_Radial=(Curvature_P1.*cos(PHI).^2+Curvature_P2.*sin(PHI).^2);
Curvature_Radial(round(size(Curvature_Radial,1)/2)+1,:)=(Curvature_Radial(round(size(Curvature_Radial,1)/2)+2,:)+Curvature_Radial(round(size(Curvature_Radial,1)/2),:))/2;
Curvature_Radial(:,round(size(Curvature_Radial,1)/2)+1)=(Curvature_Radial(:,round(size(Curvature_Radial,1)/2)+2)+Curvature_Radial(:,round(size(Curvature_Radial,1)/2)))/2;
Curvature_Radial(isnan(Curvature_Radial))=0;



Curvature_Radial_TEST=(Curvature_P1.*sin(PHI).^2+Curvature_P2.*cos(PHI).^2);
Curvature_Radial_TEST(round(size(Curvature_Radial_TEST,1)/2)+1,:)=(Curvature_Radial_TEST(round(size(Curvature_Radial_TEST,1)/2)+2,:)+Curvature_Radial_TEST(round(size(Curvature_Radial_TEST,1)/2),:))/2;
Curvature_Radial_TEST(:,round(size(Curvature_Radial_TEST,1)/2)+1)=(Curvature_Radial_TEST(:,round(size(Curvature_Radial_TEST,1)/2)+2)+Curvature_Radial_TEST(:,round(size(Curvature_Radial_TEST,1)/2)))/2;
Curvature_Radial_TEST(isnan(Curvature_Radial_TEST))=0;

%imagesc(Curvature_Radial);
%Curvature_center=Curvature_Radial(round(size(Curvature_Radial,1)/2),round(size(Curvature_Radial,2)/2));
%caxis([Curvature_center/2 Curvature_center*1.5]);
%axis equal

n_cornea=1.3375;
Diopter_Radial=(n_cornea-1).*1E6.*Curvature_Radial;
Diopter_Radial_TEST=(n_cornea-1).*1E6.*Curvature_Radial_TEST;

Diopter_Principal_1=(n_cornea-1).*1E6.*Curvature_P1;
Diopter_Principal_2=(n_cornea-1).*1E6.*Curvature_P2;

subplot(2,1,1)
imagesc(Diopter_Radial,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial,1)-1));
caxis([35.23 43.61]);
colorbar 
axis equal

subplot(2,1,2)
imagesc(Diopter_Radial_TEST,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial_TEST,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial_TEST,1)-1));
caxis([35.23 43.61]);
colorbar 
axis equal

%imagesc(Diopter_Principal_2);
%Curvature_center=Diopter_Radial(round(size(Curvature_Radial,1)/2),round(size(Curvature_Radial,2)/2));
%caxis([35 45]);
%axis equal


%X_Radial=X-X_Center;
%Y_Radial=Y-Y_Center;
%Norm_Radial=(X_Radial.^2+Y_Radial.^2).^0.5;
%X_Radial=X_Radial./Norm_Radial;
%Y_Radial=Y_Radial./Norm_Radial;
%Z?

%X_Radial_resam=X_Radial(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
%Y_Radial_resam=Y_Radial(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);

%quiver3(X_resam,Y_resam,Z_resam,X_P2_resam,Y_P2_resam,Z_P2_resam);


%% trans to polar