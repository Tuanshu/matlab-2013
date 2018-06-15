clear all

% larger image index = larger OPD => deeper
%Too bad axial resolution
fitting=0;
new_bit_depth=16;

Sample_Path='"D:\Users\TuanShu\121106_Lens2\';
Sample_Path_Cscan='D:\Users\TuanShu\121106_Lens2\Cscan';
if exist(Sample_Path_Cscan,'file')  ~= 7
    mkdir('D:\Users\TuanShu\121106_Lens2\Cscan');
end
image_index=200:450;
cd(Sample_Path);
ROI=[1 768;129 896];      %up, down, left, right
Lateral_Resolutiion=5;      %micron
Lowpass=10;
Longpass=150;
 for p=1:length(image_index)

        Center=[393 386];

        Image_Temp=imread(sprintf('%i.png',image_index(p)));   
        Image_Temp=Image_Temp(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2));
        
       
        [r c] = size(Image_Temp);
        [X Y] = meshgrid(1:c,1:r);
        
        
        NNN=360;
        
        r_wish=repmat([1:360/NNN:360],NNN,1);
        theta_wish=repmat([2*pi/NNN:2*pi/NNN:2*pi]',1,NNN);
        
        X_wish=r_wish.*cos(theta_wish)+Center(2);
        Y_wish=r_wish.*sin(theta_wish)+Center(1);

        image_Rtheta=interp2(X,Y,double(Image_Temp),X_wish,Y_wish);
        
        
        FFT_R=fft(image_Rtheta,[],2);
        FFT_R=FFT_R-mean(FFT_R(:,round(size(FFT_R,2)/2)));
        FFT_R(:,1:Lowpass)=0;
        FFT_R(:,Longpass:end)=0;    
        image_Rtheta_New=abs(ifft(FFT_R,[],2));
        image_Rtheta_New(:,320:end)=0;
       
        X_grid=repmat(1:768,768,1);
        Y_grid=repmat([1:768]',1,768);
        
        [r_rt c_rt] = size(Image_Temp);
        [X_rt Y_rt] = meshgrid(1:c_rt,1:r_rt);
        
        r_XY=((X-Center(2)).^2+(Y-Center(1)).^2).^0.5;
        theta_XY=atan2((Y-Center(1)),(X-Center(2)))+pi;
        
        image_XY_New=interp2(r_wish,theta_wish,image_Rtheta_New,r_XY,theta_XY);
        image_XY_New(isnan(image_XY_New))=0;
        image_volume(p,:,:)=image_XY_New;
        disp(p);
 end
%%
plot(real(FFT_R(132,:)));


%%
image_enface(:,:)=image_volume(1,:,:);
imagesc(image_enface);


%%
        cd(Sample_Path);

        p=383;        
        Center=[393 386];
        Lowpass=20;
        Longpass=250;
        Image_Temp=imread(sprintf('%i.png',(p)));   
        Image_Temp=Image_Temp(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2));
        
       
        [r c] = size(Image_Temp);
        [X Y] = meshgrid(1:c,1:r);
        
        NNN=360;
        
        r_wish=repmat([1:360/NNN:360],NNN,1);
        theta_wish=repmat([2*pi/NNN:2*pi/NNN:2*pi]',1,NNN);
        
        X_wish=r_wish.*cos(theta_wish)+Center(2);
        Y_wish=r_wish.*sin(theta_wish)+Center(1);

        image_Rtheta=interp2(X,Y,double(Image_Temp),X_wish,Y_wish);
        
        
        FFT_R=fft(image_Rtheta,[],2);
        FFT_R=FFT_R-mean(FFT_R(:,round(size(FFT_R,2)/2)));
        FFT_R(:,1:Lowpass)=0;
        FFT_R(:,Longpass:end)=0;    
        image_Rtheta_New=abs(ifft(FFT_R,[],2));
        
        image_Rtheta_New(:,1:10)=0;
        
                imagesc(image_Rtheta_New);

        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca,'XColor','white');
        set(gca,'YColor','white');
        
        image_Rtheta_New(:,320:end)=0;

       
        X_grid=repmat(1:768,768,1);
        Y_grid=repmat([1:768]',1,768);
        
        [r_rt c_rt] = size(Image_Temp);
        [X_rt Y_rt] = meshgrid(1:c_rt,1:r_rt);
        
        r_XY=((X-Center(2)).^2+(Y-Center(1)).^2).^0.5;
        theta_XY=atan2((Y-Center(1)),(X-Center(2)))+pi;
        
        image_XY_New=interp2(r_wish,theta_wish,image_Rtheta_New,r_XY,theta_XY);
        
        
        imagesc(image_XY_New);
        axis equal
        caxis([3 12]);
        colormap(gray)
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca,'XColor','white');
        set(gca,'YColor','white');

        %imagesc(Image_Temp);
        %axis equal
        %caxis([3 12]);
        %colormap(gray)
        %set(gca, 'XTick', []);
        %set(gca, 'YTick', []);
        %set(gca,'XColor','white');
        %set(gca,'YColor','white');

%%

    image_volume=image_volume/max(max(max(image_volume)));
    %Image_Stack_New_Env=Image_Stack_New_Env/max(abs(Image_Stack_New_Env(:)));
        cd(Sample_Path_Cscan);
        %dlmwrite('Max_Value_Sam.txt',Max_Value_of_The_Volume,'delimiter','\t','newline','pc');
        for rr=1:size(image_volume,1)
            Temp_Bscan(:,:)=(image_volume(rr,:,:));
            imwrite(Temp_Bscan,sprintf('Enface_%i.bmp',rr));%,'Bitdepth',new_bit_depth);
            %dlmwrite(sprintf('Bscan_Sam_%i.txt',rr),Temp_Bscan,'delimiter','\t','newline','pc');
            disp(rr);
        end