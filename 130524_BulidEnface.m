clear all

% larger image index = larger OPD => deeper
%Too bad axial resolution
fitting=0;
new_bit_depth=16;

Sample_Path='D:\Users\TuanShu\130531_MK\COA\';
Sample_Path_Bscan='D:\Users\TuanShu\130531_MK\COA\Enface';

set_peak_index=250;
Averaging_Factor=1;
Sampling_Rate=30;
PZT_Speed=0.5;         %micron/sec

TD_OCT_point_Spacing=PZT_Speed/(Sampling_Rate/Averaging_Factor);   
Binning_Factor=1;

%Lateral_Spacing=Pixel_Size/Porjection_Lens_Focal_Length*Objective_Focal_Length;
Lateral_Spacing=2.8846*Binning_Factor/10*4;

image_index=1009:2610;

ROI=[1 768;1 1024];      %up, down, left, right

Low_Pass=150;   %pixel
High_Pass=450;

%ROI=[43 342;363 662];      %up, down, left, right
%ROI=[1 384;257 768];      %up, down, left, right
%%ROI=[151 234;457 568];      %up, down, left, right
Image_Stack(1:ceil(length(image_index)/Averaging_Factor),1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;

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
            Image_Temp_2=Image_Temp(ROI(1,1):ROI(1,2)/Binning_Factor,ROI(2,1):ROI(2,2)/Binning_Factor)+Image_Temp_2;
            Image_Stack(ceil(p/Averaging_Factor),:,:)=Image_Temp_2;
            Image_Temp_2=0;
        else
            Image_Temp_2=Image_Temp(ROI(1,1):ROI(1,2)/Binning_Factor,ROI(2,1):ROI(2,2)/Binning_Factor)+Image_Temp_2;
        end
        %Image_Stack(ceil(p/Averaging_Factor),:,:)=Image_Temp;
        %Image=imread(sprintf('1.jpg',image_index(p)));    
        disp(p);
    end
    Image_Stack(size(Image_Stack,1),:,:)=Image_Stack(size(Image_Stack,1)-1,:,:);
    Position_micron=(1:size(Image_Stack,1))*TD_OCT_point_Spacing;

    clear DC_volume
    [Max_value Max_index]=max(Image_Stack(:,250,250));
    %Image_Stack=circshift(Image_Stack,set_peak_index-Max_index);

    plot(real(Image_Stack(:,250,250)));
    
    FFT_Stack=fft(Image_Stack,[],1);
    plot(real(FFT_Stack(:,473,286)));

    FFT_Stack(1:Low_Pass,:,:)=0;
    FFT_Stack(round(size(FFT_Stack,1)/2):end,:,:)=0;
    FFT_Stack(High_Pass:end,:,:)=0;

    clear Image_Stack
    Image_Stack_New=(ifft(FFT_Stack,[],1));
    
    clear FFT_Stack
    %Max_Value_of_The_Volume=max(abs(Image_Stack_New(:)));
    %Image_Stack_New_Normalized=Image_Stack_New./max(abs(Image_Stack_New(:)));
    %Image_Stack_New_Env=log10(abs(Image_Stack_New)/max(abs(Image_Stack_New(:)))*10);
    Image_Stack_New_Env=(abs(Image_Stack_New)/max(abs(Image_Stack_New(:))));
    plot(real(Image_Stack_New_Env(:,250,250)));

    [maxvalue maxindex]=max(Image_Stack_New_Env);
    Height_Map(:,:)=(maxindex-maxindex(1,1,1))*TD_OCT_point_Spacing;
    %%
    P1=[40 990];
    
    P2=[20 50];
    
    P3=[750 990];
    
    Z1=Height_Map(P1(1),P1(2))-Height_Map(P2(1),P2(2));
    Z2=Height_Map(P1(1),P1(2))-Height_Map(P3(1),P3(2));
    
    X1=P1(1)-P2(1);
    X2=P1(1)-P3(1);
    
    Y1=P1(2)-P2(2);
    Y2=P1(2)-P3(2);
    
    a=(Z1*Y2-Y1*Z2)/(X1*Y2-Y1*X2);
    
    b=(X1*Z2-Z1*X2)/(X1*Y2-Y1*X2);
    
    X=repmat([1:size(Height_Map,1)]',1,size(Height_Map,2));
    
    Y=repmat(1:size(Height_Map,2),size(Height_Map,1),1);
       imagesc(Height_Map,'xdata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(2,2),'ydata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(1,2));
    xlabel('(micron)');
    ylabel('(micron)');
    axis equal
    %Height_Map_1=Height_Map-(a.*X+b.*Y);

    
    
 
    
    imagesc(Height_Map,'xdata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(2,2),'ydata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(1,2));
    xlabel('(micron)');
    ylabel('(micron)');
%%

    
    %Image_Stack_New_Env=Image_Stack_New_Env/max(abs(Image_Stack_New_Env(:)));
        cd(Sample_Path_Bscan);
        %dlmwrite('Max_Value_Sam.txt',Max_Value_of_The_Volume,'delimiter','\t','newline','pc');
        for rr=1:size(Image_Stack_New_Env,1)
            Temp_Bscan(:,:)=(Image_Stack_New_Env(rr,:,:));
            imwrite(Temp_Bscan,sprintf('Enface_Sam_AVE%i_BIN%d_%i.bmp',Averaging_Factor,Binning_Factor,rr));%,'Bitdepth',new_bit_depth);
            %dlmwrite(sprintf('Bscan_Sam_%i.txt',rr),Temp_Bscan,'delimiter','\t','newline','pc');
            disp(rr);
        end
imagesc(Height_Map,'xdata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(2,2),'ydata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(1,2));
xlabel('(micron)');
ylabel('(micron)');
axis equal

dlmwrite('Height_Map.txt',Height_Map,'delimiter','\t','newline','pc');
imwrite(Height_Map,sprintf('Height_Map.bmp',rr));%,'Bitdepth',new_bit_depth);


%%
%plot(Position_micron,Image_Stack_New_Env(:,120,160));
%Frame(:,:)=abs(Image_Stack_New(433,:,:));
%imagesc(Frame,'xdata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(2,2),'ydata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(1,2));
%xlabel('(micron)');
%ylabel('(micron)');
%axis equal
