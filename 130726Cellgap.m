clear all; close all; clc

% larger image index = larger OPD => deeper
%Too bad axial resolution
fitting=0;
new_bit_depth=16;

Sample_Path='D:\Users\MengKo\130726_Cellgap335-18';
Sample_Path_Save='D:\Users\MengKo\130829\Cellgap';

set_peak_index=250;
Averaging_Factor=1;
Sampling_Rate=30*4;
PZT_Speed=0.5*4;         %micron/sec $ 1/4partial

TD_OCT_point_Spacing=PZT_Speed/(Sampling_Rate/Averaging_Factor);   
Binning_Factor=1;

%Lateral_Spacing=Pixel_Size/Porjection_Lens_Focal_Length*Objective_Focal_Length;
Lateral_Spacing=2.8846*Binning_Factor/10*4;

image_index=0:4257;

ROI=[1 136;1 1024];      %up, down, left, right

Low_Pass=200;   %pixel
High_Pass=400;

%ROI=[43 342;361 662];      %up, down, left, right
%ROI=[1 384;257 763];      %up, down, left, right
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

    plot(real(Image_Stack(:,50,510)));
    clear DC_volume
    [Max_value Max_index]=max(Image_Stack(:,50,510));
    Image_Stack=circshift(Image_Stack,set_peak_index-Max_index);

    FFT_Stack=fft(Image_Stack,[],1);
   figure; plot(real(FFT_Stack(:,50,510)));

    FFT_Stack(1:Low_Pass,:,:)=0;
    FFT_Stack(round(size(FFT_Stack,1)/2):end,:,:)=0;
    FFT_Stack(High_Pass:end,:,:)=0;

    Image_Stack_New=(ifft(FFT_Stack,[],1));
    %Max_Value_of_The_Volume=max(abs(Image_Stack_New(:)));
    %Image_Stack_New_Normalized=Image_Stack_New./max(abs(Image_Stack_New(:)));
    %Image_Stack_New_Env=log10(abs(Image_Stack_New)/max(abs(Image_Stack_New(:)))*10);
    clear Image_Stack
    a=Image_Stack_New;
    b=a(:,103,576);
    c=abs(b);
    figure;plot(c)
    clear a b c
    Image_Stack_New_Env=(abs(Image_Stack_New)/max(abs(Image_Stack_New(:))));
    a=Image_Stack_New_Env(3900:end,:,:);
    b=vertcat(a,Image_Stack_New_Env);
    plot(b(:,103,576))
 %%   
    pixel_low=1;
    pixel_high=1300;
    cd (Sample_Path_Bscan)
    for tt=1:136
        Txt_Bscan=Image_Stack_New_Env(pixel_low:pixel_high,tt,:);
        dlmwrite(sprintf('Bscan_%i.txt',tt),Txt_Bscan,'delimiter','\t','newline','pc');
        display(tt)
    end    
    
    %%
    pixel_low=1;
    pixel_high=700;
%     tempp=Image_Stack_New_Env(pixel_low:pixel_high,:,:);
    tempp=b(pixel_low:pixel_high,:,:);
    [ma mb]=max(tempp(100:400,:,:));
    [mc md]=max(tempp(400:700,:,:));
%     [me mf]=max(tempp(600:800,:,:));
%     [mg mh]=max(tempp(800:900,:,:));
    toto=(md+400)-(mb+100);
    t1=toto*TD_OCT_point_Spacing;
    cellgap1=permute(t1,[2 3 1]);
%     toto2=(mf+600)-(md+400);
%     t2=toto2*TD_OCT_point_Spacing;
%     cellgap2=permute(t2,[2 3 1]);
%     toto3=(mh+800)-(mf+600);
%     t3=toto3*TD_OCT_point_Spacing;
%     cellgap3=permute(t3,[2 3 1]);
%     tmp1=zeros(700,1020);
%     tmp2=zeros(700,1020);
%     for pp=1:136
% %     tmp=dlmread(sprintf('Bscan_%i.txt',pp)); %cellgap
%     tmp=tempp()
%     tmp1(200:500,:)=tmp(200:500,:);
%     [ma mb]=max(tmp1);
%     tmp2(600:800,:)=tmp(600:800,:);
%     [mc md]=max(tmp2);
%     cellgap(pp,:)=(md-mb)*TD_OCT_point_Spacing;
%     display(pp)
%     end
figure; imagesc(cellgap1)
 cd(Sample_Path_Save)
 dlmwrite('Cellgap1_OPD-330-13.txt',cellgap1,'delimiter','\t','newline','pc');
%  dlmwrite('Cellgap2_OPD-912-2.txt',cellgap2,'delimiter','\t','newline','pc');
%  dlmwrite('Cellgap3_OPD-912-2.txt',cellgap3,'delimiter','\t','newline','pc');
 
 %%
    
%     [maxvalue maxindex]=max(Image_Stack_New_Env);
%     Height_Map(:,:)=(maxindex-maxindex(1,1,1))*TD_OCT_point_Spacing;
%     %%
%     P1=[40 990];
%     
%     P2=[20 50];
%     
%     P3=[750 990];
%     
%     Z1=Height_Map(P1(1),P1(2))-Height_Map(P2(1),P2(2));
%     Z2=Height_Map(P1(1),P1(2))-Height_Map(P3(1),P3(2));
%     
%     X1=P1(1)-P2(1);
%     X2=P1(1)-P3(1);
%     
%     Y1=P1(2)-P2(2);
%     Y2=P1(2)-P3(2);
%     
%     a=(Z1*Y2-Y1*Z2)/(X1*Y2-Y1*X2);
%     
%     b=(X1*Z2-Z1*X2)/(X1*Y2-Y1*X2);
%     
%     X=repmat([1:size(Height_Map,1)]',1,size(Height_Map,2));
%     
%     Y=repmat(1:size(Height_Map,2),size(Height_Map,1),1);
%        imagesc(Height_Map,'xdata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(2,2),'ydata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(1,2));
%     xlabel('(micron)');
%     ylabel('(micron)');
%     axis equal
%     Height_Map_1=Height_Map-(a.*X+b.*Y);

    
    
 
    
%     imagesc(Height_Map_1,'xdata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(2,2),'ydata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(1,2));
%     xlabel('(micron)');
%     ylabel('(micron)');
%%

    
    %Image_Stack_New_Env=Image_Stack_New_Env/max(abs(Image_Stack_New_Env(:)));
        cd(Sample_Path_Bscan);
        %dlmwrite('Max_Value_Sam.txt',Max_Value_of_The_Volume,'delimiter','\t','newline','pc');
        for rr=1:size(Image_Stack_New_Env,1)
            Temp_Bscan(:,:)=(Image_Stack_New_Env(rr,:,:));
            imwrite(Temp_Bscan,sprintf('Enface_Sam_%i.bmp',rr));%,'Bitdepth',new_bit_depth);
            %dlmwrite(sprintf('Bscan_Sam_%i.txt',rr),Temp_Bscan,'delimiter','\t','newline','pc');
            disp(rr);
        end
% imagesc(Height_Map,'xdata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(2,2),'ydata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(1,2));
% xlabel('(micron)');
% ylabel('(micron)');
% axis equal
% 
% dlmwrite('Height_Map.txt',Height_Map,'delimiter','\t','newline','pc');
% imwrite(Height_Map,sprintf('Height_Map.bmp',rr));%,'Bitdepth',new_bit_depth);


%%
%plot(Position_micron,Image_Stack_New_Env(:,120,160));
%Frame(:,:)=abs(Image_Stack_New(433,:,:));
%imagesc(Frame,'xdata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(2,2),'ydata',Lateral_Spacing:Lateral_Spacing:Lateral_Spacing*ROI(1,2));
%xlabel('(micron)');
%ylabel('(micron)');
%axis equal
