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
Lowpass=20;
Longpass=350;
 for p=1:length(image_index)
        Image_Temp=imread(sprintf('%i.png',image_index(p)));   
        Image_Temp=Image_Temp(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2));
        FFT_X=fft(Image_Temp,[],2);
        FFT_X(:,1:Lowpass)=0;
        FFT_X(:,Longpass:end)=0;    
        New_X=ifft(FFT_X,[],2);
        FFT_Y=fft(Image_Temp,[],1);
        FFT_Y(1:Lowpass,:)=0;
        FFT_Y(Longpass:end,:)=0;    
        New_Y=ifft(FFT_Y,[],1);
        New_XY(p,:,:)=abs(abs(New_X).^2+abs(New_Y).^2).^0.5;
 end
%%
plot(1:768,FFT_X(332,:))

%%
image_enface(:,:)=New_XY(146,:,:);
imagesc(image_enface);

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
