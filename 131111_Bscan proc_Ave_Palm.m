clear all

cd('D:\Users\TuanShu\131015\');


index_Ref=100:130;
Data1_2D=0;
for p=1:length(index_Ref)
    Data1_2D_Temp=importdata(sprintf('Palm_%i',index_Ref(p)));
    Data1_2D=Data1_2D+Data1_2D_Temp;
    disp(p);
end
Data1_2D=Data1_2D/length(index_Ref);


Data1_2D=rot90(Data1_2D);


imagesc((Data1_2D));
caxis([0 250]);
colormap(gray);
xlabel('X Pixel Number');
ylabel('Y Pixel Number');

%%
index=1:90;
Data2_2D=0;
ROI=[605 655 360 440];
TH=0;   %50
pp=0;
for p=1:length(index)
    Data2_2D_Temp=importdata(sprintf('Palm_%i',index(p)));
    if mean(mean(Data2_2D_Temp(ROI(1):ROI(2),ROI(3):ROI(4))))>TH
        Data2_2D=Data2_2D+Data2_2D_Temp;
        pp=pp+1;        
        disp(pp);
    end
end
Data2_2D=Data2_2D/pp;
    
%OPD1=61.605919

Data2_2D=rot90(Data2_2D);


imagesc((Data2_2D));
caxis([0 250]);
colormap(gray);
xlabel('X Pixel Number');
ylabel('Y Pixel Number');

%%

cut=100;

Data1_2D_Cut=Data1_2D(1:(length(Data1_2D)-cut+1),:);
Data2_2D_Cut=Data2_2D(1:(length(Data2_2D)-cut+1),:);
Data_Diff=Data2_2D_Cut-Data1_2D_Cut;
Data_Diff(Data_Diff<1)=1;
Data_Diff=Data_Diff(size(Data_Diff,1):-1:1,:);
imagesc(log10(Data_Diff));
colormap(gray);
xlabel('X Pixel Number');
ylabel('Y Pixel Number');

imagesc(log10(Data_Diff));
colormap(gray);
xlabel('X Pixel Number');
ylabel('Y Pixel Number');


imagesc((Data_Diff),'xdata',4.65*[1:size(Data_Diff,2)],'ydata',0.4*[1:size(Data_Diff,1)]);
set(gca,'YDir','normal')
xlim([1000 4000]);
ylim([0 300]);
caxis([0 100]);
colormap(gray);
xlabel('X (micron)');
ylabel('Y (micron)');




imagesc(10*log10(Data_Diff),'xdata',4.65*[1:size(Data_Diff,2)],'ydata',0.4*[1:size(Data_Diff,1)]);
set(gca,'YDir','normal')
xlim([1000 4000]);
ylim([0 300]);
caxis([12 20]);
colormap(gray);
xlabel('X (micron)');
ylabel('Y (micron)');