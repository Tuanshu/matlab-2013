clear all

cd('D:\Users\TuanShu\131015\');

Data1_2D=importdata('_122');
%OPD4=47.769897

Data1_2D=rot90(Data1_2D);

index=60;
Data2_2D=0;
for p=1:length(index)
    Data2_2D_Temp=importdata(sprintf('_%i',index(p)));
    Data2_2D=Data2_2D+Data2_2D_Temp;
    disp(p);
end
Data2_2D=Data2_2D/length(index);
    
%OPD1=61.605919

Data2_2D=rot90(Data2_2D);

cut=200;

Data1_2D=Data1_2D(1:(length(Data1_2D)-cut+1),:);
Data2_2D=Data2_2D(1:(length(Data1_2D)-cut+1),:);
Data_Diff=Data2_2D-Data1_2D;
Data_Diff(Data_Diff<1)=1;
imagesc(log10(Data_Diff));
colormap(gray);
xlabel('X Pixel Number');
ylabel('Y Pixel Number');

imagesc(log10(Data_Diff));
colormap(gray);
xlabel('X Pixel Number');
ylabel('Y Pixel Number');


imagesc((Data_Diff));
caxis([0 250]);
colormap(gray);
xlabel('X Pixel Number');
ylabel('Y Pixel Number');