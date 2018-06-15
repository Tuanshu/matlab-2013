clear all

cd('D:\Users\TuanShu\131014\');

Data1_2D=importdata('30');
%OPD4=47.769897

Data1_2D=rot90(Data1_2D);

Data2_2D=importdata('1568');
%OPD1=61.605919

Data2_2D=rot90(Data2_2D);

cut=100;

Data1_2D=Data1_2D(cut:end,:);
Data2_2D=Data2_2D(cut:end,:);
Data_Diff=Data2_2D-Data1_2D;
Data_Diff(Data_Diff<1)=1;
imagesc(log10(Data_Diff));
colormap(gray);
caxis([0 250]);
xlabel('X Pixel Number');
ylabel('Y Pixel Number');