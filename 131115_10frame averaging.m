clear all


Frame_index=1:530;

Starting_N=[11 11 10 17 13 12 15 19 17 8];

for p=1:length(Frame_index)
    Image_Sum(1:488,1:648)=0;
    for q=1:length(Starting_N)
        cd(sprintf('D:\\Users\\TuanShu\\131114\\%i\\',q));
        Image_Temp=double(imread(sprintf('%08i.png',Starting_N(q)+Frame_index(p)-1)));
        Image_Sum=Image_Sum+Image_Temp(:,:,1);
    end
    %Image_Sum=(Image_Sum/length(Starting_N));
    
    cd('D:\Users\TuanShu\131114\sum\');
    imwrite(uint16(Image_Sum),sprintf('%08i.png',Frame_index(p)),'png','bitdepth',16);
    disp(p);
end
