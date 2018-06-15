clear all
filter=[-1 0 1;0 1 2; 1 2 3];
for p=3:3
image=imread(sprintf('D:\\Users\\TuanShu\\130820\\Ring %i.png',p));%,'Bitdepth',new_bit_depth);
image=abs(image-129);
end

imagesc(image);
colormap(gray);
axis equal