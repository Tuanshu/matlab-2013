clear all
Axial_resolution=10;
R_Ref=7800;
R=5000:100:13000;      %micron
MAXMAX=6;

r=0:5500;    %micron




for p=1:length(R)
    r=0:5500;    %micron

    Cornea_function=abs((R(p)-(R(p)^2-r.^2).^0.5)-(R_Ref-(R_Ref^2-r.^2).^0.5));
    index=find(Cornea_function>6,1,'first');
    if isempty(index)
        FOV(p)=5500;
    else
        FOV(p)=r(find((Cornea_function)>6,1,'first'));
    end
end
plot(R/1000,2*FOV/1000);
xlabel('Averaged Radius of Curvature of Cornea (mm)');
ylabel('FOV (mm)');
