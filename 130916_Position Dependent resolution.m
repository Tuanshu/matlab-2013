clear all
Axial_resolution=10;

r=0:0.1:5500;    %micron
R=7800;      %micron

Cornea_function=(R^2-r.^2).^0.5;

Height=max(Cornea_function)-min(Cornea_function);

for j=1:length(r)
    if isempty(find(Cornea_function<(Cornea_function(j)-Axial_resolution),1,'first'))
        Lateral_resolution(j)=Lateral_resolution(j-1);
    else
        Lateral_resolution(j)=r(find(Cornea_function<(Cornea_function(j)-Axial_resolution),1,'first'))-r(j);
    end
end
Cornea_function_2=R-Cornea_function;
Cornea_function_2_pre=0;
p=1;
while(isempty(find(Cornea_function_2>(Cornea_function_2_pre+36),1,'first'))==0)
    index_array(p)=find(Cornea_function_2>(Cornea_function_2_pre+36),1,'first');
    Cornea_function_2_pre=Cornea_function_2(index_array(p));
    p=p+1;
end

r_array=r(index_array);
plot(r_array, '.');
xlabel('Number of Frames (#)');
ylabel('Radius of Ring (micron)');

plot(r,Lateral_resolution);
xlabel('Radial Position (micron)');
ylabel('Fringe Spacing (micron)');