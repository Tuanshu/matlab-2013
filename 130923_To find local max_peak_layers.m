clear all


Path='D:\Users\TuanShu\130913_Skin';

cd(sprintf('%s\\',Path));

Threshold=20;


Data=imread('line.png');   

Data=Data(:,:,1);   
[maxvalue maxindex]=max(Data(:,512));
plot((1:768)-maxindex,(Data(:,512)-min(Data(:,512))));
xlabel('Pixel Number');
ylabel('Signal (a.u.)');
plot(1:1024-max(Data(:,512)),Data(:,512))

Separation_Index_1=470;
Separation_Index_2=520;


%%diff是後面減前面, check

D_Signal=diff(Signal,1);
D_Signal(size(D_Signal,1)+1,:)=D_Signal(size(D_Signal,1),:);
plot(Position,D_Signal);

q=1;

plot(Position,Signal(:,1),Position,D_Signal(:,1));
for p=1:size(D_Signal,2)    
    D_Signal_Temp=D_Signal(:,p);
    q=1;
    Start_index=1;
    while Start_index<length(D_Signal_Temp)
        Start_index=find(D_Signal_Temp(Start_index:end)>0,1,'first')+Start_index-1;
        if find(D_Signal_Temp(Start_index:end)<0,1,'first') >0
            stored_max(q,p)=find(D_Signal_Temp(Start_index:end)<0,1,'first')+Start_index-1;
            Start_index=find(D_Signal_Temp(Start_index:end)<0,1,'first')+Start_index-1;
            q=q+1;
        else
            break;
        end
    end
end
stored_max(stored_max==0)=1;
stored_max_Position=Position(stored_max);

Thickness_Array=diff(stored_max_Position);
Thickness_Array(Thickness_Array<0)=0;
Thickness_Array(size(Thickness_Array,1)+1,:)=0;
for p=1:size(stored_max,2)
    linetemp=stored_max(:,p);
    Q=linetemp(linetemp~=1);
    Number_of_Layers(p)=length(Q)-1;
    Total_Thickness(p)=max(stored_max_Position(stored_max(:,p)~=1,p))-min(stored_max_Position(stored_max(:,p)~=1,p));
    Ave_Thickness(p)=Total_Thickness(p)/Number_of_Layers(p);
    Min_Thickness(p)=min(Thickness_Array(Thickness_Array(:,p)~=0,p));
    Max_Thickness(p)=max(Thickness_Array(Thickness_Array(:,p)~=0,p));
    SD_Thickness(p)=std(Thickness_Array(Thickness_Array(:,p)~=0,p));
end
%%

N=50;

plot(Position,Signal(:,N),'-',Position(stored_max(:,N)),Signal(stored_max(:,N),N),'o');
xlabel('Thickness (micron)');
ylabel('Interference Signal');

fprintf('\n\n\n\n\nN:                 %d \n',N);
fprintf('Number_of_Layers:  %d \n',Number_of_Layers(N));
fprintf('Total_Thickness:   %d micron\n',Total_Thickness(N));
fprintf('Ave_Thickness:     %d micron\n',Ave_Thickness(N));
fprintf('Min_Thickness:     %d micron\n',Min_Thickness(N));
fprintf('Max_Thickness:     %d micron\n',Max_Thickness(N));
fprintf('SD_Thickness:      %d micron\n',SD_Thickness(N));

Output=[[1:N]' Number_of_Layers' Total_Thickness' Ave_Thickness' Min_Thickness' Max_Thickness' SD_Thickness'];
%Number of Layers
%Ave Thickness
%Min Thickness
%Max Thikcness
%SD of Thickness