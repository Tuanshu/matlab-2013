clear all

Frame_index=1:240;
unit_length=0.19;
Number_of_Layers_Min=5;

cd('D:\Users\TuanShu\131126_skin');
Volume(1:length(Frame_index),1:488,1:648)=0;
for p=1:length(Frame_index)
    Image_Temp=imread(sprintf('%08i.tif',Frame_index(p)));
    Volume(p,:,:)=Image_Temp(:,:,1);
    disp(p);
end

imagesc(Volume(:,:,300));

Separation_Index_1=470;
Separation_Index_2=520;


for w=1:size(Volume,3)
    D_Frame=diff(Volume(:,:,w),1);
    D_Frame(size(D_Frame,1)+1,:)=D_Frame(size(D_Frame,1),:);
    for p=1:size(D_Frame,2)    
        D_Signal=D_Frame(:,p);
        q=1;
        Start_index=1;
        while Start_index<length(D_Signal)
            Start_index=find(D_Signal(Start_index:end)>0,1,'first')+Start_index-1;
            if find(D_Signal(Start_index:end)<0,1,'first') >0
                stored_max(q,p,w)=find(D_Signal(Start_index:end)<0,1,'first')+Start_index-1;
                Start_index=find(D_Signal(Start_index:end)<0,1,'first')+Start_index-1;
                q=q+1;
            else
                break;
            end
        end
        
    end
    disp(w);
end
stored_max(stored_max==0)=1;

%%
Threshold=30;
stored_max_threshold=stored_max;

for w=1:size(stored_max_threshold,3)
    for p=1:size(stored_max_threshold,2)
        for q=1:size(stored_max_threshold,1)
            if Volume(stored_max_threshold(q,p,w),p,w)<Threshold
                stored_max_threshold(q,p,w)=1;
            end
        end
    end
    disp(w);
end
%%
N=248;
imagesc(Volume(:,:,N));
hold on
plot(stored_max_threshold(:,:,N)','.');
hold off

%%

N=250;
Frame(:,:)=Volume(:,N,:);
imagesc(Frame,'ydata',[0:unit_length:size(Frame,1)*unit_length]);
hold on
stored_max_2d(:,:)=stored_max_threshold(:,N,:);
plot(stored_max_2d'*unit_length,'.');
hold off

%%
stored_max_Position=stored_max_threshold*unit_length;
stored_max_Position(stored_max_Position==unit_length)=0;
Thickness_Array=diff(stored_max_Position,[],1);
Thickness_Array(Thickness_Array<0)=0;
Thickness_Array(size(Thickness_Array,1)+1,:,:)=0;

for w=1:size(stored_max_threshold,3)
    for p=1:size(stored_max_threshold,2)
        linetemp=stored_max_threshold(:,p,w);
        max_temp=linetemp(linetemp~=1);
        max_position_temp=unit_length*max_temp;
        if length(max_temp)>Number_of_Layers_Min
            Number_of_Layers(p,w)=length(max_temp)-1;
            Total_Thickness(p,w)=max(max_position_temp)-min(max_position_temp);
            Ave_Thickness(p,w)=Total_Thickness(p,w)/Number_of_Layers(p,w);
            Min_Thickness(p,w)=min(Thickness_Array(:,p,w));
            Max_Thickness(p,w)=max(Thickness_Array(:,p,w));
            SD_Thickness(p,w)=std(Thickness_Array(:,p,w));
        else            
            Number_of_Layers(p,w)=-1;
            Total_Thickness(p,w)=-1;
            Ave_Thickness(p,w)=-1;
            Min_Thickness(p,w)=-1;
            Max_Thickness(p,w)=-1;
            SD_Thickness(p,w)=-1;
        end
    end
    disp(w);
end

%%
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