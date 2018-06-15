clear all

level=0:10;

total_cost=0;
total_percent=0;
cost_to_overcome_the_level(1:length(level))=0;
cost2_to_overcome_the_level(1:length(level))=0;
for p=1:length(level)
        if level(p)==0
            cost=64953;
            cost2=10;
            perpcent=1;
        elseif level(p) ==1
            cost=64953;
             cost2=20;
            perpcent=1;
        elseif level(p) ==2
            cost=64953;
            cost2=30;
            perpcent=1;
        elseif level(p) ==3
            cost=64953;
            cost2=40;

            perpcent=1;        
        elseif level(p) ==4
            cost=77943;
            cost2=50;
            perpcent=0.9;       
            penalty=1;
        elseif level(p) ==5
            cost=90934;
            cost2=60;
            perpcent=0.8;               
            penalty=1;
        elseif level(p) ==6
            cost=103924;
            cost2=70;
            perpcent=0.6;                       
            penalty=1;
        elseif level(p) ==7
            cost=116915;
            cost2=80;
            perpcent=0.8;                       
            penalty=3;
        elseif level(p) ==8
            cost=129906;
            cost2=90;
            perpcent=0.9;                 
            penalty=3;
        elseif level(p) ==9
            cost=129906;
            cost2=100;
            perpcent=0.5;                               
            penalty=3;
        end
        if perpcent==1
            cost_to_overcome_the_level(p)=cost;
            cost2_to_overcome_the_level(p)=cost2;

        elseif penalty==1
            cost_to_overcome_the_level(p)=cost+(1-perpcent)/perpcent*(cost+cost_to_overcome_the_level(p-1));
            cost2_to_overcome_the_level(p)=cost2+(1-perpcent)/perpcent*(cost2+cost2_to_overcome_the_level(p-1));

        elseif penalty==3
            cost_to_overcome_the_level(p)=cost+(1-perpcent)/perpcent*(cost+cost_to_overcome_the_level(p-1)+cost_to_overcome_the_level(p-2)+cost_to_overcome_the_level(p-3));
            cost2_to_overcome_the_level(p)=cost2+(1-perpcent)/perpcent*(cost2+cost2_to_overcome_the_level(p-1)+cost2_to_overcome_the_level(p-2)+cost2_to_overcome_the_level(p-3));

        end
        total_cost(p)=sum(cost_to_overcome_the_level(1:p));
        total_cost2(p)=sum(cost2_to_overcome_the_level(1:p));

end
total_cost=total_cost';
total_cost2=total_cost2';
cost_to_overcome_the_levelcost_to_overcome_the_level=cost_to_overcome_the_level';
cost2_to_overcome_the_levelcost_to_overcome_the_level=cost2_to_overcome_the_level';