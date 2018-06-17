function [net,gene]=INCM_Module(Cancer_Type);

load(['Data_mat\INCM_Simulation\',Cancer_Type,'/1'])
load(['Data_mat\INCM_Significant\',Cancer_Type,'\Significant_Pair.mat'])
load(['Data_mat\Gene_Distance\',Cancer_Type,'.mat'])

GGS(:,4)=GGS(:,4)/2;
m=find(GGS(:,4)<0.01);
ggs=GGS(m,:);
ggs(:,6)=0;
Gene=unique(Net(:));
[a,b]=ismember(ggs(:,1:2),Gene);
for i=1:length(ggs)
    ggs(i,6)=Distance(b(i,1),b(i,2));
end
mg=unique(Mu(:,2));
for i=1:length(mg)
    mg(i,2)=length(find(Mu(:,2)==mg(i)));
end
[a,b]=ismember(ggs(:,1:2),mg(:,1));
ggs(:,7:8)=[mg(b(:,1),2) mg(b(:,2),2)];

ggs=sortrows(ggs,[-3,4]);
ggs((find(ggs(:,7)==1 | ggs(:,8)==1)),:)=[];
lll=100;
nodes=unique(ggs(1:lll,1:2));
[a,b]=ismember(Net,nodes);
c=find(a(:,1).*a(:,2)==1);
[LG,L1]=liuchuang_largest_component([Net(c,:);ggs(1:lll,1:2)]);
net=LG{find(L1==max(L1))};
gene=unique(net(:));
