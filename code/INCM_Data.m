function [Net,Mu,PL]=INCM_Data(Cancer_Type);
% Cancer_Type='LUAD';

a=textread(['Data/Genetic gene-gene network/',Cancer_Type,'_tumor_co-expression_GGIN.txt'],'%s');
a(1:20)=[];
l=length(a)/7;
Net=zeros(l,5);
for i=1:length(Net)
    Net(i,1)=str2num(a{7*(i-1)+1});
    Net(i,2)=str2num(a{7*(i-1)+2});
    Net(i,3)=str2num(a{7*(i-1)+5});
    Net(i,4)=str2num(a{7*(i-1)+6});
    Net(i,5)=str2num(a{7*(i-1)+7});
end

[a1,a2]=textread(['Data/Mutation_Individual/TCGA.',Cancer_Type,'.mutect.NonSilent.txt'],'%s%d');%,'headerlines',1);
a3{length(a1),1}=[];
for i=1:length(a1)
    str=a1{i};
    m=strfind(a1{1},'-');
    str=str(1:m(3)-1);
    a3{i,1}=str;
end  

m=unique(a3);
[c1,c2]=ismember(a3,m);
% for i=1:length(c2)
%     c2(i,2)=str2num(a2{i});
% end
c2(:,2)=a2;
c2=unique(c2,'rows');
c2=sortrows(c2,1);
Mu=c2;
PL=unique(a3);
