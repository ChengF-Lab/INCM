function Gene_Distance_Cluster(Cancer_Type);

[Net,Mu]=INCM_Data(Cancer_Type);
Net=Net(find(Net(:,4)<0.05),:);
Net=Net(:,1:2);
Net(find(Net(:,1)-Net(:,2)==0),:)=[];
[LG,L]=largest_component(Net);
Net=LG{find(L==max(L))};
Gene=unique(Net(:));
Distance{length(Gene),1}=[];
matlabpool local 2
parfor II=1:length(Gene)
    Distance{II}=Node_Distance(Gene(II),Gene',Net);
end
Distance=cell2mat(Distance);
save(['Data_mat/Gene_Distance/',Cancer_Type],'Distance','Net')
matlabpool close


function D=Node_Distance(n1,n2,Net);
D=zeros(1,length(n2));
n2=n2';
n2(:,2)=(1:length(n2))';
n2(find(n2(:,1)==n1),:)=[];
k=1;
while ~isempty(n2)
    a=ismember(Net,n1);
    n=unique([Net(a(:,1),2);Net(a(:,2),1)]);
    [m,m1,m2]=intersect(n2(:,1),n);
    D(n2(m1,2))=k;
    n2(m1,:)=[];
    n1=unique([n1;n]);
    Net(find(a(:,1)+a(:,2)~=0),:)=[];
    k=k+1;
end


