 function INCM_Simulation_Cluster(Cancer_Type);

[Net,Mu,PL]=INCM_Data(Cancer_Type);
Net=Net(find(Net(:,4)<0.05),:);
Net=Net(:,1:2);
Net(find(Net(:,1)-Net(:,2)==0),:)=[];
[LG,L]=largest_component(Net);
Net=LG{find(L==max(L))};
Mu=unique(Mu,'rows');
mu=unique(Mu(:,1));
for i=1:length(mu)
    mu(i,2)=length(find(Mu(:,1)==mu(i)));
end
x=4*mean(mu(:,2));
m=find(mu(:,2)>=x);
if ~isempty(m)
    for i=1:length(m)
        x=find(Mu(:,1)==mu(m(i),1));
        Mu(x,:)=[];
    end
end
x=ismember(Mu(:,2),unique(Net));
Mu=Mu(x,:);
GM=unique(Mu(:,2));
for i=1:length(GM)
    GM(i,2)=length(find(Mu(:,2)==GM(i)));
end
GM(:,3)=cumsum(GM(:,2))/sum(GM(:,2));
mu=unique(Mu(:,1));

matlabpool local 12

parfor I=1:10000    
    gravity_model(I,Net,Mu,GM,Cancer_Type,PL);
end
matlabpool close

function gravity_model(I,Net,Mu,GM,Cancer_Type,PL);
load(['Data_mat/Gene_Distance/',Cancer_Type]);
rand('state',sum(clock)*I);
Patient=PL;
Gene=unique(Net(:));
PG=sparse(zeros(length(Gene),length(Gene)));
mu=unique(Mu(:,1));
MP{length(mu),1}=[];
if I~=1
    for j=1:length(mu)
        l=length(find(Mu(:,1)==mu(j)));
        lx=rand(l,1);
        lrg=zeros(l,1);
        for jj=1:length(lx)
            mm=find(GM(:,3)>=lx(jj),1,'first');
            lrg(jj)=GM(mm,1);
        end
        MP{j,1}=lrg;
    end    
else
    for j=1:length(mu)
        x=find(Mu(:,1)==mu(j));
        MP{j}=Mu(x,2);
    end
end
m=cell2mat(MP);
GM1=unique(m);
for i=1:length(GM1)
    GM1(i,2)=length(find(m==GM1(i)));
end

for k=1:length(MP)
    pg=sparse(zeros(length(Gene),length(Gene)));
    [a,b]=ismember(MP{k},Gene);
    pg(b,b)=1;
    pg=triu(pg);
    PG=PG+pg;
end
PG(logical(eye(size(PG))))=0;
g=Gene;
g(:,2)=0;
[a,b]=ismember(g(:,1),GM1(:,1));
g(a,2)=GM1(b(a),2);
g=sqrt(g(:,2));
gg=g*g';
gg=sparse(triu(gg));
Distance=sparse(triu(Distance));
PG=PG./gg./Distance.^2;
PG(isnan(PG))=0;

if I==1
    save(['Data_mat/INCM_Simulation/',Cancer_Type,'/1'],'Net','PG','mu','Mu','Patient')
else
    save(['Data_mat/INCM_Simulation/',Cancer_Type,'/',num2str(I)],'PG')
end



        