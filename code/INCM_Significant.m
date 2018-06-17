function Genetic_Gravity_Significant(Cancer_Type);
tic;
load(['Data_mat/INCM_Simulation/',Cancer_Type,'/1'])
Gene=unique(Net(:));
[a,b]=find(PG~=0);
m=find(PG~=0);
GPG=PG(m);
GGS=[Gene(a) Gene(b) GPG];

if mod(length(m),10000)==0
    bin=[0:10000:length(m)];
else
    bin=[0:10000:length(m) length(m)];
end
ll=length(bin)-1;
BGGS{length(bin)-1,1}=[];

matlabpool local 12
parfor I=1:ll
    significant_block_calculation(I,m,bin,GPG,GGS,Cancer_Type);
end

GN_S=zeros(length(Gene),10000);
for j=1:10000
    load(['Data_mat/INCM_Simulation/',Cancer_Type,'/',num2str(j)])
    pg=PG+PG';
    GN_S(:,j)=sum(pg,2);
end

save(['Data_mat/INCM_Significant/',Cancer_Type,'/Gravity_Gene'],'GN_S','Net')

matlabpool close
toc;
    
function BGGS=significant_block_calculation(i,m,bin,GPG,GGS,Cancer_Type);
x=m(bin(i)+1:bin(i+1));
rgp=zeros(length(x),10000);
rgp(:,1)=GPG(bin(i)+1:bin(i+1));
for j=2:10000
    load(['Data_mat/Genetic_Gravity_Simulation/',Cancer_Type,'/',num2str(j)])
    rgp(:,j)=PG(x);
end
for j=1:length(x)
    p(j,1)=(length(find(rgp(j,2:end)>rgp(j,1)))+0.5*length(find(rgp(j,2:end)==rgp(j,1))))/9999;
    z(j,1)=(rgp(j,1)-mean(rgp(j,2:end)))/std(rgp(j,2:end));
end
BGGS=[GGS((bin(i)+1:bin(i+1)),:) p z];
load(['Data_mat/INCM_Simulation/',Cancer_Type,'/1'])
save(['Data_mat/INCM_Significant/',Cancer_Type,'/',num2str(i)],'BGGS','Net','mu','Patient')
    