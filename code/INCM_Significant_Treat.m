function INCM_Significant_Treat(Cancer_Type);

GGS=[];
filename=dir(['Data_mat/INCM_Significant/',Cancer_Type]);
filename(1:2)=[];
for i=1:length(filename)-1
    load(['Data_mat/INCM_Significant/',Cancer_Type,'/',num2str(i)])
    GGS=[GGS;full(BGGS)];
end
load(['Data_mat/INCM_Simulation/',Cancer_Type,'/1'])
save(['Data_mat/INCM_Significant/',Cancer_Type,'/Significant_Pair'],'GGS','Net','mu','Mu','Patient')