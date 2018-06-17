The data and code for paper “Individualized genetic network analysis reveals new therapeutic tumor vulnerabilities”.

The main code is INCM_Simulation_Cluster.m, which is used to calculate the INCM value. 

Data Description:
Genetic gene-gene network: 
The raw gene-gene interaction data for 14 cancer types.

Mutation_Individual:
The mutation data for individual sample for 14 cancer types.

Gene Set:
List of the functional gene sets

Drug_Gene:
Drug-Gene interactions.

Wiki Pathway validation:
Pathway enrichment data for 14 cancer types.


Code Description:
INCM_Data.m    
Treat the genetic network data and the individual mutation data.

Gene_Distance_Cluster.m   
Calculate the Gene-Gene distance on the genetic network for each cancer type.

largest_component.m  
Calculate the largest component of the network.

INCM_Simulation_Cluster.m  
Calculate the INCM value for the raw data as well as the simulation data. 

INCM_Significant.m   
Calculate the significant of INCM value for each gene-gene pair.

INCM_Significant_Treat_Treat.m   
Treat the significant result.

INCM_Module.m  
Calculate the significant module.


It should be noted that we save and recall many intermediate results in our code for the computational complexity for this project. And we also provide the intermediate results in the folder ‘Data_mat’. 

Data_mat:
Gene_Distance:
The distance between gene pair in the genetic network. 

INCM_Simulation:
The INCM value for each gene pair. The files are too large to upload here. 

INCM_Significant:
The significant for each gene pair. 
