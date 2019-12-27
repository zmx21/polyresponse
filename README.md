## Polygenic modelling of treatment effect heterogeneity
The standard workflow should be:
1. Calculate the quality metric of each SNP
    - The scripts assume that genotype files are in bgen format, and organized by chromosome (as from UKBB). 
    - [CalcVariantQuality.R](./Variant_Quality/CalcVariantQuality.R) in [Variant_Quality](./Variant_Quality) is the function to calculate quality metircs of each SNP. This will save text files which would contain quality metrics (Info, MAF, etc) and will be read in subsequent analysis. Ensure that QCTOOL is installed. 
2. Identify proxy SNPs (pharmacomimetic variants)
    - [GenePhenotypeAssociation.R](./Gene_Phenotype_Association/GenePhenotypeAssociation.R) in [Gene_Phenotype_Association](Gene_Phenotype_Association) is the function to identify proxy SNPs. This will identify a set of independent SNPs, ranked by their association to phenotype. Output will contain the rsids of each SNP, as well as the beta coefficient. 
3. Calculate genomewide GxG Interactions (individual SNPs against proxy SNP(s))
    - [GxG_Interactions.R](./Target_Gene_Interactions/GxG_Interactions.R) in [Target_Gene_Interactions](./Target_Gene_Interactions) is the main script which calculates GxG interactions. Proxy SNPs should be specified, seperated by commas in a string. Filtering criteria for SNPs to analyze should also be included. 
4. Construct Random Forest of Interaction Trees. 
    - [RandomForest.R](./Recursive_Partitioning/RandomForest.R) in [Recursive_Partitioning](./Recursive_Partitioning) is the main function which constructs random forest. Arguments can specify the dataset to use, the p-value threshold for predictors to include, and the minimum node size of the interactions trees. Output are saved as .rds files. 
