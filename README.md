# Polygenic modelling of treatment effect heterogeneity
Xu, ZM, Burgess, S. Polygenic modelling of treatment effect heterogeneity. Genetic Epidemiology. 2020; 1– 12. https://doi.org/10.1002/gepi.22347
## Data Availability 
Summary statistics of genome-wide interaction study and RFIT files availiable at: https://doi.org/10.5281/zenodo.3975035
## Usage
1. Calculate the quality metric of each SNP
    - The scripts assume that genotype files are in bgen format, and organized by chromosome (as from UKB) 
    - [CalcVariantQuality.R](./Variant_Quality/CalcVariantQuality.R) in [Variant_Quality](./Variant_Quality) calculate quality metrics for each SNP. Text files containing quality metrics (Info, MAF, etc) will be saved and read in subsequent analysis. Ensure that QCTOOL is installed 
2. Evaluate pharmacomimetic variants (genetic proxy for treatment)
    - [GenePhenotypeAssociation.R](./Gene_Phenotype_Association/GenePhenotypeAssociation.R) in [Gene_Phenotype_Association](Gene_Phenotype_Association) calculates the marginal effect of each pharmacomimetic variant for their association to the specified phenotype. Output will contain the rsids of each SNP and beta coefficient. 
3. Conduct a genome-wide GxG Interaction Study (effect moderating variant against pharmacomimetic variant/score)
    - [GxG_Interactions.R](./Target_Gene_Interactions/GxG_Interactions.R) in [Target_Gene_Interactions](./Target_Gene_Interactions) is the main script which searches genome-wide for GxG interactions between candidate effect moderating variants and the pharmacomimetic variant/score. Pharmacomimetic variants should be specified, and seperated by commas in a string if multiple variants are specified. If multiple variants are specified, a polygenic pharmacomimetic score would be conducted. Filtering criteria for candidate effect moderating variants to be included in the analysis should also be specified
4. Construct Random Forest of Interaction Trees (composite polygenic effect modifier)
    - [RandomForest.R](./Recursive_Partitioning/RandomForest.R) in [Recursive_Partitioning](./Recursive_Partitioning) is the main function which constructs the RFITs. Arguments can specify the dataset to use, the p-value threshold for predictors to include, and the minimum node size of the interactions trees. Output are saved as .rds files and stored in seperate directories (for each p-value and minimum node size setting combination)
