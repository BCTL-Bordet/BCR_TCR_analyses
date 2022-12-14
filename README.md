# BCR_TCR_analyses
Functions to compute BCR/TCR diversity measures, with examples, from MiXCR outputs. 
The scripts have been tested on R (v4.0.5 and v4.2.1).
The two examples represent the toy samples created to produce Supplementary Figure 1 in PMID:...

First, download the files list_examples.RDS and tot_number_reads.RDS.
- list_examples.RDS is a list with 2 toy samples containing, in this case, only IG information, and it shows how MiXCR output should be formatted to run the script (explanations for the columns are available in /R/scripts/script_for_example.R)
- tot_number_reads.RDS is dataframe with the total number of reads mapping to genes for the two samples, required for the normalization of BCR/TCR number of reads


The functions to compute BCR/TCR measures are stored in /R/functions.

The script to run the example is stored in /R/scripts and requires the content of /R/functions, "list_examples.RDS" and "tot_number_reads.RDS".

The files "measures.RDS" and "df_measures.RDS" represent the expected outputs created with script_for_example.R and can be checked to verify the results obtained with the given examples.


Total runtime for script_for_example.R on a computer with 32GB RAM, Apple M1 Max CPU is ~1.6-2 seconds.
