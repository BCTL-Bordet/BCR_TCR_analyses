# BCR_TCR_analyses
Functions to compute BCR/TCR diversity measures from MiXCR outputs. 
The scripts have been tested on the R software (v4.0.5 and v4.2.1).
The two simulated samples were used to produce Supplementary Figure 1 in Rediti M, et al. Immunological and clinicopathological features predict HER2-positive breast cancer prognosis in the neoadjuvant NeoALTTO and CALGB 40601 randomized trials. Nat Commun. 2023. doi: 10.1038/s41467-023-42635-2. PMID: 37923752.

First, download the files list_examples.RDS and tot_number_reads.RDS.
- list_examples.RDS is a list with 2 simulated samples containing, in this case, only IG information, and it shows how MiXCR output should be formatted to run the script (explanations for the columns are available in /R/scripts/script_for_example.R)
- tot_number_reads.RDS is dataframe with the total number of reads mapping to genes for the two samples, required for the normalization of BCR/TCR number of reads


The functions to compute BCR/TCR measures are stored in /R/functions.

The script to run the example is stored in /R/scripts and requires the content of /R/functions, "list_examples.RDS" and "tot_number_reads.RDS".

The files "measures.RDS" and "df_measures.RDS" represent the expected outputs created with script_for_example.R and can be checked to verify the results obtained with the given examples.


Total runtime for script_for_example.R on a computer with 32GB RAM, Apple M1 Max CPU is ~1.6-2 seconds.


**Please cite our work when using this script.**

Citation:
Rediti M, Fernandez-Martinez A, Venet D, Rothé F, Hoadley KA, Parker JS, Singh B, Campbell JD, Ballman KV, Hillman DW, Winer EP, El-Abed S, Piccart M, Di Cosimo S, Symmans WF, Krop IE, Salgado R, Loi S, Pusztai L, Perou CM, Carey LA, Sotiriou C. Immunological and clinicopathological features predict HER2-positive breast cancer prognosis in the neoadjuvant NeoALTTO and CALGB 40601 randomized trials. Nat Commun. 2023 Nov 3;14(1):7053. doi: 10.1038/s41467-023-42635-2. PMID: 37923752.
