# Rcode_Weibull_measurement_error_and_misclassification
The folders and the R code help to reproducing Tables 1-5 and Tables S1-S3 in the article "Weibull regression with both measurement error and misclassification in covariates" by Zhiqiang Cao and Man Yu Wong (under review).

For questions or comments about the code, please contact Zhiqiang Cao zcaoae@connect.ust.hk.
You will need to change the directory to use the example code in script R code. This folder includes the following functions:

1.main_simulation_for_tables_1-4 and S1-S2.R can reproduce simulation results in Tables 1-4 of main text and Tables S1-S2 of supplementary materials;

2.main_simulation_table5.R can reproduce simulation results in Table 5 of main text;

3.supp_simulation_table_S3.R can produce similar results in Table S3 of supplementary materials;

4.real_data_male_table6.R can produce similar results in the left panel of Table 6 of main text using masked EPIC-InterAct data;

5.real_data_female_table6.R can produce similar results in the right panel of Table 6 of main text using masked EPIC-InterAct data;

6.score_fisher_setup_new.R are source programs used in simulation studies;

7.score_fisher_setup_real_data_ana.R are source programs used in real data analysis;

8.epic_mask.RData is the masked EPIC-InterAct data used in real_data_male_table6.R and real_data_female_table6.R to produce similar estimation results in Table 6 of main text.

