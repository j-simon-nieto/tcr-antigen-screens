library(data.table)
library(readxl)
jv_twist=read_xlsx("PEPITOPE/s5_create_libraries/Download_TWIST_Oligopool JV B cell Screen Pt1-5_v0501.xlsx", skip = 11)
my_list=fread("PEPITOPE/s5_create_libraries/new_rna_library_all_patients_duplo_order.txt")
identical(jv_twist$`Oligo Seq`, my_list$CONCATENATED_SEQ)
identical(jv_twist$`Oligo Name`, my_list$variant_name)

table(stringr::str_split_i(jv_twist$`Oligo Name`,"_",1))/2
