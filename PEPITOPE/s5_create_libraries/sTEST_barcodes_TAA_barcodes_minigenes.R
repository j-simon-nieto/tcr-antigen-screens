library(data.table)
library(ggVennDiagram)
library(readxl)

# small library
taa=fread("PEPITOPE/s5_create_libraries/resources_zm/v3_minigen_antigens_4k.csv", header = F)
neo=fread("PEPITOPE/s5_create_libraries/new_rna_library_all_patients_duplo.txt")

intersect_list=list(neo=neo$BARCODE, taa=taa$V3)
pdf("PEPITOPE/s5_create_libraries/v3_library_intersect.pdf", width = 4, height=4)
ggVennDiagram(intersect_list, label_color = "white")
dev.off()

intersect(taa$V3,neo$BARCODE)
nchar(taa$V3[1])
nchar(neo$BARCODE[1])

# big library 1
taa_big1=fread("PEPITOPE/s5_create_libraries/resources_zm/Large_TAA2.0_Oligo library - shared antigens V3 (93nt minigenes)_CSV_FINAL_LIBRARY.csv", header = T, skip = 4)

intersect_list=list(neo=neo$BARCODE, taa_big1=taa_big1$`Oligo barcode`)
pdf("PEPITOPE/s5_create_libraries/large_TAA_1_library_intersect.pdf", width = 4, height=4)
ggVennDiagram(intersect_list, label_color = "white")
dev.off()

intersect(taa_big1$`Oligo barcode`,neo$BARCODE)
nchar(taa_big1$`Oligo barcode`[1])
nchar(neo$BARCODE[1])

# big library 2
taa_big2=readxl::read_xlsx("PEPITOPE/s5_create_libraries/resources_zm/Large TAA_Barcode list_TAA2.0.xlsx", col_names = F)
names(taa_big2)=c("Oligo barcode", "x")
taa_big2=taa_big2[taa_big2$x=="Shared library",]

intersect_list=list(neo=neo$BARCODE, taa_big2=taa_big2$`Oligo barcode`)
pdf("PEPITOPE/s5_create_libraries/large_TAA_2_library_intersect.pdf", width = 4, height=4)
ggVennDiagram(intersect_list, label_color = "white")
dev.off()

intersect(taa_big2$`Oligo barcode`,neo$BARCODE)
nchar(taa_big2$`Oligo barcode`[1])
nchar(neo$BARCODE[1])