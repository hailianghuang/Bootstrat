europe <- read.table("../pop_strat/cont_data.bim")

america <- read.table("../american/nimh_controls_qc_fordbgap_euro.bim")

merged <- merge(europe, america, by="V2")

check_pos <- merged$V4.x == merged$V4.y & merged$V1.x == merged$V1.y
sel_chr <- merged$V1.x<23

ATGC <- (merged$V5.x=="A" & merged$V6.x=="T" ) | (merged$V5.x=="G" & merged$V6.x=="C" ) |(merged$V5.x=="T" & merged$V6.x=="A" ) |(merged$V5.x=="C" & merged$V6.x=="G" )

ATGC <- ATGC | (merged$V5.y=="A" & merged$V6.x=="T" ) | (merged$V5.y=="G" & merged$V6.x=="C" ) |(merged$V5.y=="T" & merged$V6.x=="A" ) |(merged$V5.y=="C" & merged$V6.x=="G" )

snp <- (merged$V5.x=="A" | merged$V5.x=="T" |merged$V5.x=="G" |merged$V5.x=="C") & (merged$V5.y=="A" | merged$V5.y=="T" |merged$V5.y=="G" |merged$V5.y=="C")

selected <- check_pos & sel_chr & !ATGC & snp

write(as.character(merged$V2[selected]), file="snp_2_merge.txt")
merged <- merged[selected, ]
flip <- as.character(merged$V5.x) != as.character(merged$V5.y)
write(as.character(merged$V2[flip]), file="snp_2_flip.txt")
