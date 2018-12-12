##########EM Complex analysis
baf_complex <- c("SMARCA4", "ARID2", "PBRM1", "SMARCC1", "SMARCC2", "SMARCD1", "SMARCD2", 
                 "SMARCD3", "SMARCB1", "SMARCE1", "ACTL6A", "ACTL6B", "ACTB")

pbaf_complex <- c("SMARCA2", "SMARCA4", "ARID1A", "SMARCC1", "SMARCC2", "SMARCD1", "SMARCD2", 
                  "SMARCD3", "SMARCB1", "SMARCE1", "ACTL6A", "ACTL6B", "ACTB")

chrach_complex <- c("SMARCA5", "BAZ1A", "POLE3", "CHRAC1")

nurf_complex <- c("SMARCA1", "BPTF", "RBBP7", "RBBP4")

nurd_complex <- c("CHD3", "CHD4", "MBD2", "MBD3", "MTA1", "MTA2", "MTA3", "HDAC1", "HDAC2",
                  "GATAD2A", "GATAD2B", "KDM1A", "RBBP4", "RBBP7")

ino80_complex <- c("INO80", "ACTR5", "ACTR8", "INO80B", "INO80C", "RUVBL1", "RUVBL2", "ACTL6A", 
                   "TFPT", "NFRKB", "MCRS1", "UCHL5", "INO80E", "INO80D")

srcap_complex <- c("SRCAP", "RUVBL1", "RUVBL2", "ACTL6A", "ACTR6", "ZNHIT1", "YEATS4", "DMAP1", "VPS72")

trrap_complex <- c("EP400", "RUVBL1", "RUVBL2", "ACTL6A", "YEATS4", "DMAP1", "VPS72", "BRD8", 
                   "TRRAP", "KAT5", "MORF4L1", "MORF4L2", "MEAF6", "MRGBP", "EPC1", "EPC2", "ING3")

staga_complex <- c("KAT2B", "TRRAP", "TADA1", "TADA3", "SUPT3H", "SUPT7L", "TAF5L", "TAF6L", "TAF9", "TAF10", "TAF12",
                   "SGF29", "STAF46", "STAF55", "STAF60", "SAP130")

pcaf_complex <- c("KAT2B", "TRRAP", "TADA2A", "TADA3", "SPTY2D1", "TAF5L", "TAF6L", "TAF9", "TAF10", "TAF12")

tftc_complex <- c("KAT2B", "TRRAP", "TADA3", "SUPT3H", "TAF2", "TAF4", "TAF5", "TAF5L", "TAF6", "TAF6L", "TAF9", "TAF10", "TAF12", 
                  "SAP130")

kmt2ab_complex <- c("KMT2A", "KMT2B", "ASH2L", "RBBP5", "WDR5", "DPY30", "MEN1", "HCFC1", "HCFC2")

kmt2cd_complex <- c("KMT2C", "KMT2D", "ASH2L", "RBBP5", "WDR5", "DPY30", "PAXIP1", "PAGR1", "NCOA6", "KDM6A")

kmt2fg_complex <- c("SETD1A", "SETD1B", "ASH2L", "RBBP5", "WDR5", "DPY30", "CXXC1", "WDR82")

prc2_complex <- c("EZH1", "EZH2", "EED", "SUZ12", "RBBP4", "RBBP7", "AEBP2", "PHF1", "MTF2", "PHF19", "JARID2")

prc1_complex <- c(paste0("PCGF", c(1,2,3,5,6)), "BMI1", paste0("CBX", c(2,4,6,7,8)), paste0("PHC", c(1,2,3)), "RING1", "RNF2")

prc1_complex_non_canonical <- c("KDM2B", "RYBP", "YAF2","RING1", "RNF2", paste0("PCGF", c(1,2,3,5,6)), "BMI1", "L3MBTL2")

corest_complex <- c("KDM1A", "RCOR1", "ZNF217", "HDAC1", "HDAC2", "RBBP4", "RBBP7")

swi_ind_3_complex <- c("SAP30", "SAP18", "SUDS3", "SIN3A", "SIN3B", "HDAC1", "HDAC2", "RBBP4", "RBBP7")



########
all_complex_subunits <- c(baf_complex, pbaf_complex, chrach_complex, nurf_complex, nurd_complex, 
                          ino80_complex, srcap_complex, trrap_complex, staga_complex, pcaf_complex, 
                          tftc_complex, kmt2ab_complex, kmt2cd_complex, kmt2fg_complex, prc1_complex, 
                          prc1_complex_non_canonical, prc2_complex, corest_complex, swi_ind_3_complex)

all_accessory_subunits <- unique(all_complex_subunits[-which(all_complex_subunits %in% epiGenesDF$Gene_name)])
all_em_subunits <- unique(all_complex_subunits[which(all_complex_subunits %in% epiGenesDF$Gene_name)])


em_plis <- getpLIs(em_excluding_tf, EMinput = TRUE)
tf_plis <- getpLIs(tf_excluding_em, EMinput = FALSE)

accessory <- exactab$gene[which(exactab$isEM %in% c(FALSE) & exactab$isTF %in% c(FALSE) 
                                & exactab$gene %in% all_accessory_subunits)]
accessory_plis <- getpLIs(accessory, EMinput = FALSE)

quartz(file = "em_accessory_plis_vs_tf.pdf", height = 3, width = 3, type = "pdf")
plot(density(em_subunit_plis, from = 0, to = 1), 
     col = 2, ylab = "Density", 
     xlab = "pLI", 
     xaxt ='n', yaxt = 'n', cex.lab = 1.4, main = "", bty = 'l', lwd = 2.5)
lines(density(accessory_plis, 
              from = 0, to = 1), 
      col = 1, lwd = 2.5)
lines(density(tf_plis, from = 0, to = 1), col = "forestgreen", lwd = 2.5)
axis(2, at = c(0.5, 1.5), cex.axis = 1.4)
axis(1, at = c(0.1, 0.5, 0.9), cex.axis = 1.4)
legend <- legend("top", legend = c("EM genes", "TF genes", "EM accessory\nsubunits"), 
                 col = c(2, 3, 1), lty = "solid", cex = 0.5, bty = 'n', lwd = 1.4)
polygon(c(0.9, 0.9, 1.05, 1.05), c(0, 2.6, 2.6, 0), col = adjustcolor("gray80", alpha.f = 0.5), lty = 0)


dev.off()

quartz(file = "em_pli_complexes.pdf", width = 5, height = 5, type = "pdf")
plot(1, col = rgb(1,1,1), xlim = c(2, 38), ylab = "pLI", xlab = "", xaxt = 'n',
     yaxt = 'n', main = "", bty = 'u', ylim = c(-0.1,1.1))
polygon(c(0, 0, 42, 42), c(0.9, 1.2, 1.2, 0.9), 
        col = adjustcolor("gray80", alpha.f = 0.5), lty = 0)
####non EM subunits
plotpLIasPoints(baf_complex[-which(baf_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 2, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(pbaf_complex[-which(pbaf_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 4, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(chrach_complex[-which(chrach_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 6, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(nurf_complex[-which(nurf_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 8, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(nurd_complex[-which(nurd_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 10, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(ino80_complex[-which(ino80_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 12, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(srcap_complex[-which(srcap_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 14, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(trrap_complex[-which(trrap_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 16, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(staga_complex[-which(staga_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 18, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(pcaf_complex[-which(pcaf_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 20, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(tftc_complex[-which(tftc_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 22, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(kmt2ab_complex[-which(kmt2ab_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 24, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(kmt2cd_complex[-which(kmt2cd_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 26, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(kmt2fg_complex[-which(kmt2fg_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 28, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(prc2_complex[-which(prc2_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 30, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(prc1_complex[-which(prc1_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 32, rgb(0, 0, 0, 0.7), EM = FALSE)


plotpLIasPoints(prc1_complex_non_canonical[-which(prc1_complex_non_canonical %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 34, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(corest_complex[-which(corest_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 36, rgb(0, 0, 0, 0.7), EM = FALSE)

plotpLIasPoints(swi_ind_3_complex[-which(swi_ind_3_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 38, rgb(0, 0, 0, 0.7), EM = FALSE)

####EM subunits
plotpLIasPoints(baf_complex[which(baf_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 2, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(pbaf_complex[which(pbaf_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 4, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(chrach_complex[which(chrach_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 6, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(nurf_complex[which(nurf_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 8, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(nurd_complex[which(nurd_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 10, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(ino80_complex[which(ino80_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 12, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(srcap_complex[which(srcap_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 14, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(trrap_complex[which(trrap_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 16, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(staga_complex[which(staga_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 18, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(pcaf_complex[which(pcaf_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 20, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(tftc_complex[which(tftc_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 22, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(kmt2ab_complex[which(kmt2ab_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 24, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(kmt2cd_complex[which(kmt2cd_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 26, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(kmt2fg_complex[which(kmt2fg_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 28, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(prc2_complex[which(prc2_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 30, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(prc1_complex[which(prc1_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 32, rgb(1,0,0,0.8), EM = TRUE)


plotpLIasPoints(prc1_complex_non_canonical[which(prc1_complex_non_canonical %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 34, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(corest_complex[which(corest_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 36, rgb(1,0,0,0.8), EM = TRUE)

plotpLIasPoints(swi_ind_3_complex[which(swi_ind_3_complex %in% epiGenesDF$Gene_name)], 
                insert = TRUE, 38, rgb(1,0,0,0.8), EM = TRUE)

axis(1, at = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38), 
     labels = c("BAF", "PBAF", "CHRACH", "NURF", "NURD", "INO80",
                "SRCAP", "TRRAP", "STAGA", "PCAF", "TFTC", "KMT2A/B", "KMT2C/D", "KMT2F/G", "PRC2", 
                "cPRC1", "ncPRC1", "COREST", "SWI IND 3"),
     las = 2, cex.axis = 0.9)
axis(2, at = c(0.1, 0.9), cex.axis = 1.2)
abline(v=3, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=5, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=7, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=9, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=11, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=13, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=15, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=17, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=19, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=21, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=23, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=25, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=27, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=29, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=31, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=33, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=35, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=37, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))

#plot(1, col = rgb(1,1,1), bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
#legend <- legend("topleft", legend = c("EM subunits", "accessory subunits"), 
#                 col = c(rgb(1,0,0,0.7), rgb(0, 0.8 , 1, 0.8)), pch = 19, bty = 'n')
#par(op)
dev.off()


############disease candidates within accessory complex subunits
omim <- read.delim("mdem/epigenetic_list/data/allgenesomim.txt")

accessory <- exactab$gene[which(exactab$gene %in% accessory & (!(exactab$chr %in% c("X", "Y"))))]

accessory_plis <- sapply(accessory, function(xx) exactab$pLI[which(exactab$gene == xx)])

accessory_pli_greater_than_0.9 <- accessory[which(accessory_plis > 0.9)]

accessory_subunit_disease_associations <- lapply(accessory_pli_greater_than_0.9, function(xx) {
  omim_entry <- grep(xx, omim$Gene.Symbols)
  omim$Phenotypes[omim_entry]
})



accessory_subunit_has_omim_entry <- rep(FALSE, length(accessory_subunit_disease_associations))
accessory_subunit_has_omim_entry[c(1, 3, 19, 30, 37, 45)] <- TRUE #based on the above list and after manually curating OMIM


######
complex_list <- list(baf_complex, pbaf_complex, chrach_complex, nurf_complex, nurd_complex, 
                     ino80_complex, srcap_complex, trrap_complex, staga_complex, pcaf_complex, 
                     tftc_complex, kmt2ab_complex, kmt2cd_complex, kmt2fg_complex, prc1_complex, 
                     prc1_complex_non_canonical, prc2_complex, corest_complex, swi_ind_3_complex)


complexes_accessory_plis <- lapply(complex_list, function(xx) 
  getpLIs(xx[-which(xx %in% epiGenesDF$Gene_name)], EMinput = FALSE))

complexes_accessory_plis_greaterThan0.9 <- lapply(complexes_accessory_plis, 
                                                  function(xx) round(length(which(xx > 0.9))/length(xx), 2))

complexes_EM_plis <- lapply(complex_list, function(xx) 
  getpLIs(xx[which(xx %in% epiGenesDF$Gene_name)], EMinput = FALSE))

complexes_EM_plis_greaterThan0.9 <- lapply(complexes_EM_plis, 
                                           function(xx) round(length(which(xx > 0.9))/length(xx), 2))

median(unlist(complexes_accessory_plis_greaterThan0.9))

median(unlist(complexes_EM_plis_greaterThan0.9))




########
EM_complex_subunits <- data.frame(gene_name = all_complex_subunits)
EM_complex_subunits$complex_name <- NA
EM_complex_subunits$subunit_type <- NA

EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% baf_complex)] <- "BAF"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% pbaf_complex)] <- "PBAF"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% chrach_complex)] <- "CHRACH"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% nurf_complex)] <- "NURF"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% nurd_complex)] <- "NURD"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% nurd_complex)] <- "NURD"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% ino80_complex)] <- "INO80"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% srcap_complex)] <- "SRCAP"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% trrap_complex)] <- "TRRAP"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% staga_complex)] <- "STAGA"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% pcaf_complex)] <- "PCAF"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% tftc_complex)] <- "TFTC"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% kmt2ab_complex)] <- "KMT2A/B"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% kmt2cd_complex)] <- "KMT2C/D"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% kmt2fg_complex)] <- "KMT2C/D"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% kmt2fg_complex)] <- "KMT2F/G"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% prc1_complex)] <- "canonical PRC1"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% prc1_complex_non_canonical)] <- "non canonical PRC1"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% prc2_complex)] <- "PRC2"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% corest_complex)] <- "CoREST"
EM_complex_subunits$complex_name[which(EM_complex_subunits$gene_name %in% swi_ind_3_complex)] <- "Swi Independent 3"

EM_complex_subunits$subunit_type[which(EM_complex_subunits$gene_name 
                                       %in% all_accessory_subunits)] <- "accessory"

EM_complex_subunits$subunit_type[which(EM_complex_subunits$gene_name 
                                       %in% all_em_subunits)] <- "EM"


write_csv(EM_complex_subunits, "data/EM_complexes.csv")



