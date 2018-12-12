load(file = "epigenetic_list/objects/epiGenesDF.rda")
load(file = "epigenetic_list/objects/exactab.rda")



###define the following functions for plotting pLI densities and pLIs as points, respectively
getpLIs <- function(genes, EMinput = c(TRUE, FALSE), name_to_use = c("exac", "gtex_rpkm", "gtex_rpm")){
  if (EMinput %in% c(TRUE)){
    if (name_to_use %in% c("exac")){
      plis <- epiGenesDF$exac.pLI[which(epiGenesDF$Gene_name_exac %in% genes & (!(epiGenesDF$Chr %in% c("chrX", "chrY"))))]
    } else if (name_to_use %in% c("gtex_rpkm")){
      plis <- epiGenesDF$exac.pLI[which(epiGenesDF$Gene_name_gtex_rpkm %in% genes & (!(epiGenesDF$Chr %in% c("chrX", "chrY"))))]
    } else {
      plis <- epiGenesDF$exac.pLI[which(epiGenesDF$Gene_name_gtex_rpm %in% genes & (!(epiGenesDF$Chr %in% c("chrX", "chrY"))))]
    }
  } else {
    plis <- exactab$pLI[which(exactab$gene %in% genes & (!(exactab$chr %in% c("X", "Y"))))]
  }
  if(length(which(is.na(plis))) > 0){
    plis <- plis[-which(is.na(plis))]
  }
  plis
}



plotpLIasPoints <- function(gene_names, insert = c(TRUE, FALSE), x_axis_position, 
                            plot_col, x_axis_lim, name_source, EM = c(TRUE, FALSE), yaxis_noise = 0, alpha){
  plis <- getpLIs(gene_names, EMinput = EM, name_to_use = "exac")
  if (insert %in% c(FALSE)){
    plot(rep(x_axis_position, length(plis)) + rnorm(length(plis), 0, 0.3), plis + rnorm(length(plis), 0, yaxis_noise),  
         ylim = c(-0.1, 1.1), xlim = x_axis_lim, yaxt = 'n', xaxt = 'n', pch = 19, col = plot_col,
         main = "", xlab = "", ylab = "pLI", cex.lab = 1.2, cex = 0.3, bty = 'u')
  } else {
    points(rep(x_axis_position, length(plis)) + rnorm(length(plis), 0, 0.3), plis + rnorm(length(plis), 0, yaxis_noise), 
           pch = 19, col = plot_col, cex = 0.3)
  }
}


##for the pli categories figure:
writers <- epiGenesDF$Gene_name_exac[which(epiGenesDF$Epi_function %in% c("Writer"))]
erasers <- epiGenesDF$Gene_name_exac[which(epiGenesDF$Epi_function %in% c("Eraser"))]
readers <- epiGenesDF$Gene_name_exac[which(epiGenesDF$Epi_function %in% c("Reader"))]
remodelers <- epiGenesDF$Gene_name_exac[which(epiGenesDF$Epi_function %in% c("Remodeler"))]
dualEMfunction <- epiGenesDF$Gene_name_exac[which(epiGenesDF$Epi_function %in% c("Writer;Reader", "Eraser;Reader", "Remodeler;Reader"))]
EMandTF <- exactab$gene[which(exactab$isEM %in% c(TRUE) & exactab$isTF %in% c(TRUE))]

##for the pli subcategories figure:
hmts <- epiGenesDF$Gene_name_exac[grep("HMT", epiGenesDF$Enzymatic.activity)]
hats <- epiGenesDF$Gene_name_exac[grep("HAT", epiGenesDF$Enzymatic.activity)]
dnmts <- epiGenesDF$Gene_name_exac[grep("DNMT", epiGenesDF$Enzymatic.activity)]
hdms <- epiGenesDF$Gene_name_exac[grep("HDM", epiGenesDF$Enzymatic.activity)]
hdacs <- epiGenesDF$Gene_name_exac[grep("HDAC", epiGenesDF$Enzymatic.activity)]
dnmerasers <- epiGenesDF$Gene_name_exac[grep("DNM_ERASE", epiGenesDF$Enzymatic.activity)]


readerdf <- epiGenesDF[which(is.na(epiGenesDF$Enzymatic.activity)), , drop = FALSE]
hmrs <- readerdf$Gene_name_exac[grep("HMR", readerdf$Reading.activity)]
hars <- readerdf$Gene_name_exac[grep("HAR", readerdf$Reading.activity)]
dnmrs <- readerdf$Gene_name_exac[grep("DNMR", readerdf$Reading.activity)]
dnumrs <- readerdf$Gene_name_exac[grep("DNUMR", readerdf$Reading.activity)]
hmrsonly <- hmrs[-which(hmrs %in% c(hars, dnmrs, dnumrs))]
harsonly <- hars[-which(hars %in% c(hmrs, dnmrs, dnumrs))]
dnmrsonly <- dnmrs[-which(dnmrs %in% c(hars, hmrs, dnumrs))]
dnumrsonly <- dnumrs[-which(dnumrs %in% c(hars, dnmrs, hmrs))]
dualreaders <- unique(c(hmrs[which(hmrs %in% c(hars, dnmrs, dnumrs))], hars[which(hars %in% c(hmrs, dnmrs, dnumrs))], 
                        dnmrs[which(dnmrs %in% c(hars, hmrs, dnumrs))], dnumrs[which(dnumrs %in% c(hars, dnmrs, hmrs))]))


em_excluding_tf <- exactab$gene[which(exactab$isEM %in% c(TRUE) & exactab$isTF %in% c(FALSE))]
tf_excluding_em <- exactab$gene[which(exactab$isTF %in% c(TRUE) & exactab$isEM %in% c(FALSE))]
allother <- exactab$gene[-which(exactab$isEM %in% c(TRUE) | exactab$isTF %in% c(TRUE))]

em_plis <- getpLIs(em_excluding_tf, EMinput = TRUE, name_to_use = "exac")
tf_plis <- getpLIs(tf_excluding_em, EMinput = FALSE)
allother_plis <- getpLIs(allother, EMinput = FALSE)
