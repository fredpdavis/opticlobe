################################################################################
#
# analyzeOpticLobeExpr.R - analyze TAPIN-seq measurements of gene expression
#                          in individual cell types.
#
# Fred P. Davis, fredpdavis@gmail.com. 2015-2018
#
################################################################################

main <- function( data,
                  loadPreprocData = FALSE,
                  returnData      = FALSE,
                  makeFigs        = FALSE ) {

# Purpose: Run all code to load data and generate figures and tables


   if (missing(data) & loadPreprocData) {
# Load preprocessed data including all raw data and inferState results

      data <- readRDS(paste0(dat$specs$outDir,"/full_dataset.rds"));

   } else {
# Load raw data and call inferState 

      if (missing(data)){
         data <- list()
         data$specs <- setSpecs()
      }

      if (! "expr" %in% names(data)) {
         data$expr  <- loadExpr( data$specs,
                                 printTxExprMat    = TRUE,
                                 printTxCleanMat   = TRUE,
                                 printGeneExprMat  = TRUE,
                                 printGeneCleanMat = TRUE)
         data$specs$sampleLists <- makeSampleLists(data)
         if (returnData) return(data)
      }

      if (! "rawpon" %in% names(dat$expr)) {
         data$expr$rawpon <- list()
         for (mode in c("cells", "drivers", "samples")) {
            data$expr$rawpon[[mode]] <- prepRawPonMats(data, mode=mode) }
         if (returnData) return(data)
      }
   }

   if (makeFigs) { makeFigs(data = data) }

   return(data)

}



makePaperNumbers <- function( dat ) {
# Purpose: generate all numbers used in the paper.

   print("*********************************************************")
   print("R3.1. OVERVIEW OF LIBRARIES")
   print("---------------------------")

   print(paste0("Control samples: n=",
        sum(dat$specs$sampleInfo$driverType == "control")))

   print(paste0("Real samples: n=",
        sum(dat$specs$sampleInfo$driverType != "control")))

   print(paste0("INTACT samples: n=",
        sum(dat$specs$sampleInfo$driverType != "control" &
            dat$specs$sampleInfo$intact.protocol == "INTACT")))

   print(paste0("TAPIN samples: n=",
        sum(dat$specs$sampleInfo$driverType != "control" &
            dat$specs$sampleInfo$intact.protocol == "TAPIN")))

   print(paste0("Dissected samples: n=",
        sum(dat$specs$sampleInfo$driverType == "dissected")))



   qcPass <- dat$specs$sampleInfo[
               dat$specs$sampleInfo$sample_name %in%
               dat$specs$sampleLists$good$allCriteria,]

   print(paste0("TOTAL NUMBER OF QC PASS DRIVERS: ",
         length(unique(qcPass$driver))))

   qcPassType <- list()
   for (type in c("pureCellType", "comboCellTypes", "broad", "dissected")) {
      qcPassType[[type]] <- unique(qcPass$celltype[qcPass$driverType == type])

      print(paste0("Number of qcPass ",type," CELLS n=",
                   length(qcPassType[[type]])))

      print(paste0("Number of qcPass ",type," SAMPLES n=",
         nrow(qcPass[qcPass$driverType == type,])))

      curNumDriv <- length(unique(qcPass$driver[qcPass$driverType == type]))
      print(paste0("-> from ",curNumDriv," drivers"))

      if (type %in% c("pureCellType", "comboCellTypes")) {
         t.muscle <- qcPassType[[type]][grepl("uscle", qcPassType[[type]])]

         t.mb <- qcPassType[[type]][grepl("^KC", qcPassType[[type]]) |
                                    grepl("^MBON", qcPassType[[type]]) |
                                    grepl("^PAM", qcPassType[[type]])]

         t.ccx <- qcPassType[[type]][grepl("^PB", qcPassType[[type]])]

         print(paste0("-> muscle: n=", length(t.muscle)))
         print(t.muscle)
         print(paste0("-> mb: n=", length(t.mb)))
         print(t.mb)
         print(paste0("-> ccx: n=", length(t.ccx)))
         print(t.ccx)
         print(paste0("-> rest: n=", length(setdiff(qcPassType[[type]],
                                                    c(t.muscle, t.mb, t.ccx)))))
         print(setdiff(qcPassType[[type]], c(t.muscle, t.mb, t.ccx)))
      }
   }

   print("\n\nSUMMARY OF SUB-OPTIMAL SAMPLES-------------- ")

   qcFail <- dat$specs$sampleInfo[
               dat$specs$sampleInfo$driverType != "control" &
               ! dat$specs$sampleInfo$sample_name %in%
               dat$specs$sampleLists$good$allCriteria,]
   print(paste0("Sub-optimal samples n=", nrow(qcFail)))
   tx <- unique(qcFail[,c("driverType", "driver")])
   print("-> from # drivers:")
   print(table(tx$driverType))

   qcFailType <- list()
   print("-> from # cells:")
   for (type in c("pureCellType", "comboCellTypes", "broad", "dissected")) {
      qcFailType[[type]] <- unique(qcFail$celltype[qcFail$driverType == type])

      print(paste0("Number of qcFail ",type," n=",length(qcFailType[[type]])))
   }

   print("p(on) VAChT in Kdm2:")
   print(round(exp(dat$expr$rawpon$cells$p["VAChT","Kdm2"]),digits=2))


   print(paste0("avg ninaE in R1-6: ",
      round(dat$expr$rawpon$cells$raw["ninaE","R1-6"],digits=2),
      " TPM"))

   print("highest ninaE")
   print(sort(dat$expr$rawpon$cells$raw["ninaE",]))

   print("TPM(TfAP-2) in T4")
   print(dat$expr$rawpon$cells$raw["TfAP-2","T4"])
   print("TPM(TfAP-2) in T5")
   print(dat$expr$rawpon$cells$raw["TfAP-2","T5"])

   print("TPM(Gad1) in Mi9")
   print(dat$expr$rawpon$cells$raw["Gad1","Mi9"])

   print("TPM(Gad1) in Mi4")
   print(dat$expr$rawpon$cells$raw["Gad1","Mi4"])

   print("TPM(Gad1) in Gad1+ pure cells")
   t.gad1pos <- colnames(dat$expr$rawpon$cells$p)[
                           exp(dat$expr$rawpon$cells$p["Gad1",]) >= 0.8]
   t.gad1pos <- intersect(t.gad1pos,
                          dat$specs$sampleInfo$celltype[
                             dat$specs$sampleInfo$driverType == "pureCellType"])
   print(t.gad1pos)
   print(dat$expr$rawpon$cells$raw["Gad1",t.gad1pos])
   print("mean TPM(Gad1) in Gad1+ pure cells")
   print(mean(unlist(dat$expr$rawpon$cells$raw["Gad1",t.gad1pos])))

   print("TPM(AstC):"); print(dat$expr$rawpon$cells$raw["AstC",
                                                        c("C2","L4", "L5")])
   print("pon(AstC):"); print(exp(dat$expr$rawpon$cells$p["AstC",
                                                          c("C2","L4","L5")]))

   print("TPM(AstC-R1):"); print(dat$expr$rawpon$cells$raw["AstC-R1","T1"])
   print("pon(AstC-R1):"); print(exp(dat$expr$rawpon$cells$p["AstC-R1","T1"]))

   print("TPM(AstA):"); print(dat$expr$rawpon$cells$raw["AstA",c("Pm3","Tm2")])
   print("pon(AstA):"); print(exp(dat$expr$rawpon$cells$p["AstA",
                                                          c("Pm3","Tm2")]))

   print("TPM(AstA-R1):"); print(dat$expr$rawpon$cells$raw["AstA-R1",
                                    c("Mi1","Tm2", "Mi15", "Dm9")])
   print("pon(AstA-R1):"); print(exp(dat$expr$rawpon$cells$p["AstA-R1",
                                     c("Mi1","Tm2", "Mi15", "Dm9")]))


   print("TPM(fkh) range for expr (pon >= 0.9) cells:");
   print(range(dat$expr$rawpon$cells$raw["fkh",
                  exp(dat$expr$rawpon$cells$p["fkh",]) >= 0.9]))
   
   print("TPM(fkh):"); print(dat$expr$rawpon$cells$raw["fkh",c("Tm4","Dm9")])
   print("TPM(Ets65A):"); print(dat$expr$rawpon$cells$raw["Ets65A",
                                 c("Tm20","Glia_Eg","Dm3")])

   return(1)




   print("*********************************************************")
   print("R3.2. BENCHMARK NUMBERS")
   print("-----------------------")

   evaluatePcalls(dat, summaryOnly = TRUE)

   print(paste0("avg para in R1-6: ",
      round(dat$expr$rawpon$cells$raw["para","R1-6"],digits=2),
      " TPM"))

   print("highest para")
   print(sort(dat$expr$rawpon$cells$raw["para",]))


   for (gene in names(dat$specs$benchmarkExp)) {
      print("-------------------")
      print(gene)
      t.posCells <- dat$specs$benchmarkExp[[gene]]$pos
      t.negCells <- dat$specs$benchmarkExp[[gene]]$neg
      allCells <- c(t.posCells, t.negCells)
      print(paste0("ALL CELLS:",length(allCells)))
      print(paste0("POSITIVES:",t.posCells))
      print(paste0("NEGATIVES:",t.negCells))

      print("Missing columns:")
      print(setdiff(allCells, colnames(dat$expr$rawpon$cells$raw)))
      curRes <- data.frame(
         cell        = allCells,
         experiment  = c(rep("+",length(t.posCells)),
                         rep("-",length(t.negCells))),
         tpm         = unlist(dat$expr$rawpon$cells$raw[gene,allCells]),
         p           = unlist(exp(dat$expr$rawpon$cells$p[gene,allCells]))
      )
      curRes$modelCall <- "A"
      curRes$modelCall[curRes$p <= 0.2] <- "-"
      curRes$modelCall[curRes$p >= 0.8] <- "+"
      curRes$comparison <- "MISMATCH"

      curRes$comparison[curRes$modelCall == curRes$experiment] <- "match"

      print(curRes)
   }

   return(1)

}


makeFigs <- function( data,
                      figs = "all" ) {

   library(calibrate) #for textxy()
   library(dplyr)

   print("Making all figures")

   if (("all" %in% figs) | ("2" %in% figs)) {
      print("Figure 2 and associated tables")
      print("- Fig 2C")
      tx <- makeTAPINyieldFig(dat, figName="msFig2C")

      print("- Fig 2J,D")
      tx <- plotExprVsYield( data,
                       figName.expr_yield = "msFig2J",
                       figName.yield_yield = "msFig2D",
                       geneList = c("ninaE"),
                       returnData=FALSE, runFit=FALSE)


      print("- Fig 2E")
      tx <- plotDetectionVsYield(data, figName = "msFig2E")

      print("- Fig 2G")
      tx <- plotRepCorVsYield(   data, figName = "msFig2G")

      print("- Fig 2I")
      tx <- plotHmapPmat(dat, mode = "gc",
                          figName = "msFig2I.tan15_markers",
                          geneList = c("svp", "bab2", "erm", "ap",
                                       "pdm3", "pros", "sens"),
                          exprMode          = "raw",
                          tpmOverlay = TRUE,
                          rawNL = "fracmax",
                          plotTitle="",
                          legend=TRUE,
                          plotWidth         = 2,
                          plotHeight        = 2,
                          geneLabel         = TRUE,
                          cluster_cols      = FALSE,
                          cluster_rows      = FALSE,
                          cellGroups = list( c("L1","L2","L3","L4","L5",
                                               "R7","R8_Rh5","R8_Rh6"))
      )

      print("- Table S1")
      tx <- makeSamplesTable(dat, tabName = "msTableS1")

      print("- Data Tables")
      tx <- writeDataTables(dat, tabName="dataTable")




      print("- Fig 2F")
      tx <- plotTxMolVsYield(    data, figName = "msFig2F")

      print("- Fig 2H")
      tx <- plotScatterTPM(data,
                     mode = "all.genes",
                     figName = "msFig2H.T4T5_male_vs_female",
                     xSamples = c("T4.T5_d1_male_rep1",
                                  "T4.T5_d1_male_rep2"),
                     ySamples = c("T4.T5_d1_female_rep1",
                                  "T4.T5_d1_female_rep2"),
                     xlab = "T4T5 male expr (TPM+1)",
                     ylab = "T4T5 female expr\n(TPM+1)")

   }

   if (("all" %in% figs) | ("S2" %in% figs)) {
      print("Figure S2")

      tx <- plot.readCoverage(   data, figName = "msFigS2D", type="maxVsYield")

      tx <- plot.readCoverage(   data, figName = "msFigS2C", type="txPosition")
      tx <- plotReplicateCorrHistos(data, figName="msFigS2E")
   }

   if (("all" %in% figs) | ("3" %in% figs)) {

      print("Figure 3")

      tx <- plotExprBreadthHisto( data, figName.histo = "msFig3F",
                                        figName.genegroups = "msFig4C",
                                        figName.cdf   = "msFig3G" )

      print(tx$genegroups[c("DPR-INTERACTING PROTEINS","BEAT FAMILY",
                            "DEFECTIVE PROBOSCIS EXTENSION RESPONSE")])

      print("TOP 5 groups:",head(tx$genegroups,n=5))
      print("BOTTOM 5 groups:",tail(tx$genegroups,n=5))

      tx <- vizPfit(data,
                    geneList          = c("VAChT"),
                    figName.tpmHisto  = "msFig3A",
                    figName.tpm_pon   = "msFig3B",
                    labelCol          = c("black"),
                    labelCells        = c("Kdm2_d1") ,
                    labelDirectly     = TRUE,
                    labelLegend       = c("Kdm2") )

      tx <- vizPfit(data,
                    geneList          = c("ninaE"),
                    figName.tpmHisto  = "msFig3C",
                    figName.tpm_pon   = "msFig3D",
                    labelCol          = c("black"),
                    labelCells        = c("R1-6_d1") ,
                    labelDirectly     = TRUE,
                    labelLegend       = "R1-6" )

      tx <- makeCoHeatmap(data, figName="msFig3E")


      tx <- plotOnRangeHisto(            data, figName = "msFigS3B" )

      tx <- plotMuOnOffHistos( data,
                               figName.muOnHisto = "msFig3H",
                               figName.muOffHisto = "msFigS3C",
                               figName.dmuOnOffHisto = "msFig3I" )

      detailedBmGene <- list( fkh = "msFig4J", Ets65A = "msFig4L")

      for (gene in names(detailedBmGene)) {
         t.posCells <- data$specs$benchmarkExp[[gene]]$pos
         t.negCells <- data$specs$benchmarkExp[[gene]]$neg

         tx <- vizPfit(data, geneList = gene,
                       figName.tpmHisto = detailedBmGene[[gene]],
                       singlePlot = TRUE,
                       mode="cells",
                       pointCol = c(NA),
                       labelCells = c(t.posCells,t.negCells),
                       labelCol = c(rep("forestgreen",length(t.posCells)),
                                   rep("black",1,length(t.negCells))),
                       labelPch = c(rep(20,length(t.posCells)),
                                   rep(20,length(t.negCells))),
                       labelDirectly = TRUE)

      }

      tx <- analyzeRepPonConcordance(data, figName = "msFigS4A")

      bmRes <- evaluatePcalls(data, figName="msFigS3A")
      print(bmRes$mismatches)
      for (gene in names(bmRes$mismatches)) {
         vizPfit(data, geneList = gene,
                 figName.tpmHisto = paste0("msFigS3_benchmark_mismatch_", gene),
                 figName.tpm_pon =  paste0("msFigS3_benchmark_mismatch_", gene),
                 labelDirectly = TRUE,
                 singlePlot = TRUE,
                 labelCells = unique(dat$pFit$inDat$driver[dat$pFit$inDat$cell
                                       %in% bmRes$mismatches[[gene]]])  ,
                 labelCol = "red",
                 labelPtCex = 1.2,
                 pointCex = 0.3,
                 pointPch = 20,
                 labelLegend = "")
      }

   }


   if (("all" %in% figs) | ("4" %in% figs)) {
      print("Figure 4")
      tx <- plotHmapPmat.markers( data, figName.qc_pass = "msFig4B",
                                  tabName = "msTabS4",
                                  plotHeight=5, plotWidth=2.5,
                                  cellGroups.marker = list(
         pr = c("R1-6", "R7","R8_Rh5","R8_Rh6"),
         glia = c("Glia_Eg", "Glia_Psg", "Glia_Mg"),
         muscle = c("Muscles_App", "Muscles_Head"),
         pigment = c("Tm20", "Dm1")
                                                    ),
                                  cellGroups.heatmap= c(
                                    dat$specs$cellGroups,
                                    list(c( "Muscles_Head", "Muscles_App")))
                                  )


      dat$pFit.scoreTPM <- inferState.scoreTPM(dat,
         samples = dat$specs$sampleInfo$sample_name[
                      dat$specs$sampleInfo$driverType == "dissected"],
         geneList = sort(names(dat$pFit$geneEfit)))

      tx <- dat$pFit.scoreTPM
      paste0(rownames(tx$pon_gc)[apply(exp(tx$pon_gc),1,max) >
               apply(exp(dat$expr$rawpon$cells$p),1,max)],collapse=", ")


      tx <- plotScatterTPM(data,
                     mode = "coding",
                     figName = "msFig4E.T4_vs_T5",
                     celltype1 = "T4",
                     celltype2 = "T5",
                     labelGenes.right = c("TfAP-2"),
                     labelGenes.left = c(),
                     colorGenes = loadTFs.flybase(data),
                     colorGenes.legend = "TF",
                     xlab = "T4 expr (TPM+1)",
                     ylab = "T5 expr (TPM+1)")

      tx <- makeCoTree(data,
                       figName.cotree   = "extra_msFigS4X_cotree",
                       figName.allgenes = "msFig4A",
                       figName.TFgenes  = "extra_msFigS4X_tfgenetree")



      tx <- plotScatterTPM(data,
                     mode = "coding",
                     figName = "msFig4H.T5_d1_vs_T5_d2",
                     driver1 = "T5_d1",
                     driver2 = "T5_d2",
                     labelGenes.right = c("klg","bi","Con"),
                     labelGenes.left = c("dac","dysf","Dscam3"),
                     colorGenes = loadTFs.flybase(data),
                     colorGenes.legend = "TF",
                     xlab = "T5_d1 expr (TPM+1)",
                     ylab = "T5_d2 expr (TPM+1)")

   }

   if (("all" %in% figs) | ("S5" %in% figs)) {
      print("Figure S5")

      tx <- plotHmapPmat.cam(data, figName = "msFigS5A", mode="sparse" )
      tx <- plotHmapPmat.trueMarkers(data, figName = "msFigS4B")

      for (exprMode in c("raw", "p")) {
         tx <- plotHmapPmat(dat, mode = "gc",
                          figName = paste0("msFigS5B.beat_",exprMode),
                          flybaseGeneGroup = "beat",
                          exprMode          = exprMode,
                          plotTitle="",legend=FALSE,
                          plotWidth         = 3.5,
                          plotHeight        = 6.5,
                          geneLabel         = TRUE,
                          cluster_cols      = FALSE,
                          cluster_rows      = FALSE,
                          colType           = "gene",
                          cellGroups        = dat$specs$cellGroups
         )

         tx <- plotHmapPmat(dat, mode = "gc",
                          figName = paste0("msFigS5C.dip_",exprMode),
                          flybaseGeneGroup = "dip",
                          exprMode          = exprMode,
                          plotTitle="",legend=FALSE,
                          plotWidth         = 3.5,
                          plotHeight        = 6.5,
                          geneLabel         = TRUE,
                          cluster_cols      = FALSE,
                          cluster_rows      = FALSE,
                          colType           = "gene",
                          cellGroups        = dat$specs$cellGroups
         )

         tx <- plotHmapPmat(dat, mode = "gc",
                          figName = paste0("msFigS5D.dpr_",exprMode),
                          flybaseGeneGroup = "dpr",
                          exprMode          = exprMode,
                          plotTitle="",legend=FALSE,
                          plotWidth         = 3.5,
                          plotHeight        = 6.5,
                          geneLabel         = TRUE,
                          cluster_cols      = FALSE,
                          cluster_rows      = FALSE,
                          colType           = "gene",
                          cellGroups        = dat$specs$cellGroups
         )
      }

      tx <- plotExprOzkanPairs(dat, figName="msFigS5E")

   }

   if (("all" %in% figs) | ("5" %in% figs) | ("6" %in% figs)) {
      print("Figure 5,6")

      tx <- plotHmapPmat(dat, mode = "gc",
                          figName = "msFig5B.np",
                          geneList = c("AstA","AstA-R1","AstA-R2",
                                       "AstC","AstC-R1","AstC-R2",
                                       "Pdf","Pdfr"),
                          gaps_genes = c(3,6),
                          plotTitle = "", legend=FALSE,
                          plotWidth         = 2,
                          plotHeight        = 6.5,
                          geneLabel         = TRUE,
                          cluster_cols      = FALSE,
                          cluster_rows      = FALSE,
                          colType           = "gene",
                          cellGroups        = dat$specs$cellGroups
      )

      tx <- plotHmapPmat.nt_and_rec(     data,
                                         figName1 = "msFig5A",
                                         figName2 = "msFig6A")

      tx <- findReg.NT(                  data, figName = "msFigS6")

   }

   if (("all" %in% figs) | ("7" %in% figs)) {
      tx <- plotHmapPmat(dat, mode = "gc",
                          figName = "msFig7B",
                          geneList = c("ort", "HisCl1"),
                          tpmOverlay = TRUE,
                          plotTitle = "", legend=FALSE,
                          plotWidth         = 1,
                          plotHeight        = 2.25,
                          geneLabel         = TRUE,
                          cluster_cols      = FALSE,
                          cluster_rows      = FALSE,
                          colType           = "gene",
                          cellGroups        = list(c("R7","R8_Rh5","R8_Rh6",
                                                     "L1","Dm9", "Tm20",
                                                     "Mi1","Mi4","Mi15"))

      )

      tx <- plotSynapseReceptors(dat,
                                 presynaptic.name = "R8 111",
                                 genes = c("ort", "HisCl1"),
                                 figName="msFig7XA") #defaults to R8 111


      tx <- plotSynapseReceptors(dat,
                                 presynaptic.name = "R7 205",
                                 genes = c("ort", "HisCl1"),
                                 figName="msFig7XB")


      tx <- plotSynapseReceptors(dat,
                                 presynaptic.name = "C2 214",
                                 genes = c("Grd", "Rdl"),
                                 figName="msFig7XC")

      tx <- plotSynapseReceptors(dat,
                                 presynaptic.name = "C3 103",
                                 genes = c("Grd", "Rdl"),
                                 figName="msFig7XD")


      tx <- plotHmapPmat(dat, mode = "gc",
                          figName = "msFig7.GABA_A",
                          geneList = c("CG8916","Grd","Rdl","Lcch3"),
                          plotTitle = "", legend=FALSE,
                          plotWidth         = 1.5,
                          plotHeight        = 3,
                          geneLabel         = TRUE,
                          tpmOverlay = TRUE,
                          cluster_cols      = FALSE,
                          cluster_rows      = FALSE,
                          colType           = "gene",
                          cellGroups        = list(c("R1-6","R7",
                                                     "R8_Rh5", "R8_Rh6",
                                                "L1","L2","L3","L4","L5",
                                                "C2","C3","Lai","Lawf1","Lawf2",
                                                "T1",
                                                "Glia_Eg","Glia_Mg","Glia_Psg",
                                                "Dm8","Dm9","Mi1","Mi4","Mi9"))
      )

      tx <- plotHmapPmat(dat, mode = "gc",
                          figName = "msFig7H",
                          geneList = c("Ekar", "GluClalpha","CG3822","Eaat1"),
                          tpmOverlay = TRUE,
                          plotTitle = "", legend=FALSE,
                          plotWidth         = 1.5,
                          plotHeight        = 3,
                          geneLabel         = TRUE,
                          cluster_cols      = FALSE,
                          cluster_rows      = FALSE,
                          colType           = "gene",
                          cellGroups        = list(c("R1-6", "R7", "R8_Rh5",
                                                 "R8_Rh6",
                                                 "L1", "L2", "L3", "L4", "L5",
                                                 "C2","C3",
                                                 "Lawf1", "Lawf2", "Lai", "T1",
                                                 "Glia_Eg", "Glia_Psg",
                                                 "Glia_Mg",
                                                 "Dm8", "Dm9"))
      )


   }

   if (("all" %in% figs) | ("S7" %in% figs)) {
      tx <- plotHmapPmat(dat, mode = "gc",
                          figName = "msFigS7B",
                          geneList = c("Tk", "TkR86C","TkR99D","sNPF",
                                       "Nos", "shakB"),
                          tpmOverlay = TRUE,
                          plotTitle = "", legend=FALSE,
                          plotWidth         = 1.5,
                          plotHeight        = 2,
                          geneLabel         = TRUE,
                          cluster_cols      = FALSE,
                          cluster_rows      = FALSE,
                          colType           = "gene",
                          cellGroups        = list(c("C3", "L1", "L2", "Mi1",
                                                     "Mi4", "Tm1", "Tm4"))
      )
   }

   return(1)
}


getMsTPM <- function( dat, gene, celltype ) {

# Purpose: retrieve TPM for particular gene, celltype

   allSamples <- dat$specs$sampleInfo$sample_name[
                  dat$specs$sampleInfo$celltype %in% c(celltype)]

   curSamples <- intersect(paste0("tpm.",allSamples),
                           colnames(dat$pFit$inDat$O))

   print(paste0(curSamples,collapse=", "))
   print(round(mean(dat$pFit$inDat$O[gene,curSamples]),digits=2))

}


setSpecs <- function(){
# Purpose: centralized routine to specify parameters

   specs <- list()

# SETUP PARAMETERS

# - Paths

   specs$baseDir <- paste0("/groups/eddy/home/davisf/work/",
                           "/genomics/eddy_seq_data_analysis/opticlobe")

   specs$timeStamp<-format(Sys.time(), "%Y%m.%H%M")
   specs$outDir <- paste0(specs$baseDir,"/analysis/20180601.fig_edits")

   if (!file.exists(specs$outDir)) dir.create(specs$outDir, recursive = TRUE)


# LOAD GLOBALLY ACCESSIBLE INFORMATION

   specs$benchmarkLitFn <- paste0(specs$baseDir,
      "/data/benchmark/ol_expr_benchmark.txt")

   specs$benchmarkExp <- list(
      fkh = list(
         pos = c("LC4","LPLC1", "T4","T5", "Dm8", "Dm9", "Dm12", "Lai",
                  "LPLC2", "Dm11", "Dm1", "Dm4"),
         neg = c("Tm4","Pm4", "Dm10","L1","L2","L3","L4","L5",
                 "Glia_Eg","Glia_Mg", "Glia_Psg",
                 "T1", "LC10a", "Dm3", "C3", "LC6")),

      Ets65A = list(
         pos = c("Mi15", "T1", "C2", "C3", "Dm3"),
         neg = c("Mi1", "LPLC2", "Lai", "Tm9", "Glia_Eg", "Tm20")
      )
   )

   specs$sampleFn <- paste0(specs$baseDir,
      "/metadata/opticlobe_ms_rnaseq_samples.txt")

   specs$sampleGeneticsFn <- paste0(specs$baseDir,
      "/metadata/fred20180123c_tableS1B_all_drivers_in_progress.csv")

   specs$kallistoDir <- paste0(specs$baseDir,
      "/results/RNAseq/kallisto.BDGP6.91")

   specs$picardDir <- paste0(specs$baseDir,
      "/results/RNAseq/picard_stats.star.BDGP6.91")

# LOGIC

   specs$runModel <- list(
      numTasks          = 200,
      outDir            = paste0(specs$outDir,"/inferState_splits_20180130")
   )

   specs$inferState <- list(
      kfoldvalidation.k         = 10,
      minMaxTPM                 = 10,
      genes.testPanel           = c("Hdc","Gad1","VGlut","VAChT","Vmat",
                                    "alphaTub84B", "Act42A", "brp",
                                    "Tdc2","ninaE","Lsp2"),
      level1.bim.stanFn         = "~/ol/src/STAN/level1_ordered.stan",
      level1.bim.cv.stanFn      = "~/ol/src/STAN/level1_ordered_CV.stan",
      level1.uni.stanFn         = "~/ol/src/STAN/level1_unimodal.stan",
      level1.uni.cv.stanFn      = "~/ol/src/STAN/level1_unimodal_CV.stan"
   )


   specs$sampleInfo <- loadSamples(specs)

# - Define 'celltypes' that are combination, control, or single drivers
   specs$controlDrivers <- unique(c(
      specs$sampleInfo$celltype[grep("mocks", specs$sampleInfo$celltype)],
      specs$sampleInfo$celltype[grep("beads", specs$sampleInfo$celltype)],
      specs$sampleInfo$celltype[grep("input", specs$sampleInfo$celltype)]))

# Load file: Transcript info
   specs$transcriptInfoFn <- paste0(specs$baseDir,
      "/data/misc_files.BDGP6.91/BDGP6.91.ERCC.INTACT.transcript_info.txt")
   specs$transcriptInfo <- read.table( specs$transcriptInfoFn,
                                       quote     = "",
                                       header    = TRUE,
                                       sep       = "\t",
                                       as.is     = TRUE )

# Load file: ERCC concentrations
   specs$erccConcFn <- paste0(specs$baseDir,
                              "/data/external/ERCC/cms_095046.txt")

# Load benchmark expression patterns
   specs$bmExpr <- loadBenchmark(specs)

# Load files: GO annotation
   specs$goFn <- paste0(specs$baseDir,
      "/data/external/flybase/flybase.fbgn_goterm.2014_05.txt.gz")
   specs$goInfo <- loadGO(specs)

# Load files: InterPro annotation
   specs$interproFn <- paste0(specs$baseDir,
      "/data/external/flybase/fbgn_interpro.dmel-r6.02.20160607.txt.gz")
   specs$interproInfo <- loadInterpro(specs)

   specs$flyBaseTfsFn <- paste0(specs$baseDir,
      "/data/external/flybase/flybase_TFs_genegroup_FBgg0000745_20180221.txt")

   specs$fullFlyBaseGroupsFn <- paste0(specs$baseDir,
      "/data/external/flybase/gene_group_data_fb_2018_02.tsv.gz")

   specs$flyBaseGroupsFn <-  list(
      beat = paste0(specs$baseDir,
             "/data/external/flybase/beat_family_FBgg0000596.txt"),
      dip = paste0(specs$baseDir,
             "/data/external/flybase/dip_FBgg0000530.txt"),
      dpr = paste0(specs$baseDir,
             "/data/external/flybase/dpr_FBgg0000529.txt")
   )

# Load files: Ozkan13 PPI
   specs$ozkan13 <- list()
   specs$ozkan13$proteinInfoFn <- paste0(specs$baseDir,
      "/data/papers/ozkan13/ozkan_tableS1.txt")

   specs$ozkan13$proteinPairsFn <- paste0(specs$baseDir,
      "/data/papers/ozkan13/ozkan_tableS2.txt")

# Set file locations: Riveraalba11 connectome
   specs$riveraalba11 <- list()
   specs$riveraalba11$conMatrixFn <- paste0(specs$baseDir,
      "/data/papers/riveraalba11/riveraalba2011_tableS2.txt")

# Set file locations: Takemura13 connectome
   specs$takemura13 <- list()
   specs$takemura13$synapsesFn <- paste0(specs$baseDir,
      "/data/papers/takemura13/nature12450-s3.xls")

# Load files: Konstaninides18 expression data
   specs$konstantinides18 <- list()
   specs$konstantinides18$sampleInfoFn<- paste0(specs$baseDir,
      "/metadata/konstantinides18_rnaseq_samples.txt")

   specs$konstantinides18$markersFn <- paste0(specs$baseDir,
      "/data/papers/konstantinides18/",
      "konstantinides18_tableS1_cluster_markers.txt")

   specs$konstantinides18$scClusterAssFn  <- paste0(specs$baseDir,
      "/data/papers/konstantinides18/",
      "konstantinides18_scope_singlecell_cluster_assignment.csv")

   specs$konstantinides18$scClusterLabelsFn  <- paste0(specs$baseDir,
      "/data/papers/konstantinides18/",
      "konstantinides18_fig3a_cluster_labels.txt")

# Load files: manual neurotransmitter/neuropeptide classes
   specs$ntSystemsFn<-paste0(specs$baseDir,
      "/data/external/manual_nt_list/fly_neurotransmission_manuallist.txt")

   specs$ntSystems <- read.table(specs$ntSystemsFn,
                                 header = TRUE,
                                 sep    = "\t",
                                 as.is  = TRUE)

   specs$cellGroups <- setCellGroups()
   specs$cellGroups.all <- setCellGroups.all()

   return(specs)
}


setCellGroups.all <- function() {

   cellGroups.all <- list(

c( "R1-6",
   "R7_Rh3",
   "R7",
   "R8_Rh5",
   "R8_Rh6"),

c( "C2",
   "C2.C3",
   "C3",
   "L1",
   "L1.L2",
   "L2",
   "L3",
   "L4",
   "L5",
   "Lawf1",
   "Lawf2",
   "Lai",
   "Lat",
   "T1",
   "Glia_Eg",
   "Glia_Mg",
   "Glia_Psg",
   "lamina"),


c( "Dm1",
   "Dm10",
   "Dm11",
   "Dm12",
   "Dm3",
   "Dm4",
   "Dm8",
   "Dm9",
   "Mi1",
   "Mi15",
   "Mi4",
   "Mi9",
   "Pm3",
   "Pm4",
   "T4",
   "T4.T5",
   "T5",
   "Tm1",
   "Tm2",
   "Tm20",
   "Tm3",
   "Tm4",
   "Tm9",
   "Tm29",
   "TmY3",
   "TmY5a",
   "opticlobe"),

c( "LC10a",
   "LC10b",
   "LC10d",
   "LC10bc",
   "LC4",
   "LC6",
   "LC16",
   "LLPC1",
   "LPC1",
   "LPi-34",
   "LPLC1",
   "LPLC2",
   "LPTC_HS.VS"),

c( "KC_ab_c",
   "KC_ab_p",
   "KC_ab_s",
   "KC_ab_c.p.s",
   "KC_apbp_ap",
   "KC_apbp_m",
   "KC_gd",
   "MBON_bp1",
   "MBON_g1pedc",
   "PAM_1",
   "PAM_2",
   "PAM_3",
   "PAM_4",
   "PAM_5",
   "PAM_6",
   "PAM_7",
   "PAM_8"),

c( "PB_1",
   "PB_2",
   "PB_3",
   "PB_4",
   "PB_5",
   "PB_6"),

c( "ChAT",
   "Gad1",
   "VGlut"),

c( "CCAP",
   "Crz",
   "Ilp2",
   "Dsk",
   "Kdm2",
   "NPF",
   "Pdf",
   "lLNv"),


c( "Muscles_Head",
   "Muscles_App")
   )

   return(cellGroups.all)
}


setCellGroups <- function() {

   cellGroups = list(
              c( "R1-6", "R7", "R8_Rh5", "R8_Rh6") ,
              c(
      "C2", "C3",
      "L1", "L2", "L3", "L4", "L5",
      "Lawf1", "Lawf2",
      "Lai",
      "T1",
      
      "Glia_Eg",
      "Glia_Mg",
      "Glia_Psg"
      ),
      
      c(
      "Dm1",
      "Dm3",
      "Dm4",
      "Dm8",
      "Dm9",
      "Dm10",
      "Dm11",
      "Dm12",

      "Mi1",
      "Mi4",
      "Mi9",
      "Mi15",

      "Pm3",
      "Pm4",
      
      "T4",
      "T5",
      "Tm1",
      "Tm2",
      "Tm3",
      "Tm4",
      "Tm9",
      "Tm20",
      "Tm29",
      "TmY3",
      "TmY5a"
      ),
      
      c(
      "LC4",
      "LC6",
      "LC10a",
      "LC10b",
      "LC16",
      "LLPC1",
      "LPC1",
      "LPLC1",
      "LPLC2",
      "LPi-34"
      ),
      
      c(
      "KC_ab_c",
      "KC_ab_p",
      "KC_ab_s",
      "KC_gd",
      "PAM_1",
      "PAM_3",
      "PAM_4"
      ),
      
      c(
   "PB_1",
   "PB_2",
   "PB_3",
   "PB_4",
   "PB_5"
      )

   )

   return(cellGroups)
}


loadSamples <- function( specs ){
# Purpose: Load information about optic lobe RNA-seq samples

   print("Load sample information")
   sampleInfo <- read.table( specs$sampleFn,
                             header       = TRUE,
                             sep          = "\t",
                             as.is        = TRUE,
                             stringsAsFactors=FALSE,
                             comment.char = '#')

   sampleInfo$yield.cdna.ng_uL[sampleInfo$yield.cdna.ng_uL == "unk"] <- NA
   sampleInfo$yield.cdna.ng_uL <- as.numeric(sampleInfo$yield.cdna.ng_uL)

   sampleInfo$yield.nuclei.thousands[sampleInfo$yield.nuclei.thousands == "unk"] <- NA
   sampleInfo$yield.nuclei.thousands <- as.numeric(sampleInfo$yield.nuclei.thousands)

# Add easy sample names:
   sampleInfo$sample_name <- paste0(sampleInfo$driver, "_rep",
          sapply(1:length(sampleInfo$driver),
                 function(x) sum(sampleInfo$driver[1:x] ==
                                 sampleInfo$driver[x])))


   print("-> Loading sample genetic info")
   geneticInfo <- read.csv( specs$sampleGeneticsFn )
   geneticInfo[,c("X","X.1","X.2","X.3","X.4")] <- NULL

   sampleInfo <- merge(sampleInfo, geneticInfo,
                       by.x="driver", by.y="driverID", all.x=TRUE)
                   
   return(sampleInfo)
}



loadExpr <- function( specs,
                      sampleInfoFn, #if want to pass custom sampleInfo table
                      printTxExprMat    = FALSE,
                      printTxCleanMat   = FALSE,
                      printGeneExprMat  = FALSE,
                      printGeneCleanMat = FALSE,
                      loadExons         = FALSE ) {

# Purpose: load KALLISTO tables of transcript abundance,
#          generates merged tables of transcript and gene estimates,
#          generates tables of exon coverage counts
#          plots rRNA content vs yield

   library(dplyr)

# Load BDGP transcript information

   print("Loading transcript info")
   transcriptInfo <- specs$transcriptInfo

   if (missing(sampleInfoFn)) {
      sampleInfo <- specs$sampleInfo
   } else {
      print(paste0("Reading in custom sample information: ", sampleInfoFn))
      sampleInfo <- read.table(sampleInfoFn,
                               header       = TRUE,
                               sep          = "\t",
                               as.is        = TRUE,
                               stringsAsFactors=FALSE,
                               comment.char = '#')
   }



# Build gene_id <-> gene_name map
   geneInfo <- unique(transcriptInfo[,c("gene_id","gene_name")])

   nSamples <- nrow(sampleInfo)
   print(paste("Will process ",nSamples, "samples"))

   print("Loading KALLISTO estimates")

   txExprMat <- NULL
   for (i in 1:nSamples){
      print(paste("Reading expression from sample # ", i,
                  ": ", sampleInfo$sample_name[i], " ", date()))

#target_id       length  eff_length      est_counts      tpm
#FBtr0082757     1844    1595    9       0.718348

      abundFn <- paste0( specs$kallistoDir,"/",
                         sampleInfo$sampleID[i],
                         "/abundance.tsv")

      if (! file.exists(abundFn)) {
         print(paste0("-> skipping sample, abundance file not found: ",
                      abundFn))
         next;
      }

      curAbund <- read.table(abundFn,
                             header     = TRUE,
                             sep        = "\t",
                             colClasses = c("character", "NULL",
                                            "NULL", "numeric", "numeric"))

      colnames(curAbund) <- c( "transcript_id",
                               paste0("est_counts.",
                                      sampleInfo$sample_name[i]),
                               paste0("tpm.", sampleInfo$sample_name[i]))

      if (is.null(txExprMat)) {
         txExprMat <- merge( transcriptInfo,
                             curAbund,
                             all.y      = TRUE,
                             by         = "transcript_id" )
      } else {
         txExprMat$length <- NULL
         txExprMat <- merge( txExprMat,
                             curAbund,
                             by = c("transcript_id"))
      }

#      2L      7528    8116    FBgn0031208     11
#      2L      8192    9484    FBgn0031208     5

      print(paste("    -> DONE: ",date()))
   }


## Prepare transcript expression matrices for output

   allCols <- colnames(txExprMat)
   tpmCols <- allCols[grep("tpm.",allCols)]
   estCountCols <- allCols[grep("est_counts.",allCols)]

   # Reorder transcript matrix columns
   newTxExprMat <- txExprMat[, c("transcript_id", "transcript_name",
                               "gene_id", "gene_name", "gene_biotype",
                               tpmCols, estCountCols)]
   # Format to max 3 decimal digits
   if (printTxExprMat) {
      x <- newTxExprMat
      x$length <- NULL
      for (i in c(tpmCols, estCountCols)) {
         x[,i] <- sprintf("%.3f", x[,i])}

      gz1 <- gzfile(paste0(specs$outDir,
                           "/transcriptExpressionMatrix.txt.gz"),"w")
      write.table(x,
                  file          = gz1,
                  col.names     = TRUE,
                  row.names     = FALSE,
                  quote         = FALSE,
                  sep           = "\t")
      close(gz1)
      x <- NULL
   }


##   biotype content
   biotypeExpr <- aggregate( txExprMat[,c(tpmCols)],
      by=list(gene_biotype = txExprMat[["gene_biotype"]]),
      FUN = sum )
   rownames(biotypeExpr)<-biotypeExpr$gene_biotype
   biotypeExpr$gene_biotype<-NULL

# Keep only protein-coding genes

   cleanTxExprMat <- newTxExprMat[newTxExprMat$gene_biotype %in% c("protein_coding"),]

   print("Renormalizing expression table after keeping just protein coding genes")
   for (i in (1:length(tpmCols))) {
      cleanTxExprMat[, tpmCols[i]] <- (1E6 * cleanTxExprMat[, tpmCols[i]] / 
                                      sum(cleanTxExprMat[, tpmCols[i]]))
   }
   print("-> DONE")

   # Format to max 3 decimal digits
   if (printTxCleanMat) {
      x <- cleanTxExprMat
      x$length <- NULL
      for (i in c(tpmCols, estCountCols)) x[, i] <- sprintf("%.3f", x[, i])

      gz1 <- gzfile(paste0(specs$outDir,
                "/transcriptExpressionMatrix.no_RRNA_ERCC_INTACT.txt.gz"), "w")
      write.table(x,
                  file          = gz1,
                  col.names     = TRUE,
                  row.names     = FALSE,
                  quote         = FALSE,
                  sep           = "\t")
      close(gz1)
      x <- NULL
   }

   print("Calculating gene expression matrix")

   geneExprMat <- aggregate( txExprMat[, c(tpmCols, estCountCols)],
                             by = list(gene_id    = txExprMat[["gene_id"]],
                                       gene_name  = txExprMat[["gene_name"]]),
                             FUN = sum )

   # Format to max 3 decimal digits
   if (printGeneExprMat) {
      x <- geneExprMat
      for (i in c(tpmCols, estCountCols)) x[,i] <- sprintf("%.3f", x[,i])

      gz1 <- gzfile(paste0(specs$outDir, "/geneExpressionMatrix.txt.gz"), "w")
      write.table(x,
                  file          = gz1,
                  col.names     = TRUE,
                  row.names     = FALSE,
                  quote         = FALSE,
                  sep           = "\t")
      close(gz1)
      x <- NULL
   }


   cleanGeneExprMat <- aggregate( cleanTxExprMat[,c(tpmCols, estCountCols)],
                           by = list(gene_id   = cleanTxExprMat[["gene_id"]],
                                     gene_name = cleanTxExprMat[["gene_name"]]),
                           FUN = sum )

   # Format to max 3 decimal digits
   if (printGeneCleanMat) {
      x <- cleanGeneExprMat
      for (i in c(tpmCols, estCountCols)) {
         x[,i] <- sprintf("%.3f", x[,i]) }

      gz1 <- gzfile(paste0(specs$outDir,
                     "/geneExpressionMatrix.no_RRNA_ERCC_INTACT.txt.gz"),"w")
      write.table(x,
                  file          = gz1,
                  col.names     = TRUE,
                  row.names     = FALSE,
                  quote         = FALSE,
                  sep           = "\t")
      close(gz1)
      x <- NULL
   }


   exprMat.forDE <- cleanGeneExprMat[, c("gene_id", estCountCols)]
   rownames(exprMat.forDE) <- exprMat.forDE$gene_id
   exprMat.forDE$gene_id <- NULL


# Set up bycell gene expression matrix
   print("Setup gene expression table by cell type")
   cleanGeneExprMat.bycell <- cleanGeneExprMat[,c("gene_id","gene_name")]
   rownames(cleanGeneExprMat.bycell) <- cleanGeneExprMat.bycell$gene_name

   for (celltype in sort(sampleInfo$celltype)){

         curSamples <- sampleInfo$sample_name[
                           sampleInfo$celltype == celltype]

         curTpmCols <- paste0("tpm.", curSamples)
         newTpmCol <- paste0("tpm.", celltype)

         if (length(curTpmCols) == 1) {
            curAvgTPM <- cleanGeneExprMat[, curTpmCols[1]]
         } else {
            curAvgTPM <- apply(cleanGeneExprMat[, curTpmCols], 1, mean)
         }
         cleanGeneExprMat.bycell[, newTpmCol] <- curAvgTPM

   }


   return(list(
      biotypeExpr       = biotypeExpr,
      geneExpr          = cleanGeneExprMat,
      geneExpr.withERCC = geneExprMat,
      geneExpr.bycell   = cleanGeneExprMat.bycell,
      txExpr            = cleanTxExprMat,
      txExpr.withERCC   = txExprMat,
      exprCounts        = exprMat.forDE,
      transcriptInfo    = transcriptInfo,
      geneInfo          = geneInfo,
      sampleInfo        = sampleInfo
   ))

}


loadGO <- function( specs ){
# Purpose: laod GO terms

   if (missing(specs)) specs <- setSpecs()

   goInfo <- read.table(gzfile(specs$goFn),
                        as.is   = TRUE,
                        sep     = "\t",
                        header  = TRUE,
                        comment = "",
                        quote   = "")
   return(goInfo)
}



loadInterpro <- function( specs ) {

   if (missing(specs)) specs <- setSpecs()

   interproInfo <- read.table(gzfile(specs$interproFn),
                        as.is   = TRUE,
                        sep     = "\t",
                        header  = TRUE,
                        comment = "",
                        quote   = "")
   return(interproInfo)
}


loadFlyBaseGeneGroups <- function( specs,
                                   fn,
                                   terminalGroupsOnly = FALSE) {

# Purpose: Load FlyBase Gene Groups

   if (missing(specs)) specs <- setSpecs()

   if (missing(fn)) fn <- specs$fullFlyBaseGroupsFn

   fbgg <- read.table(gzfile(fn),
                        as.is   = TRUE,
                        sep     = "\t",
                        header  = FALSE,
                        quote   = "")
   colnames(fbgg) <-  c(
      "FB_group_id", "FB_group_symbol", "FB_group_name", "Parent_FB_group_id",
      "Parent_FB_group_symbol", "Group_member_FB_gene_id",
      "Group_member_FB_gene_symbol")

   if (!terminalGroupsOnly) {

      fbggGroupInfo <- unique(fbgg[,c("FB_group_id",
                                      "FB_group_name",
                                      "FB_group_symbol")])
      group2parent <- c()
      for (i in 1:nrow(fbgg)) {
         if (fbgg$Parent_FB_group_symbol[i] == "") {next;}
         group2parent[fbgg$FB_group_id[i]] <- fbgg$Parent_FB_group_id[i]
      }
   
      group2allparents <- list()
      for (fbgroup in names(group2parent)) {
         group2allparents[[fbgroup]] <- c()
         curNode <- fbgroup
   
         parRows <- fbgg[fbgg$FB_group_id %in% curNode &
                         fbgg$Parent_FB_group_id != "",]
   
         while (nrow(parRows) > 0) {
   
            curParent <- unique(parRows$Parent_FB_group_id)
            curNode <- curParent
            group2allparents[[fbgroup]] <- c(group2allparents[[fbgroup]], curNode)
   
            parRows <- fbgg[fbgg$FB_group_id %in% curNode &
                           fbgg$Parent_FB_group_id != "",]
            if (nrow(parRows) == 0) {break;}
   
         }
         group2allparents[[fbgroup]] <- unique(group2allparents[[fbgroup]])
      }
   
      groupid2parent <- data.frame(childGroupID = unlist(
                                    lapply(names(group2allparents),
                                           function(x,tx) { rep(x, length(tx[[x]]))},
                                           tx=group2allparents)),
                                   parentGroupID = unlist(group2allparents),
                                   stringsAsFactors = FALSE)
      rownames(groupid2parent) <- NULL

# expand the terminal group -- gene assignments to all parent groups

# using: FB_group_id, FB_group_name, Group_member_FB_gene_id
      geneRows <- fbgg[fbgg$Group_member_FB_gene_id != "",
                       c("Group_member_FB_gene_id",
                         "Group_member_FB_gene_symbol",
                         "FB_group_id")]

      expRows <-  merge(geneRows, groupid2parent,
                        by.x="FB_group_id",
                        by.y="childGroupID")

      expRows <- na.omit(expRows)

      expRows$FB_group_id <- expRows$parentGroupID
      expRows$parentGroupID <- NULL

      mergeRows <- rbind(geneRows, expRows)

      mergeRows<- merge(mergeRows, fbggGroupInfo, by = "FB_group_id")
      mergeRows <- unique(mergeRows)

      fbgg <- mergeRows

   }

   return(fbgg)
}


loadTFs.flybase <- function( data ){

# Purpose: Load TF list from FlyBase gene group 

   tx <- read.table(data$specs$flyBaseTfsFn,
                    header = TRUE,
                    sep="\t",
                    as.is  = TRUE)

   return(intersect(tx$gene_name, data$expr$geneExpr$gene_name))

   return(tx$gene_name)
}



loadOzkanPPI <- function( specs ) {
# Purpose: load Ozkan et al., 2013 extracellular protein interactions

   outdat <- list()

   proteinInfo <- read.table(specs$ozkan13$proteinInfoFn,
                             header       = TRUE,
                             sep          = "\t",
                             comment.char = "$",
                             as.is        = TRUE)

   proteinPairs <- read.table(specs$ozkan13$proteinPairsFn,
                              header       = TRUE,
                              sep          = "\t",
                              comment.char = "#",
                              as.is        = TRUE)

   outdat$proteinInfo <- proteinInfo
   outdat$proteinPairs <- proteinPairs

# get pairs in CG codes.
   pairs <- proteinPairs[,c("Interactor.1", "Interactor.2")]
   names(pairs) <- c("symbol1", "symbol2")
   translateTable <- proteinInfo[,c("Symbol","Annotation","Flybase.Gene.ID")]
   names(translateTable) <- c("symbol","geneID","FBgn")

   pairs <- merge(pairs, translateTable, by.x="symbol1", by.y="symbol",all.x=TRUE)
   pairs$geneID1 <- pairs$geneID
   pairs$geneID <- NULL

   pairs$FBgn1 <- pairs$FBgn
   pairs$FBgn <- NULL


   pairs <- merge(pairs, translateTable, by.x="symbol2", by.y="symbol",all.x=TRUE)
   pairs$geneID2 <- pairs$geneID
   pairs$geneID <- NULL
   pairs$FBgn2 <- pairs$FBgn
   pairs$FBgn <- NULL

# Update gene symbols using FBgn's and current dat$specs$transcriptInfo

   pairs$symbol1.orig <- pairs$symbol1
   pairs$symbol2.orig <- pairs$symbol2

   pairs$symbol1 <- dat$specs$transcriptInfo$gene_name[
      match(pairs$FBgn1, dat$specs$transcriptInfo$gene_id)]

   pairs$symbol2 <- dat$specs$transcriptInfo$gene_name[
      match(pairs$FBgn2, dat$specs$transcriptInfo$gene_id)]

   outdat$genePairs <- unique(pairs)
   return(outdat)

}


compareReplicates <- function( data,
                               plotRepScatter          = TRUE,
                               returnDatOnly           = FALSE,
                               mode                    = "bioReps",
                               selectFrom ){
# Purpose: make replicate scatterplots and quantify correlation

   allSamples <- data$specs$sampleInfo[data$specs$sampleInfo$driverType != "control",]

   if (!missing(selectFrom)) {
      allSamples <- allSamples[allSamples$sample_name %in% selectFrom,] }

   if ( mode == "bioReps" ) {
      print("Computing biological replicate correlations")

      bioRepCorMat.sample1<-c()
      bioRepCorMat.sample2<-c()
      bioRepCorMat.cor<-c()
      for (driver in sort(unique(allSamples$driver))){
         curSamples <- allSamples[allSamples$driver == driver, ]
         if (nrow(curSamples) == 1) next;
         print(paste0("driver ",driver," n=",nrow(curSamples)," replicates"))

         curRepR <- c()
         for (i in (1:(nrow(curSamples) - 1))) {

            colX <- paste0("tpm.", curSamples$sample_name[i])
            curX <- data$expr$geneExpr[, colX]

            for (j in ((i + 1):nrow(curSamples))) {

               if (curSamples$intact.protocol[i] !=
                   curSamples$intact.protocol[j]) next;

               colY <- paste0("tpm.", curSamples$sample_name[j])
               curY <- data$expr$geneExpr[, colY]

               curR <- cor(log2(1 + curX), log2(1 + curY))
               curRepR <- c(curRepR, curR)

               bioRepCorMat.sample1 <- c(bioRepCorMat.sample1,
                                         curSamples$sample_name[i])
               bioRepCorMat.sample2 <- c(bioRepCorMat.sample2,
                                         curSamples$sample_name[j])
               bioRepCorMat.cor <- c(bioRepCorMat.cor, curR)

               if (plotRepScatter) {
                  curX <- data$expr$geneExpr[, colX]
                  curY <- data$expr$geneExpr[, colY]

                  pdf(paste0(data$specs$outDir, "/scatter_bio_replicate_",
                             driver, "_",
                             curSamples$sample_name[i], "_vs_",
                             curSamples$sample_name[j], ".pdf"))

                  par(cex=1.3)

                  plot(1 + curX,
                       1 + curY,
                       log = "xy",
                       main= "",
                       xlab= paste0(curSamples$sample_name[i], " (TPM + 1)"),
                       ylab= paste0(curSamples$sample_name[j], " (TPM + 1)"),
                       pch = 20,
                       cex = 0.5)

                  legend("bottomright",
                         legend= paste0("pearson r=", sprintf("%.2f", curR)),
                         bty   = "n")

                  dev.off()
               }
            }
         }
      }
      range.repR <- range(bioRepCorMat.cor)
      print(paste("Biological Replicate r range: ",
                  paste0(range.repR, collapse=" - ")))

      bioRepCorMat<-data.frame( rep1 = bioRepCorMat.sample1,
                                rep2 = bioRepCorMat.sample2,
                                cor  = bioRepCorMat.cor,
                                stringsAsFactors = FALSE)
      if (returnDatOnly) {return(bioRepCorMat)}

      outFn<-paste0(data$specs$outDir,
                       "/table_biological_replicate_correlation.txt")
      write.table(bioRepCorMat,
                  outFn,
                  col.names = TRUE,
                  row.names = FALSE,
                  quote     = FALSE,
                  sep       = "\t")

      pdf(paste0(data$specs$outDir,
                 "/hist_biological_replicate_correlations.pdf"))
      par(cex=1.5, mar=c(5, 5, 1, 1))
      hist(bioRepCorMat.cor,
           lwd  = 2,
           xlab = "Biological replicate correlation (pearson)",
           main = "")
      dev.off()
   }

## ALTERNATE DRIVERS
   if ( mode == "drivers" ) {

      driverRepCorMat.sample1<-c()
      driverRepCorMat.sample2<-c()
      driverRepCorMat.cor<-c()

      print(paste0("Computing alternate driver correlation "))
      for (celltype in c(unique(allSamples$celltype))){

         curSamples <- allSamples[allSamples$celltype == celltype, ]

         curRepR <- c()
         if (length(unique(allSamples$driver[allSamples$celltype == celltype])) == 1) next;

         if (nrow(curSamples) == 1) next;

         nSamplePairs <- 0
         for (i in (1:(nrow(curSamples) - 1))) {
            colX <- paste0("tpm.", curSamples$sample_name[i])
            curX <- data$expr$geneExpr[, colX]

            for (j in ((i + 1):nrow(curSamples))) {

               if (curSamples$intact.protocol[i] != curSamples$intact.protocol[j] |
                   curSamples$driver[i] == curSamples$driver[j]) next;

               nSamplePairs <- nSamplePairs + 1
               colY <- paste0("tpm.", curSamples$sample_name[j])
               curY <- data$expr$geneExpr[, colY]

               print(paste("comparing ", colX, " to ", colY))

               curR <- cor(log2(1 + curX), log2(1 + curY))
               curRepR <- c(curRepR, curR)

               driverRepCorMat.sample1 <- c(driverRepCorMat.sample1,
                                          curSamples$sample_name[i])
               driverRepCorMat.sample2 <- c(driverRepCorMat.sample2,
                                          curSamples$sample_name[j])
               driverRepCorMat.cor <- c(driverRepCorMat.cor, curR)

               if (plotRepScatter) {
                  curX <- data$expr$geneExpr[, colX]
                  curY <- data$expr$geneExpr[, colY]

                  pdf(paste0(data$specs$outDir, "/scatter_driver_replicate_",
                             celltype, "_",
                             curSamples$sample_name[i], "_vs_",
                             curSamples$sample_name[j], ".pdf"))

                  par(cex=1.3)

                  plot(1 + curX,
                       1 + curY,
                       log = "xy",
                       main= "",
                       xlab= paste0(curSamples$sample_name[i], " (TPM + 1)"),
                       ylab= paste0(curSamples$sample_name[j], " (TPM + 1)"),
                       pch = 20,
                       cex = 0.5)

                  legend("bottomright",
                         legend= paste0("pearson r=", sprintf("%.2f", curR)),
                         bty   = "n")

                  dev.off()
               }
            }
         }
         print(paste0(" celltype ",celltype," n=",nSamplePairs," sample pairs"))
      }
      range.repR <- range(driverRepCorMat.cor)
      print(paste("Driver replicate r range: ",
                  paste0(range.repR, collapse=" - ")))

      driverRepCorMat<-data.frame( rep1 = driverRepCorMat.sample1,
                                 rep2 = driverRepCorMat.sample2,
                                 cor  = driverRepCorMat.cor,
                                 stringsAsFactors = FALSE)
      if (returnDatOnly) {return(driverRepCorMat)}
      outFn<-paste0(data$specs$outDir,
                       "/table_driver_replicate_correlation.txt")
      write.table(driverRepCorMat,
                  outFn,
                  col.names = TRUE,
                  row.names = FALSE,
                  quote     = FALSE,
                  sep       = "\t")

      pdf(paste0(data$specs$outDir,
                 "/hist_driver_replicate_correlations.pdf"))
      par(cex=1.5, mar=c(5, 5, 1, 1))
      hist(driverRepCorMat.cor,
           lwd=2,
           xlab="Driver replicate correlation (pearson)",
           main="")
      dev.off()
   }


## DIFFERENT PROTOCOLS: 2013 v 2015 for same driver
   if ( mode == "protocols" ) {
      allSamples <- data$specs$sampleInfo
      allSamples <- allSamples[
         grep("no_driver", allSamples$driver, invert=TRUE), ]
      allSamples <- allSamples[
         grep("aFLAG", allSamples$driver, invert=TRUE), ]

      protRepCorMat<-NULL
      protRepCorMat.sample1<-c()
      protRepCorMat.sample2<-c()
      protRepCorMat.cor<-c()

      print(paste0("Computing alternate protocol correlations "))

      for (driver in c(unique(allSamples$driver))){

         curSamples <- allSamples[allSamples$driver == driver, ]
         if (length(unique(curSamples$intact.protocol)) == 1) next;

         curRepR <- c()
         nSamplePairs <- 0
         for (i in (1:(nrow(curSamples) - 1))) {

            colX <- paste0("tpm.", curSamples$sample_name[i])
            curX <- data$expr$geneExpr[, colX]

            for (j in ((i + 1):nrow(curSamples))) {

               if (curSamples$intact.protocol[i] ==
                   curSamples$intact.protocol[j]) next;

               nSamplePairs <- nSamplePairs + 1

               colY <- paste0("tpm.", curSamples$sample_name[j])
               curY <- data$expr$geneExpr[, colY]

               curR <- cor(log2(1 + curX), log2(1 + curY))
               curRepR <- c(curRepR, curR)

               protRepCorMat.sample1 <- c(protRepCorMat.sample1,
                                          curSamples$sample_name[i])
               protRepCorMat.sample2 <- c(protRepCorMat.sample2,
                                          curSamples$sample_name[j])
               protRepCorMat.cor <- c(protRepCorMat.cor, curR)

               if (plotRepScatter) {
                  curOut<-paste0(data$specs$outDir,
                                 "/scatter_protocol_replicate_",
                                 driver, "_",
                                 curSamples$sample_name[i], "_vs_",
                                 curSamples$sample_name[j], ".pdf" )
                  pdf(curOut)
                  plot(1 + curX,
                       1 + curY,
                       log   = "xy",
                       main  = "",
                       xlab  = paste0(curSamples$sample_name[i], " (TPM + 1)"),
                       ylab  = paste0(curSamples$sample_name[j], " (TPM + 1)"),
                       pch   = 20,
                       cex   = 0.5)

                  legend("bottomright",
                         legend= paste0("pearson r=", sprintf("%.2f", curR)),
                         bty   = "n")

                  contamMat <- data$expr$geneExpr[
                     data$expr$geneExpr$gene_name %in% c(
                        "ninaE", "Arr2","Gad1", "VAChT", "VGlut"),
                     c(colX,colY,"gene_name")]
                  library(calibrate)
                  textxy(contamMat[,colX],
                         contamMat[,colY],
                         contamMat$gene_name,
                         col    = "red",
                         cex    = 1.5)

                  points(1 + contamMat[, colX],
                         1 + contamMat[, colY],
                         pch=20, cex=2, col="red")

                  dev.off()
               }
            }
         }
         print(paste0("   ", sort(curRepR), collapse=", "))
         print(paste0(" driver ",driver," n=",nSamplePairs," sample pairs"))
      }
      range.repR <- range(protRepCorMat.cor)
      print(paste("Protocol Replicate r range: ", paste0(range.repR, collapse=" - ")))

      protRepCorMat<-data.frame( rep1 = protRepCorMat.sample1,
                                 rep2 = protRepCorMat.sample2,
                                 cor  = protRepCorMat.cor,
                                 stringsAsFactors = FALSE)
      if (returnDatOnly) {return(protRepCorMat)}
      outFn<-paste0(data$specs$outDir,
                       "/table_protocol_replicate_correlation.txt")
      write.table(protRepCorMat,
                  outFn,
                  col.names = TRUE,
                  row.names = FALSE,
                  quote     = FALSE,
                  sep       = "\t")

      pdf(paste0(data$specs$outDir,
                 "/hist_protocol_replicate_correlations.pdf"))
      par(cex=1.5, mar=c(5, 5, 1, 1))
      hist(protRepCorMat.cor,
           lwd=2,
           xlab="Protocol replicate correlation (pearson)",
           main="")
      dev.off()
   }

   return(TRUE)
}


makeSampleLists <- function( data ) {
# purpose: make lists of samples to use for most analysis and inferState() analysis

   sampleList <- list()

   sampleList$control <- data$specs$sampleInfo$sample_name[ data$specs$sampleInfo$driverType == "control"]
   sampleList$all <- data$specs$sampleInfo$sample_name[data$specs$sampleInfo$driverType != "control"]
   sampleList$good <- selectGoodSamples(data, selectFrom = sampleList$all)

   return(sampleList)

}


selectGoodSamples <- function( data, selectFrom ) {
# Purpose: select samples based on several quality criteria

   cutoffs <- list(
      repCor = 0.85,
      nDetGenes = 8500,
      yield.cdna.ng_uL = 100 
   ) ;

   qcList <- list()
   criteria <- list()
   { # Only keep samples with max(replicate correlation) >= cutoffs$repCor

      library(magrittr); library(dplyr);
      tx1 <- compareReplicates(data, mode="bioReps", plot=FALSE, returnDatOnly = TRUE)
      qcList$repCor <- tx1
      criteria$repCor <- unique(unlist(tx1[tx1$cor >= cutoffs$repCor,
                                           c("rep1", "rep2")]))
   }

   {# genes detected >= cutoffs$nDetGenes
      numDet <- vector()
      tpmCols <- colnames(data$expr$geneExpr)
      tpmCols <- tpmCols[grep("tpm.",tpmCols)]

      for (tpmCol in tpmCols) {
          numDet[gsub("tpm.","",tpmCol)] <-
            sum(data$expr$geneExpr[,tpmCol] > 0) }

      qcList$numDet <- numDet

      criteria$nDetGenes <- names(numDet[numDet >= cutoffs$nDetGenes])
   }

   {# cDNA yield >= cutoffs$yield.cdna.ng_uL
      criteria$yield.cdna.ng_uL <- data$specs$sampleInfo$sample_name[
         !is.na(data$specs$sampleInfo$yield.cdna.ng_uL) &
         data$specs$sampleInfo$yield.cdna.ng_uL != "unk" &
         data$specs$sampleInfo$yield.cdna.ng_uL >= cutoffs$yield.cdna.ng_uL]
   }

   if (!missing(selectFrom)) {
      for (type in names(criteria)) {
         criteria[[type]] <- intersect(criteria[[type]], selectFrom)
      }
   }

   criteria$allCriteria <- intersect(criteria$repCor, criteria$nDetGenes)
   criteria$allCriteria <- intersect(criteria$allCriteria, criteria$yield.cdna.ng_uL)

# only keeps samples if more than 1 replicate

   curCells <- data$specs$sampleInfo$driver[
         data$specs$sampleInfo$sample_name %in% criteria$allCriteria]

   cellCounts <- table(curCells)
   keepCells <- names(cellCounts[cellCounts > 1])

   criteria$allCriteria <- data$specs$sampleInfo$sample_name[
      data$specs$sampleInfo$sample_name %in% criteria$allCriteria &
      data$specs$sampleInfo$driver%in% keepCells]

   return(criteria)

}


runModel <- function( data,
                      mode = "submit",  # "merge"
                      cluster.outDir    ) {

   if (missing(cluster.outDir)) {
      cluster.outDir <- data$specs$runModel$outDir
   }

   if (mode == "submit") {

      tx <- inferState(data = data, nSplit = 1,
                       allGenes         = TRUE,
                       sampleList       = intersect(
                        data$specs$sampleList$good$allCriteria,
                        data$specs$sampleInfo$sample_name[
                           data$specs$sampleInfo$driverType != "dissected"]),
                       sampleList.minmax = intersect(
                        data$specs$sampleList$good$allCriteria,
                        data$specs$sampleInfo$sample_name[
                           data$specs$sampleInfo$driverType == "dissected"]),
                       nlMode = "raw",
                       cluster.numTasks = data$specs$runModel$numTasks,
                       cluster.submit   = TRUE,
                       cluster.outDir   = cluster.outDir) ;

   } else if (mode == "merge") {

      tx <- inferState(data = data, nSplit = 1,
                       allGenes         = TRUE,
                       cluster.merge    = TRUE,
                       sampleList       = intersect(
                        data$specs$sampleList$good$allCriteria,
                        data$specs$sampleInfo$sample_name[
                           data$specs$sampleInfo$driverType != "dissected"]),
                       sampleList.minmax = intersect(
                        data$specs$sampleList$good$allCriteria,
                        data$specs$sampleInfo$sample_name[
                           data$specs$sampleInfo$driverType == "dissected"]),
                       nlMode = "raw",
                       cluster.outDir   = cluster.outDir) ;

   } else if (mode == "test.submit") {

      tx <- inferState(data = data, nSplit = 1,
                       allGenes         = FALSE,
                       geneList         = data$specs$inferState$genes.testPanel,
                       sampleList       = intersect(
                        data$specs$sampleList$good$allCriteria,
                        data$specs$sampleInfo$sample_name[
                           data$specs$sampleInfo$driverType != "dissected"]),
                       sampleList.minmax = intersect(
                        data$specs$sampleList$good$allCriteria,
                        data$specs$sampleInfo$sample_name[
                           data$specs$sampleInfo$driverType == "dissected"]),
                       nlMode = "raw",
                       cluster.numTasks = 2,
                       cluster.submit   = TRUE,
                       cluster.outDir   = cluster.outDir) ;

   } else if (mode == "test.merge") {

      tx <- inferState(data = data, nSplit = 1,
                       allGenes         = FALSE,
                       geneList         = data$specs$inferState$genes.testPanel,
                       cluster.merge    = TRUE,
                       sampleList       = intersect(
                        data$specs$sampleList$good$allCriteria,
                        data$specs$sampleInfo$sample_name[
                           data$specs$sampleInfo$driverType != "dissected"]),
                       sampleList.minmax = intersect(
                        data$specs$sampleList$good$allCriteria,
                        data$specs$sampleInfo$sample_name[
                           data$specs$sampleInfo$driverType == "dissected"]),
                       nlMode = "raw",
                       cluster.outDir   = cluster.outDir) ;

   }

   return(tx)

}


inferState <- function(data,

# gene selection options
                       geneList,
                       allGenes         = FALSE,

# sample selection options
                       sampleList,
                       sampleList.minmax = FALSE, #samples to consider for adding min/max dummys

# STAN options
                       nIter = 500,

# Level 1 results if done
                       geneEfit,
                       geneMode = "gene", # "gene" or "transcript"

# Return options
                       return.nlMats    = FALSE,
                       return.namesOnly = FALSE,

# gene-specific fit options
                       nlMode           = "raw", # raw, qnl, batch, qnl.batch, batch.qnl
                       level1.mode      = "stan",
                       level1.bim.stanFn,
                       level1.uni.stanFn,
                       level1.bim.cv.stanFn,
                       level1.uni.cv.stanFn,

# Parallelize: multiCPU
                       nSplit = 1, #assume 1 machine = 1 x 4-chain runs

# Parallelize: cluster
                       cluster.numTasks, #number of tasks to split on cluster
                       cluster.taskNum,  #current task number for cluster split
                       cluster.manager = "LSF",
                       cluster.merge,
                       cluster.submit,
                       cluster.outDir
   ){

# Purpose: infer expression state (on or off) from TPM matrix

   library(matrixStats)
   library(mclust)
   library(limma)
   library(preprocessCore)
   library(txtplot)

# Note: don't need paralellization for mclust, works fast.
   library(parallel); options(mc.cores=1)

# Part 0. Setup O(g,s) matrix input

   if (missing(sampleList)) {
      sampleList <- data$specs$sampleList$good$allCriteria }

   inDat <- inferState.prepInput(data,
                                 sampleList        = sampleList,
                                 geneMode          = geneMode,
                                 sampleList.minmax = sampleList.minmax,
                                 minExprThresh     = data$specs$inferState$minMaxTPM)
   
   E_gs <- inDat$O

   if (missing(geneList)) {
      if (allGenes) {
         geneList <- sort(rownames(E_gs))
      } else {
         geneList <- data$specs$inferState$genes.testPanel
      }
   }


# CLUSTER: if in cluster run mode, figure out right subset to analyze
   if (! missing(cluster.taskNum)) {

# Make sure the cluster directory option exists
      if (missing(cluster.outDir)) {
         print("ERROR: must specify cluster.outDir option")
         return(1)
      } else if (! file.exists(cluster.outDir)) {
         dir.create(cluster.outDir,recursive=TRUE)
      }

      numGenes <- length(geneList)

      if (numGenes %% cluster.numTasks == 0) {
         chunkSizes <- rep( ceiling(numGenes / cluster.numTasks),
                            times=cluster.numTasks)
      } else {
         t.n1 <- numGenes - cluster.numTasks *
                            floor(numGenes / cluster.numTasks)

         chunkSizes <- c(rep( ceiling(numGenes / cluster.numTasks), times=t.n1),
                        rep( floor(numGenes / cluster.numTasks),
                              times=cluster.numTasks - t.n1))
      }

      print(paste0("Set up cluster run for ",numGenes,
                   " genes in ",cluster.numTasks," tasks.",
                   " on task # ",cluster.taskNum))

      chunkStart <- c(1, 1 + cumsum(chunkSizes[-length(chunkSizes)]))
      chunkEnd <- cumsum(chunkSizes)
      if (chunkEnd[length(chunkEnd)] < numGenes) {
         print(paste0("ERROR: only processing ",chunkEnd," of ",
                      numGenes, "genes"))
         cat(paste0(" task range ",1:cluster.numTasks,": ",
                    chunkStart,"-",chunkEnd,"\n"))
         return(1)
      }

      geneList <- geneList[chunkStart[cluster.taskNum]:chunkEnd[cluster.taskNum]]
      print(paste0("current cluster task will process ",length(geneList)," genes"))
      print(paste0(" = genes #",chunkStart[cluster.taskNum],"-",chunkEnd[cluster.taskNum]))

   }

# remove batch effect using all genes first
   print("Normalizing expression matrix")
   nlMats <- list()
   nlMats$raw           <- log1p(E_gs)
   if (nlMode != "raw") {
      nlMats$qnl           <- log1p(normalizeQuantiles(E_gs));
      nlMats$batch         <- removeBatchEffect(nlMats$raw, covariate=inDat$D);
      nlMats$qnl.batch     <- removeBatchEffect(nlMats$qnl, covariate=inDat$D) ;

      nlMats$rqnl.batch    <- removeBatchEffect(
                              log1p(normalize.quantiles.robust(
                                 E_gs,
                                 n.remove = (ncol(nlMats$raw) - 50),
                                 remove.extreme = "both")),
                              covariate=inDat$D) ;

      rownames(nlMats$rqnl.batch) <- rownames(E_gs)
      colnames(nlMats$rqnl.batch) <- colnames(E_gs)

      nlMats$batch.qnl     <- log1p(normalizeQuantiles(exp(nlMats$batch))) ;
   }

   if (return.nlMats) { return(nlMats) ; }

# logC = corrected log-scale expression matrix; raw if nlMode == raw
   logC                 <- nlMats[[nlMode]]
   rownames(logC)       <- rownames(E_gs)
   colnames(logC)       <- colnames(E_gs)

   inDat$logC           <- logC[geneList,,drop=FALSE]
#            return(inDat)

   nGenes       <- nrow(inDat$logC); genes <- rownames(inDat$logC)
   nSamples     <- ncol(inDat$logC); samples <- colnames(inDat$logC)
   nCells       <- length(unique(inDat$cell)); cells <- inDat$cell
   nDrivers     <- length(unique(inDat$driver)); drivers <- inDat$driver



   print("Running level 1 fit")
   print(paste0("nGenes  = ",nGenes))
   print(paste0("nSamples  = ",nSamples))
   print(paste0("nDrivers == ",nDrivers))
   print(paste0("nCells = ",nCells))

   if (return.namesOnly) {
      return(list(samples       = samples,
                  cells         = unique(inDat$cell))) }

# Goal is to populate these probability matrices:
   p_on_bimodal_gc      <- matrix(nrow = nGenes,
                                  ncol = nCells,
                                  dimnames = list(genes,unique(cells)))

   p_on_bimodal_gd      <- matrix(nrow = nGenes,
                                  ncol = nDrivers,
                                  dimnames = list(genes,unique(drivers)))

   p_on_bimodal_gs      <- matrix(nrow = nGenes,
                                  ncol = nSamples,
                                  dimnames = list(genes,samples))

   diff_elpd_g          <- setNames(vector(length=nGenes),genes)
   diff_elpd_se_g       <- setNames(vector(length=nGenes),genes)

   logp_on_unimodal_g   <- setNames(vector(length=nGenes),genes)


   logp_on_gs           <- matrix(nrow = nGenes,
                                  ncol = nSamples,
                                  dimnames = list(genes,samples))

   logp_on_gd           <- matrix(nrow = nGenes,
                                  ncol = nDrivers,
                                  dimnames = list(genes,unique(drivers)))

   logp_on_gc           <- matrix(nrow = nGenes,
                                  ncol = nCells,
                                  dimnames = list(genes,unique(cells)))

   print(paste0("Part 1. Fit gene-specific 2-component mixture of logC: ",
                date()))

   if (level1.mode == "mclust") {

      geneEfit <- list()
      geneEfit <- mclapply(genes,
                           function(xx, logC){ Mclust(logC[xx,], G=2)},
                           logC)
      names(geneEfit) <- genes

      for (g in genes) {
         if (is.null(geneEfit[[g]] )) {
            p_on_bimodal_gs[g,] <- 0
         } else {
            p_on_bimodal_gs[g,] <- geneEfit[[g]]$z[,2]
         }
      }

   } else {

      if (missing(level1.bim.stanFn)) {
         level1.bim.stanFn <- data$specs$inferState$level1.bim.stanFn }

      if (missing(level1.bim.cv.stanFn)) {
         level1.bim.cv.stanFn <- data$specs$inferState$level1.bim.cv.stanFn }

      if (missing(level1.uni.stanFn)) {
         level1.uni.stanFn <- data$specs$inferState$level1.uni.stanFn }

      if (missing(level1.uni.cv.stanFn)) {
         level1.uni.cv.stanFn <- data$specs$inferState$level1.uni.cv.stanFn }

      library(rstan)
      library(loo)

# Goal: parallelize this; run first one so have compiled STAN model to use

      if (missing(geneEfit)) {

         geneEfit <- list()

         if (!missing(cluster.merge)){
            print("Merging cluster results")

            geneEfit <- inferState.mergeCluster(cluster.outDir)

         } else if (!missing(cluster.submit)){
            print("Submitting to cluster")

# Write necessary input for inferState jobs on cluster to disk
            clusterOpts <- list(
               data             = data,
               geneList         = geneList,
               geneMode         = geneMode,
               sampleList       = sampleList,
               sampleList.minmax= sampleList.minmax,
               nIter            = nIter,
               nSplit           = nSplit,
               nlMode           = nlMode,
               level1.mode      = level1.mode,
               level1.bim.stanFn = level1.bim.stanFn,
               level1.uni.stanFn = level1.uni.stanFn,
               level1.bim.cv.stanFn = level1.bim.cv.stanFn,
               level1.uni.cv.stanFn = level1.uni.cv.stanFn,
               cluster.manager  = "LSF",
               cluster.numTasks = cluster.numTasks,
               cluster.outDir   = cluster.outDir
            )

            sgeJob <- inferState.submitCluster(clusterOpts = clusterOpts)
            return(sgeJob)

# then use do.call to run inferState with the right options

            # data, nSplit, geneList, cluster.numTasks, cluster.outDir


         } else {

            print(paste0("Compiling STAN models on the first gene:", genes[1]))
            print(paste0("- compiling bimodal full model"))
            stanOptions <- list()
            stanOptions$bim <- list(
               file        = level1.bim.stanFn,
               data        = list(
                  nSamples =  nSamples,
                  nCells   =  nCells,
                  nDrivers =  nDrivers,
                  logE     =  inDat$logC[genes[1],],
                  cell     = as.integer(factor(inDat$cell, levels=unique(inDat$cell))),
                  driver   = as.integer(factor(inDat$driver, levels=unique(inDat$driver)))
               ),
               iter        = nIter,
               chains      = 4,
               cores       = 4,
               verbose     = FALSE,
               control     = list(adapt_delta = 0.98)
            )

            stanOptions$bim$fit <- do.call(stan, stanOptions$bim)

            print(paste0("- compiling unimodal full model"))
            stanOptions$uni <- stanOptions$bim
            stanOptions$uni$file <- level1.uni.stanFn
            stanOptions$uni$fit <- NULL
            stanOptions$uni$fit <- do.call(stan, stanOptions$uni)

            print(paste0("ABOUT TO PARTITION k = ",nSamples))
            kpart1 <- inferState.partitionK(n           = nSamples,
                                            k           = data$specs$inferState$kfoldvalidation.k,
                                            skipThese = c( nSamples, 
                                                           nSamples - 1),
                                            labels      = as.integer(as.factor(
                                                            inDat$driver)) )
            print("---> DONE PARTITIONING!")

            print(paste0("- compiling bimodal CV model"))
            stanOptions$bim.cv <- stanOptions$bim
            stanOptions$bim.cv$file <- level1.bim.cv.stanFn
            stanOptions$bim.cv$data <- list(
                  nSamples_t = length(which(kpart1 != 1)),
                  nSamples_h = length(which(kpart1 == 1)),
                  logE_t     =  inDat$logC[genes[1], which(kpart1 != 1)],
                  logE_h     =  inDat$logC[genes[1], which(kpart1 == 1)])
            stanOptions$bim.cv$fit <- NULL
            stanOptions$bim.cv$fit <- do.call(stan, stanOptions$bim.cv)

            print(paste0("- compiling unimodal CV model"))
            stanOptions$uni.cv <- stanOptions$bim.cv
            stanOptions$uni.cv$file <- level1.uni.cv.stanFn
            stanOptions$uni.cv$fit <- NULL
            stanOptions$uni.cv$fit <- do.call(stan, stanOptions$uni.cv)

            print("Now onto the rest!")
            startTime <- Sys.time()
            geneEfit <- mclapply(
               genes[1:length(genes)],

               function(x, stanOptions, logC) {

                  curNsamples <- ncol(logC)

                  print(paste0("NOW ON ",x))
                  stanOptions$bim$data$logE <- logC[x,];
                  stanOptions$uni$data$logE <- logC[x,];

      # 2. inferState.getStanFitPars generalize so works for bim and uni
      # 3. double check elpd math
      # 4. test on ninaE.

                  # full data fit 
                  print(paste0("FULL FIT bimodal for ",x,"!"))
                  print(paste0("-> USING STAN FILE: ", stanOptions$bim$file))
                  curBimFit <- do.call(stan, stanOptions$bim) ;
                  print("-> extracting FULL FIT bimodal parameters!")
                  bimPars <- inferState.getStanFitPars(data = stanOptions$bim$data,
                                                       fit  = curBimFit,
                                                       mode  = "bimodal")
                  print("-> DONE")

                  print(paste0("FULL FIT unimodal for ",x,"!"))
                  print(paste0("-> USING STAN FILE: ", stanOptions$uni$file))
                  curUniFit <- do.call(stan, stanOptions$uni) ;
                  uniPars <- inferState.getStanFitPars(data = stanOptions$uni$data,
                                                       fit  = curUniFit,
                                                       mode ="unimodal")

# 1. paritionK code
                  numK <- data$specs$inferState$kfoldvalidation.k
                  kparts <- inferState.partitionK(n     = curNsamples,
                                                  k     = numK,
                                                  skipThese = c(
                                                      curNsamples, 
                                                      curNsamples - 1),
                                                  labels = as.integer(as.factor(
                                                            inDat$driver)) )
                  bim.elpd_i <- vector(length=curNsamples,"numeric")
                  uni.elpd_i <- vector(length=curNsamples,"numeric")

#DEBUG                  d.hLists <- list(); d.tLists <- list() ;
                  for (k in 1:numK) {
                     print(paste("on K partition : ",k," for gene ",x))

                     hList <- which(kparts == k)
                     tList <- which(kparts != k)

#DEBUG                     d.hLists[[k]] <- hList
#DEBUG                     d.tLists[[k]] <- tList

#                     next; #DEBUG 180202_0954

                     stanOptions$bim.cv$data <- list(
                        logE_h <- logC[x,hList],
                        logE_t <- logC[x,tList],
                        nSamples_h <- length(hList),
                        nSamples_t <- length(tList))

                     stanOptions$uni.cv$data <- stanOptions$bim.cv$data

                     print("calling bim.cv")
                     bimFit.cv <- do.call(stan, stanOptions$bim.cv) ;

                     print("calling uni.cv")
                     uniFit.cv <- do.call(stan, stanOptions$uni.cv) ;

                     print("extracting log lik")
                     bim.kll <- extract_log_lik(bimFit.cv,
                                                parameter_name="log_lik_h")
                     uni.kll <- extract_log_lik(uniFit.cv,
                                                parameter_name="log_lik_h")

                     for (h in 1:length(hList)) {
                        bim.elpd_i[hList[h]] <- logSumExp(bim.kll[h,]) -
                                                log(nrow(bim.kll)) ;

                        uni.elpd_i[hList[h]] <- logSumExp(uni.kll[h,]) -
                                                log(nrow(uni.kll)) ;
                     }
                  }

#DEBUG                   return(list(d.hLists=d.hLists, d.tLists=d.tLists, inDat))

                  bim.elpd <- sum(bim.elpd_i)
                  uni.elpd <- sum(uni.elpd_i)
                  diff.elpd <- bim.elpd - uni.elpd
                  print(paste0("bim.elpd=",bim.elpd))
                  print(paste0("uni.elpd=",uni.elpd))

                  bim.elpd.se  <- sqrt(curNsamples * var(bim.elpd_i))
                  uni.elpd.se  <- sqrt(curNsamples * var(uni.elpd_i))
                  diff.elpd.se <- sqrt(curNsamples * var(bim.elpd_i - uni.elpd_i))

                  return( list( bimPars      = bimPars,
                                uniPars      = uniPars,
                                bim.elpd_i   = bim.elpd_i,
                                uni.elpd_i   = uni.elpd_i,
                                bim.elpd     = bim.elpd,
                                uni.elpd     = uni.elpd,
                                bim.elpd.se  = bim.elpd.se,
                                uni.elpd.se  = uni.elpd.se,
                                diff.elpd    = diff.elpd,
                                diff.elpd.se = diff.elpd.se ) )
               },
               stanOptions      = stanOptions,
               logC             = inDat$logC,
               mc.cores         = nSplit,
               mc.allow.recursive = TRUE,
               mc.preschedule   = F
            )
            endTime <- Sys.time()
            lenTime <- endTime - startTime
            print(paste0("Fit ",length(genes)," genes in ",lenTime))
            print(paste0("   started ",startTime," and ended ",endTime))
            names(geneEfit) <- genes

#            print("DEBUG: CHECK geneEfit")
#            return(geneEfit)

            curRes <- list(geneEfit     = geneEfit,
                           geneList     = geneList,
                           inDat        = inDat)

            if (!missing(cluster.numTasks)) {
               outFn <- paste0(cluster.outDir,"/inferState_task_",cluster.taskNum,
                               "_of_",cluster.numTasks,".rds")
               saveRDS(curRes$geneEfit, file=outFn)
               return(1)
            } else {
               return(curRes)
            }
         }
      }

      for (g in 1:nGenes) {
         if (! "diff.elpd" %in% names(geneEfit[[genes[g]]])) {
            print(paste0("SKIPPING: ",geneList[genes[g]]," since missing diff.elpd esitmate"))
            next;
         }
         diff_elpd_g[g] <- geneEfit[[genes[g]]]$diff.elpd
         diff_elpd_se_g[g] <- geneEfit[[genes[g]]]$diff.elpd.se
         p_on_bimodal_gc[g,] <- exp(geneEfit[[genes[g]]]$bimPars$p_on_bimodal_gc)
         p_on_bimodal_gd[g,] <- exp(geneEfit[[genes[g]]]$bimPars$p_on_bimodal_gd)
         p_on_bimodal_gs[g,] <- exp(geneEfit[[genes[g]]]$bimPars$p_on_bimodal_gs)
      }

   }



   print(paste0("Part 2. Model global distributions of p(E|off) and p(E|on) ",
                "from bimodal genes (diff_elpd > 2 * diff_elpd_se) ",date()))
   {
      t.bimGenes <- names(diff_elpd_g)[which(diff_elpd_g > 2 * diff_elpd_se_g)]
      print(paste0("-> n bimodal genes=",length(t.bimGenes)))
      bim.logC <- logC[t.bimGenes,]
      bim.p_on_gc <- p_on_bimodal_gc[t.bimGenes,]

      coarse.p_on_bimodal_gs <- matrix(nrow = nGenes,
                                       ncol = nSamples,
                                       dimnames = list(genes,samples))
      for (s in 1:nSamples) {
         coarse.p_on_bimodal_gs[,s] <- p_on_bimodal_gc[,inDat$cell[s]] }

      bim.p_on_gs <- coarse.p_on_bimodal_gs[t.bimGenes,]

      real.onE <- as.vector(unlist(
                     lapply(1:nrow(bim.logC), function(x,t1,t2){
                         t1[x,which(t2[x,] > 0.9)]},
                         t1 = bim.logC, t2 = bim.p_on_gs)))

      txtdensity(real.onE, width=80)

      real.offE <- as.vector(unlist(
                     lapply(1:nrow(bim.logC), function(x,t1,t2){
                         t1[x,which(t2[x,] < 0.1)]},
                         t1 = bim.logC, t2 = bim.p_on_gs)))
      txtdensity(real.offE, width=80)

   }

   print("FITTING ON DISTRIBUTION")
   uniEonfit            <- Mclust(real.onE,G=1)
   print("- >DONE")

   uniEonfit$para$sd    <- sqrt(uniEonfit$para$variance$sigmasq)


   print("FITTING off DISTRIbutION")
   uniEofffit           <- Mclust(real.offE,G=1)
   print("- >DONE")
   uniEofffit$para$sd   <- sqrt(uniEofffit$para$variance$sigmasq)
   pi_on.uni            <- length(real.onE) /
                           (length(real.onE) + length(real.offE))
   logpi_on.uni         <- log(pi_on.uni)
   logpi_off.uni        <- log1p(-1 * pi_on.uni)


   print(paste0("uniEon: mean=",uniEonfit$para$mean," sd=",uniEonfit$para$sd))
   print(paste0("uniEoff: mean=",uniEofffit$para$mean," sd=",uniEofffit$para$sd))

   print(paste0("Part 2b.  Estimate p(on|unimodal) for all genes.",date()))
   i <- 1
   for (g in geneList) {
      logp_numer <- 0 ; logp_denom <- 0 ;
#      print(paste0("now on ",g))
      for (s in samples) {
         t.on <- logpi_on.uni + dnorm(logC[g,s],
                                      mean = uniEonfit$para$mean,
                                      sd   = uniEonfit$para$sd,
                                      log  = TRUE)
         logp_numer <- logp_numer + t.on

         t.off <- logpi_off.uni + dnorm(logC[g,s],
                                        mean    = uniEofffit$para$mean,
                                        sd      = uniEofffit$para$sd,
                                        log     = TRUE)
         logp_denom <- logp_denom + t.off

      }
      logp_on_unimodal_g[g] <- logp_numer - logSumExp(c(logp_denom, logp_numer))
#      print(paste0("logp_on_unimodal_g[",g,"] = ", logp_on_unimodal_g[g]))
      if (i %% 500 == 0) { print(paste0("done with ",i," genes")) }
      i <- i + 1
   }

   print(paste0("Part 3. Estimate P(on_gs) by choosing unimodal or bimodal ",
                 date()))
   i <- 1
   for (g in geneList) {

      if (diff_elpd_g[g] > 0 &
          diff_elpd_g[g] > 2 * diff_elpd_se_g[g]) { #NEW BIT 171023_1205

         geneEfit[[g]]$exprType <- "bimodal"

         logp_on_gs[g,] <- log(p_on_bimodal_gs[g,]) ;
         logp_on_gd[g,] <- log(p_on_bimodal_gd[g,]) ;
         logp_on_gc[g,] <- log(p_on_bimodal_gc[g,]) ;

      } else {

         geneEfit[[g]]$exprType <- "unimodal"

         logp_on_gs[g,] <- logp_on_unimodal_g[g] ;
         logp_on_gd[g,] <- logp_on_unimodal_g[g] ;
         logp_on_gc[g,] <- logp_on_unimodal_g[g] ;

      }
      if (i %% 500 == 0) { print(paste0("done with ",i," genes")) }
      i <- i + 1
   }
   print(paste0("Finished: ", date()))

   return(list(
      uniFit = list(uniEonfit   = uniEonfit,
                    uniEofffit  = uniEofffit,
                    pi_on.uni   = pi_on.uni),
      inDat             = inDat,
      logC              = logC,
      geneEfit          = geneEfit,
      p_on_bimodal_gc   = p_on_bimodal_gc,
      p_on_bimodal_gd   = p_on_bimodal_gd,
      p_on_bimodal_gs   = p_on_bimodal_gs,
      logp_on_unimodal_g = logp_on_unimodal_g,
      logp_on_gc        = logp_on_gc,
      logp_on_gd        = logp_on_gd,
      logp_on_gs        = logp_on_gs
   ))

}


inferState.prepInput <- function(data,
                                 geneMode = "gene",
                                 sampleList,
                                 sampleList.minmax = FALSE,
                                 minExprThresh) {

# Purpose: prepares input for inferState() run

   if (geneMode == "transcript") {
      dataMat <- "txExpr"
      nameCol <- "transcript_id"
   } else if (geneMode == "gene") {
      dataMat <- "geneExpr"
      nameCol <- "gene_name"
   }

# Expr matrix
   obs_gs <- data$expr[[dataMat]][,paste0("tpm.",sampleList)]
   obs_gs <- as.matrix(obs_gs)

   rownames(obs_gs) <- data$expr[[dataMat]][,nameCol]


# Filter out genes expressed at low levels in all samples
   print(paste0("Ignoring genes with max(TPM) < ",minExprThresh))
   print(paste0("-> original: ", nrow(obs_gs)," genes"))
   obs_gs <- obs_gs[apply(obs_gs,1,max) >= minExprThresh,]
   print(paste0("->  reduced: ", nrow(obs_gs)," genes"))

# Sample names
   sampleNames <- colnames(obs_gs)
   nSamples <- ncol(obs_gs)
   nGenes <- nrow(obs_gs)

# Add protocol info vector
   cellList <- c()
   cellListNames <- c()
   nCells <- 0

   driverList <- c()
   driverListNames <- c()
   nDrivers <- 0
   {
      driverList <- data$specs$sampleInfo$driver[
          match(gsub("tpm.","",colnames(obs_gs)),
               data$specs$sampleInfo$sample_name)]

      cellList <- data$specs$sampleInfo$celltype[
          match(gsub("tpm.","",colnames(obs_gs)),
               data$specs$sampleInfo$sample_name)]

      fullCellList <- cellList
#      print(paste("cells: ",unique(cellList)))
      cellListNames <- unique(cellList)
      cellList <- as.integer(as.factor(cellList))
      nCells<-length(unique(cellList))

      fullDriverList <- driverList
#      print(paste("drivers: ",unique(driverList)))
      driverListNames <- unique(driverList)
      driverList <- as.integer(as.factor(driverList))
      nDrivers <- length(unique(driverList))

   }

   if (!is.logical(sampleList.minmax)) {
      print("GOT IN LOGICAL!")
      minMaxMat <- data.frame(
         tpm.min = apply(data$expr[[dataMat]][
                           match(rownames(obs_gs),data$expr[[dataMat]][,nameCol]),
                           paste0("tpm.", union(sampleList, sampleList.minmax)) ], 1, min),

         tpm.max = apply(data$expr[[dataMat]][
                           match(rownames(obs_gs),data$expr[[dataMat]][,nameCol]),
                           paste0("tpm.", union(sampleList, sampleList.minmax))], 1, max)
      )
      obs_gs <- as.matrix(cbind(obs_gs, minMaxMat))
      nSamples <- nSamples + 2
      nCells <- nCells + 2
      nDrivers <- nDrivers + 2
      fullCellList <- c(fullCellList, "min", "max")
      fullDriverList <- c(fullDriverList, "min", "max")
   }

   print(paste0("------ nDrivers= ",nDrivers))
   print(paste0("------ nCells= ",nCells))
   print(paste0("------ nSamples= ",nSamples))
   print(fullCellList)
                                             

   return(list(nSamples      = nSamples,
               nCells        = nCells,
               nDrivers      = nDrivers,
               nGenes        = nGenes,
               O             = obs_gs,
               cell          = fullCellList,
               driver        = fullDriverList
   ))

}


inferState.submitCluster <- function(clusterOpts) {
# Purpose: Sets up R script, input, and SGE script to run inferState() job on cluster

# write necessary input data to rds file for loading

   print("Saving input files for cluster nodes")
   if (!file.exists(clusterOpts$cluster.outDir)) {
         dir.create(clusterOpts$cluster.outDir, recursive = TRUE) }
   fn.clusterOpts.rds <- paste0(clusterOpts$cluster.outDir,"/inferState_cluster_input.rds")
   if (!file.exists(fn.clusterOpts.rds)) {
      saveRDS(clusterOpts, file=fn.clusterOpts.rds)
   } else {
      print("Already done")
   }
   print("Finished, now submitting job:")


# make R script to load necessary data and run subjobs.
   fn.Rscript <- paste0(clusterOpts$cluster.outDir,"/inferState_cluster_job.R")
   curFh <- file(fn.Rscript)
   writeLines(c(
'source("~/ol/src/R/analyzeOpticLobeExpr.R");',
'args = commandArgs(trailingOnly=TRUE);',
'curTaskNum <- as.numeric(args[1])',
paste0('clusterOpts <- readRDS("',fn.clusterOpts.rds,'");'),
'do.call(inferState, c(clusterOpts, cluster.taskNum=curTaskNum));'
), curFh)
   close(curFh)

# Make SGE output directory
   sgeOutDir <- "SGEout"
   sgeOutDir.full <- paste0(clusterOpts$cluster.outDir,"/",sgeOutDir)
   if (!file.exists(sgeOutDir.full)) dir.create(sgeOutDir.full, recursive = TRUE)

   if (clusterOpts$cluster.manager == "SGE") {

# Write SGE script

      fn.SGEscript <- paste0(clusterOpts$cluster.outDir,"/inferState_cluster_job.SGE.sh")
      curFh <- file(fn.SGEscript)
      writeLines(c(
'#!/bin/csh',
'#$ -S /bin/csh',
'#$ -cwd',
paste0('#$ -o ',sgeOutDir),
paste0('#$ -e ',sgeOutDir),
'#$ -pe batch 4',
paste0('#$ -t 1-',clusterOpts$cluster.numTasks),
'source /sge/current/default/common/settings.csh',
paste0('hostname'),
paste0('echo /groups/eddy/home/davisf/software/R/R-3.3.1/bin/Rscript --vanilla ',
       fn.Rscript,' $SGE_TASK_ID'),
paste0('date'),
paste0('/groups/eddy/home/davisf/software/R/R-3.3.1/bin/Rscript --vanilla ',
       fn.Rscript,' $SGE_TASK_ID'),
paste0('date')
), curFh)
      close(curFh)

      tCom <- paste0("qsub -wd ",clusterOpts$cluster.outDir," ",fn.SGEscript)
      sgeRes <- system(tCom, intern=TRUE)

   } else if (clusterOpts$cluster.manager == "LSF") {

      fn.LSFscript <- paste0(clusterOpts$cluster.outDir,"/inferState_cluster_job.LSF.sh")
      print(paste0("Writing LSF script here: ", fn.LSFscript))
      curFh <- file(fn.LSFscript)
   writeLines(c(
'#!/bin/csh',
#paste0('#BSUB -o ',sgeOutDir.full,'/%J.%I.out'),
#paste0('#BSUB -e ',sgeOutDir.full,'/%J.%I.err'),
paste0('#BSUB -o ',sgeOutDir,'/%J.%I.out'),
paste0('#BSUB -e ',sgeOutDir,'/%J.%I.err'),
'#BSUB -n 4',
'#BSUB -W 10:00',
'source /misc/lsf/conf/cshrc.lsf',
paste0('hostname'),
paste0('echo /groups/eddy/home/davisf/software/R/R-3.3.1/bin/Rscript --vanilla ',
       fn.Rscript,' $LSB_JOBINDEX'),
paste0('date'),
paste0('/groups/eddy/home/davisf/software/R/R-3.3.1/bin/Rscript --vanilla ',
       fn.Rscript,' $LSB_JOBINDEX'),
paste0('date')
), curFh)
      close(curFh)

      tCom <- paste0("bsub -cwd ",clusterOpts$cluster.outDir,
                     " -J \"inferState_cluster_job[1-",
                     clusterOpts$cluster.numTasks,"]\" < ",fn.LSFscript)
      sgeRes <- system(tCom, intern=TRUE)
   }

   return(sgeRes)

}


inferState.mergeCluster <- function(outputDir) {
# Purpose: merges result of cluster-parallelized inferState() run

# assume chunks only wrote geneEfit, add on inDat 

#inferState_task_#_of_#.rds
   print(paste0("Merging subjobs in ",outputDir))
   chunkFns <- list.files(path = outputDir, pattern = "inferState_task.*.rds")
   print(chunkFns)
   chunkOrder <- gsub("_of.*", "",chunkFns)
   chunkOrder <- gsub("*_task_", "",chunkFns)
   chunkOrder <- as.numeric(chunkOrder)
   chunkOrder <- order(chunkOrder)
   for (i in 1:length(chunkFns)) {
      curFn <- paste0(outputDir,"/",chunkFns[chunkOrder[i]])
      print(paste0("- reading split: ",curFn))
      if (i == 1) {
         mergeRes <- readRDS(curFn)
      } else {
         curChunk <- readRDS(curFn)
         mergeRes <- c(mergeRes, curChunk)
      }
   }
   saveRDS(mergeRes, file=paste0(outputDir,"/merged_results.rds"))
   return(mergeRes)

}


inferState.getStanFitPars <- function(data, fit, mode) {
# Purpose: extracts STAN fit parameters

   curPars <- list()

   if (mode == "cv") {

      curPars$log_lik_t <- extract_log_lik(fit, parameter_name="log_lik_t")
      curPars$log_lik_h <- extract_log_lik(fit, parameter_name="log_lik_h")

   } else if (mode == "bimodal") {

      stansum <- summary(fit)$summary ;
#      curPars$fullSum <- stansum

      curPars$p_on_bimodal_gc           <- stansum[
         paste0("pon_gc[",1:data$nCells,"]"),"50%"]
      curPars$p_on_bimodal_gd           <- stansum[
         paste0("pon_gd[",1:data$nDrivers,"]"),"50%"]
      curPars$p_on_bimodal_gs           <- stansum[
         paste0("pon_gs[",1:data$nSamples,"]"),"50%"]
      curPars$mu_on                     <- stansum["mu_on","50%"]
      curPars$mu_off                    <- stansum["mu_off","50%"]
      curPars$sd_on                     <- stansum["sd_on","50%"]
      curPars$sd_off                    <- stansum["sd_off","50%"]
      curPars$pi_on                     <- stansum["pi_on","50%"]

      curPars$rhat.p_on_bimodal_gc      <- stansum[
         paste0("pon_gc[",1:data$nCells,"]"),"Rhat"]
      curPars$rhat.p_on_bimodal_gd      <- stansum[
         paste0("pon_gd[",1:data$nDrivers,"]"),"Rhat"]
      curPars$rhat.p_on_bimodal_gs      <- stansum[
         paste0("pon_gs[",1:data$nSamples,"]"),"Rhat"]
      curPars$rhat.mu_on                <- stansum["mu_on","Rhat"]
      curPars$rhat.mu_off               <- stansum["mu_off","Rhat"]
      curPars$rhat.sd_on                <- stansum["sd_on","Rhat"]
      curPars$rhat.sd_off               <- stansum["sd_off","Rhat"]
      curPars$rhat.pi_on                <- stansum["pi_on","Rhat"]

      curPars$neff.p_on_bimodal_gc      <- stansum[
         paste0("pon_gc[",1:data$nCells,"]"),"n_eff"]
      curPars$neff.p_on_bimodal_gd      <- stansum[
         paste0("pon_gd[",1:data$nDrivers,"]"),"n_eff"]
      curPars$neff.p_on_bimodal_gs      <- stansum[
         paste0("pon_gs[",1:data$nSamples,"]"),"n_eff"]
      curPars$neff.mu_on                <- stansum["mu_on","n_eff"]
      curPars$neff.mu_off               <- stansum["mu_off","n_eff"]
      curPars$neff.sd_on                <- stansum["sd_on","n_eff"]
      curPars$neff.sd_off               <- stansum["sd_off","n_eff"]
      curPars$neff.pi_on                <- stansum["pi_on","n_eff"]

   } else if (mode == "unimodal") {

      stansum <- summary(fit)$summary ;
      print("DEBUGGING unimodal par extraction!!")
      print(rownames(stansum))

      curPars$mu1       <- stansum["mu1","50%"]
      curPars$sd1       <- stansum["sd1","50%"]
      curPars$rhat.mu1  <- stansum["mu1","Rhat"]
      curPars$rhat.sd1  <- stansum["sd1","Rhat"]

      curPars$neff.mu1  <- stansum["mu1","n_eff"]
      curPars$neff.sd1  <- stansum["sd1","n_eff"]

   }

   return(curPars)
}


inferState.partitionK <- function(n, k, labels, skipThese) {
# Purpose: Partitions a set of points for k-fold cross-validation
#
# returns: array of k-partition assignment

   fullList <- 1:n
   if (!missing(skipThese)) {
      fullList <- setdiff(fullList, skipThese)
   } else {
      skipThese <- c()
   }

   kparts <- vector(mode="numeric",length=n)

   shuffledI <- sample(fullList, length(fullList))
   medSize <- floor(length(fullList) / k)
   for (i in 1:k) {
      t.start <- ((i - 1) * medSize) + 1
      for (j in t.start:(t.start + medSize - 1)) {
         kparts[shuffledI[j]] <- i
      }
   }

   if (length(fullList) %% k > 0) {
      for (j in (k * medSize + 1):length(fullList)) {
         kparts[shuffledI[j]] <- k
      }
   }


# If points of a particuilar label are all found in same parition, shuffle with
# another class point to reblanace

   keepLabels <- labels[setdiff(fullList, skipThese)]

   allGood <- FALSE
   while (!allGood) {
      allGood <- TRUE
      for (curLabel in unique(keepLabels)) {
         parts <- unique(kparts[labels == curLabel])
         if (length(parts) == 1) {

# swap 2 points: 1 from this label, 2 from another label, in another class.
            i <- sample(which(labels == curLabel))[1]
            j <- sample(which(labels != curLabel & kparts != parts[1]))[1]
            part.t <- kparts[i]
            kparts[i] <- kparts[j]
            kparts[j] <- part.t

            allGood <- FALSE
         } 
      }
   }

   for (i in skipThese) {
      kparts[i] <- k + 1 }

   return(kparts)

}


inferState.scoreTPM <- function(data,
                                geneList,
                                samples,
                                tpm,     # can either provide TPM list of sep. samples
                                tpm.reps = FALSE #FALSE: TPMs from different driver/cell
                                                 #TRUE: TPMs from same driver/cell
                                ) {

# Purpose: given inferState model fit, estimate p(on) sample, driver, cell for
#          other samples -- especially bulk samples not used in the model fit

   library(matrixStats)

   nGenes <- length(geneList)

# setup input as in STAN inferState model
   if (!missing(tpm)) {

      nSamples <- length(tpm)
      samples <- 1:nSamples
      logE <- log1p(tpm)

      if (tpm.reps == FALSE) {
         drivers <- 1:nSamples
         cells <- 1:nSamples
      } else {
         drivers <- rep(1, nSamples)
         cells <- rep(1, nSamples)
      }

   } else {

      nSamples <- length(samples)
      drivers <- data$specs$sampleInfo$driver[
                        match(samples, data$specs$sampleInfo$sample_name)]
      cells <- data$specs$sampleInfo$celltype[
                        match(samples, data$specs$sampleInfo$sample_name)]

      logE <- as.matrix(log1p(dat$expr$geneExpr[
                           match(geneList, dat$expr$geneExpr$gene_name),
                           paste0("tpm.",samples),
                           drop=FALSE]))
      colnames(logE) <- samples
      rownames(logE) <- geneList

   }

   drivers.uniq <- sort(unique(drivers))
   cells.uniq <- sort(unique(cells))

   nDrivers <- length(drivers.uniq)
   nCells <- length(cells.uniq)


   # goal: populate these:
   pon_gs <- matrix(nrow = nGenes, ncol=nSamples,
                    dimnames=list(geneList,samples))
   pon_gd <- matrix(nrow = nGenes, ncol=nDrivers,
                    dimnames=list(geneList,drivers.uniq))
   pon_gc <- matrix(nrow = nGenes, ncol=nCells,
                    dimnames=list(geneList,cells.uniq))

# Replicating level1_ordered.stan model
   for (gi in 1:nGenes) {
      gene <- geneList[gi]

      if (!gene %in% names(data$pFit$geneEfit)) {
         print(paste0("No inferState model found for gene ", gene))
         next;
      }

      if (data$pFit$geneEfit[[gene]]$exprType == "bimodal") {

         mu_on <- data$pFit$geneEfit[[gene]]$bimPars$mu_on
         mu_off <- data$pFit$geneEfit[[gene]]$bimPars$mu_off
         sd_on <- data$pFit$geneEfit[[gene]]$bimPars$sd_on
         sd_off <- data$pFit$geneEfit[[gene]]$bimPars$sd_off
         pi_on <- data$pFit$geneEfit[[gene]]$bimPars$pi


         logpi_on  <- log(pi_on)
         logpi_off <- log1p(-1 * pi_on)

         pon_gs_onmass  <- vector(length=nSamples,mode="numeric");
         names(pon_gs_onmass) <- samples
         pon_gs_offmass <- vector(length=nSamples,mode="numeric");
         names(pon_gs_offmass) <- samples

         pon_gd_onmass  <- vector(length=nDrivers,mode="numeric");
         names(pon_gd_onmass) <- drivers.uniq
         pon_gd_offmass <- vector(length=nDrivers,mode="numeric");
         names(pon_gd_offmass) <- drivers.uniq

         pon_gc_onmass  <- vector(length=nCells, mode="numeric");
         names(pon_gc_onmass) <- cells.uniq
         pon_gc_offmass <- vector(length=nCells, mode="numeric");
         names(pon_gc_offmass) <- cells.uniq

         for (s in 1:nSamples) {
            curSample <- samples[s]
            ps1 <- dnorm(logE[gene,curSample], mean=mu_on, sd=sd_on, log=TRUE)
            ps2 <- dnorm(logE[gene,curSample], mean=mu_off, sd=sd_off, log=TRUE)
            sumps <- logSumExp(c(logpi_on + ps1, logpi_off + ps2))

            pon_gs[gene,curSample] <- logpi_on + ps1 - sumps ;

            pon_gd_onmass[drivers[s]] <- pon_gd_onmass[drivers[s]] + ps1 ;
            pon_gd_offmass[drivers[s]] <- pon_gd_offmass[drivers[s]] + ps2 ;


            pon_gc_onmass[cells[s]] <- pon_gc_onmass[cells[s]] + ps1 ;
            pon_gc_offmass[cells[s]] <- pon_gc_offmass[cells[s]] + ps2 ;

         }

         for (ci in 1:nCells) {
            curCell <- cells.uniq[ci]

            pon_gc[gene,curCell] <- logpi_on  + pon_gc_onmass[curCell] -
                        logSumExp(c(logpi_on  + pon_gc_onmass[curCell],
                                    logpi_off + pon_gc_offmass[curCell]))
         }

         for (di in 1:nDrivers) {
            curDriver <- drivers.uniq[di]

            pon_gd[gene,curDriver] <- logpi_on  + pon_gd_onmass[curDriver] -
                          logSumExp(c(logpi_on  + pon_gd_onmass[curDriver],
                                      logpi_off + pon_gd_offmass[curDriver]))
         }

      } else { # Unimodal

         pon_gs[gene,] <- rep(data$pFit$logp_on_gs[gene,1], length=nSamples)
         pon_gd[gene,] <- rep(data$pFit$logp_on_gd[gene,1], length=nDrivers)
         pon_gc[gene,] <- rep(data$pFit$logp_on_gc[gene,1], length=nCells)

      }
   }

   return(list( pon_gs = pon_gs,
                pon_gd = pon_gd,
                pon_gc = pon_gc
               ))

}


plotExprVsYield <- function(data,
                            geneList,
                            protocolFilter,
                            biotype             = FALSE,
                            plot                = TRUE,
                            multiPlot           = FALSE,
                            score               = "tpm.",
                            dataMat             = "geneExpr",
                            runFit              = FALSE,
                            addLabels           = FALSE,
                            logX                = TRUE,
                            logY                = TRUE,
                            cdna.conc_to_ug            = .03,
                            colorProtocol       = TRUE,
                            printYieldCorrGenes = FALSE,
                            figName.expr_yield  = "expr_yield_plots",
                            figName.yield_yield = "yield_yield_plots",
                            returnData          = FALSE                 ){
# Purpose: plot expression vs yield relationship for specified genes

# alt: score="tpm.", dataMat="biotypeExpr"
# alt: score="intron_exon_fpkm_ratio", dataMat="exonIntronMat"

# alt: 

   library(calibrate)

   colorList <- c("darkOrange3","steelBlue3","purple","black")

   x <- data$specs$sampleInfo[!(data$specs$sampleInfo$celltype %in%
                                data$specs$controlDrivers),
                              c("sample_name", "yield.nuclei.thousands",
                                "yield.cdna.ng_uL", "intact.protocol")]

   x$yield.nuclei.thousands <- suppressWarnings(as.numeric(x$yield.nuclei.thousands))
   x$yield.cdna.ng_uL   <- suppressWarnings(as.numeric(x$yield.cdna.ng_uL))

   x <- x[x$sample_name %in% data$specs$sampleList$all,]

   if (!missing(protocolFilter)) x <- x[x$intact.protocol %in% protocolFilter, ]

   samplesWithYield <- sort(x$sample_name)
   yield <- x$yield.nuclei.thousands[ match(samplesWithYield, x$sample_name)]
   yield.cdna.ng_uL <- x$yield.cdna.ng_uL[ match(samplesWithYield, x$sample_name)]
   protocols <- x$intact.protocol[match(samplesWithYield, x$sample_name)]
   print("protocols with yield numbers:")
   print(paste0(unique(protocols), collapse=", "))

   protocolNames <- protocols
   protocolColors <- as.numeric(as.factor(protocols))

   if (biotype & missing(dataMat)) dataMat<-"biotypeExpr"

   tpmCols <- paste0("tpm.",samplesWithYield)
   tpm <- data$expr[[dataMat]][,tpmCols]
   if (!biotype & grepl("gene", dataMat)) {
      rownames(tpm) <- data$expr[[dataMat]]$gene_name}

   if (!biotype & grepl("tx", dataMat)) {
      rownames(tpm) <- data$expr[[dataMat]]$transcript_id}

   if (returnData) {
      return(list(
         yield = yield,
         tpm   = tpm,
         yield.cdna.ng_uL = yield.cdna.ng_uL
      ))
   }


   if (biotype) {
      rownames(tpm) <- rownames(data$expr[[dataMat]])
   } else {
      rownames(tpm) <- data$expr[[dataMat]]$gene_name
   }


   if (printYieldCorrGenes) {
      yieldScore <- apply(tpm, 1, function(a){cor.test(a, yield)$estimate})
      yieldInverseCorrGenes <- head(rownames(
         tpm[order(yieldScore, decreasing=FALSE), ]), n=100)
      yieldCorrGenes <- head(rownames(
         tpm[order(yieldScore, decreasing=TRUE), ]), n=100)

      print("Inversely yield-correlated genes:")
      print(paste(yieldInverseCorrGenes, collapse=",  "))

      print("Yield-correlated genes:")
      print(paste(yieldCorrGenes, collapse=",  "))
   }


   geneSets <- list(
      manual = c("Gad1",  "VAChT",  "VGlut", "ninaE", "VGAT",
                 "Grd", "Rdl", "bsh", "ap", "ey", "toy", "dac",
                 "Lim3", "brp", "svp", "nSyb")
   )

   if (! missing(geneList)) geneSets <- list(manual = c(geneList))

   if (multiPlot) {
      pdf(paste0(data$specs$outDir,"/",figName.expr_yield,"_expr_vs_yield.pdf"),
          height=11, width=8.5)
      par(mfrow=c(5, 3), mar=c(5, 5, 1, 6), xpd=TRUE)
   } 

   plottedYieldYield <- FALSE
   for (geneSet in names(geneSets)) {

      for (i in (1:length(geneSets[[geneSet]]))) {

         geneName <- geneSets[[geneSet]][i]

         if (! multiPlot) {
            pdf(paste0(data$specs$outDir, "/",figName.expr_yield,"_expr_vs_yield_",
                       geneSet, "_", geneName, ".pdf"),
                height = 2,
                width  = 2.5)
            par(mar=c(3,4.0,0.5,0.5),ps=10,xpd=TRUE)
         }

         curx <- as.numeric(tpm[geneName, ])
         xy <- as.data.frame(cbind(expr=curx, yield=yield, yield.cdna.ng_uL = yield.cdna.ng_uL))
         rownames(xy) <- samplesWithYield
         xy$protocolColors <- protocolColors
         xy$protocolNames <- protocolNames
         xy <- na.omit(xy)

         if (runFit) {
            print(paste0("sending gene ", geneName, " for fitContam()"))
            eyFit <- fitExprYield(expr          = xy$expr,
                                  yield         = xy$yield,
                                  geneName      = geneName,
                                  decontamMode  = "binary",
                                  includeBgLevel=FALSE,
                                  details       = TRUE)
         }

         pointCols <- "black"
         xy$protocolColors <- as.numeric(as.factor(xy$protocolNames))
         if (colorProtocol)        pointCols <- c("darkOrange3","steelBlue3")[xy$protocolColors]

         if (!plottedYieldYield) {
            plotX <- xy$yield
            plotY <- xy$yield.cdna.ng_uL * cdna.conc_to_ug
            yyPointCol <- pointCols
            yyPch <- rep(20,length(plotX))

# add mock points at position 1E-1
            mockPoints <- data$specs$sampleInfo[
                              data$specs$sampleInfo$yield.nuclei.thousands %in% c(0) &
                              data$specs$sampleInfo$intact.protocol %in%
                                 c("INTACT", "TAPIN") &
                              data$specs$sampleInfo$driver %in%
                              c("no_driver_aGFP_beads", "L1L2_aFLAG_beads",
                                "Tm1_d1_mock2", "mock2_201611"),]

            mockX <- rep(1E-1, nrow(mockPoints))
            mockY <- mockPoints$yield.cdna.ng_uL * cdna.conc_to_ug
            mockProt <- mockPoints$intact.protocol
            mockCols <- c("darkOrange3", "steelBlue3")[as.numeric(as.factor(mockProt))]
            mockPch <- rep(17, nrow(mockPoints))

            plotX <- c(plotX, mockX)
            plotY <- c(plotY, mockY)
            yyPointCol <- c(yyPointCol, mockCols)
            yyPch<- c(yyPch, mockPch)

            print("CAN WE GET AWAY WITHOUT 1+ for the YIELD-YIELD PLOT?************")
            if (0) {
            if (logX) plotX <- log10(1 + plotX)
            if (logY) plotY <- log10(1 + plotY)
            }

            if (logX) plotX <- log10(plotX)
            if (logY) plotY <- log10(plotY)
   
               pdf(paste0(data$specs$outDir, "/",figName.yield_yield,"_yield_vs_yield.pdf"),
                   height = 2,
                   width  = 2.5)
               par(mar=c(3,4.0,0.5,0.5),ps=10,xpd=TRUE)
   
            plot(plotX,
                 plotY,
                 xlab = "", ylab="",
                 main="",
                 pch          = yyPch,
                 cex          = 0.75,
                 las          = 1,
                 col          = yyPointCol,
                 yaxt         = "n",
                 xaxt         = "n")
   
            mtext(side = 1, text = "nuclear yield (thousands)", line = 1.5)
            mtext(side = 2, text = "cDNA yield (ug)", line = 2)
   
            labels <- c(1, 5, 10, 20, 50, 100, 200, 500)
            ats <- labels
            if (logY) ats <- log10(ats)
            axis(2, at=ats, labels=labels, las=1, xpd=FALSE,cex.axis=0.8)
   
            labels <- c(0.5, 1, 5, 10, 50, 100, 200, 500)

            ats <- c(0.1, labels)
            labels <- c("mock\ncontrol", labels)

            if (logX) ats <- log10(ats)
            axis(1, at=ats, labels=labels, las=1, xpd=FALSE,cex.axis=0.8,
                 mgp=c(3,0.5,0))
   
               legend("bottomright",
                     legend=c(levels(factor(xy$protocolNames))),
                     text.col=c("darkOrange3","steelBlue3"),
                     pch=c(NA,NA,20, 17),
                     bty="n")
   
            dev.off()
         }



         plotX <- xy$yield
         plotY <- xy$expr

         if (logX) plotX <- log10(1 + plotX)
         if (logY) plotY <- log10(1 + plotY)

         plot(plotX,
              plotY,
              xlab = "", ylab="",
              main="",
              pch          = 20,
              cex          = 0.75,
              las          = 1,
              col          = pointCols,
              yaxt         = "n",
              xaxt         = "n")

         mtext(side = 1, text = "nuclear yield (thousands)", line = 1.5)
         mtext(side = 2, text = paste0(geneName, " abundance (TPM+1)"),
               cex=0.9,
               line = 2.7)

         yRange <- range(xy$expr)

         labels <- c(1, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000)
         ats <- labels
         if (logY) ats <- log10(ats)
         axis(2, at=ats, labels=labels, las=1, xpd=FALSE,cex.axis=0.8,
              mgp=c(3,0.5,0))

         labels <- c(1, 2, 5, 10, 20, 50, 100, 200, 500)
         ats <- labels
         if (logX) ats <- log10(ats)

         axis(1, at=ats, labels=labels, las=1, xpd=FALSE,cex.axis=0.8,
              mgp=c(3,0.5,0))

         if (addLabels) {
            library(calibrate)
            textxy(plotX,
                   plotY,
                   sub("-rep.*", "", rownames(xy)),
                   cex     = 0.5,
                   offset  = 0,
                   col     = "black")
         }

         if (runFit) {
            for (j in c(1:eyFit$numComp)) {
               print(paste0("****** NOW ON COMPONENT: ",j))
               curCol <- "black"
               if (eyFit[["bgComp"]] == j) {
                  curCol <- "red"
               }

               print(paste0("comp ", j, "yint = ", eyFit$curves[[j]]["yint"]))
               print(paste0(" (pval) = ", eyFit$curves[[j]]["pval.yint"]))
               print(paste0("-> slope = ", eyFit$curves[[j]]["slope"]))

               abline(a   = eyFit$curves[[j]]["yint"],
                      b   = eyFit$curves[[j]]["slope"],
                      col = curCol,
                      lwd = 2,
                      xpd = FALSE)
            }

            if (multiPlot) {
               legendText <- paste0(i,".",geneName)
            } else {
               legendText <- geneName
            }
            if (eyFit[["bgComp"]] > 0) {
               legendText <- paste(legendText,
                                   paste0("slope: ",sprintf("%.2f", eyFit$curves[[eyFit[["bgComp"]]]]["slope"])),
                                   paste0("est BG: ",sprintf("%.2f",
eyFit[["bgEst.yield0"]]), " TPM"),
                                   sep="\n")
            }
            legend("topright",legend=legendText,bty="n",text.col="red",xjust=1)
         }

         if (multiPlot) {
            legend("topright", inset=c(-0.6,0),
                  legend=levels(factor(xy$protocolNames)),
                  text.col=c("darkOrange3","steelBlue3"),
                  title.col="black",
                  title="Protocols",bty="n")
         } else {
            legend("bottomleft",
                  legend=levels(factor(xy$protocolNames)),
                  text.col=c("darkOrange3","steelBlue3"),
                  bty="n")
         }
         dev.off()
      }
   }
   if (multiPlot) dev.off()

   return(TRUE)

}


fitExprYield <- function( expr,
                          yield,
                          geneName      = NULL,
                          decontamMode  = "binary",
                          debug         = TRUE,
                          details       = FALSE,
                          includeBgLevel= FALSE){

# Purpose: Fit linear regression model to expression~yield relationship

   library(flexmix)

   out <- list(geneName = geneName)

   flexmix.nrep <- 10

# order    
   y <- cbind(expr=expr, yield=yield)
   y <- as.data.frame(log10(1 + y))

# FlexMix to fit 1- and 2-component glm's to log10(1 + expr)~log10(1 + yield)


# jitter'ing to avoid errors for genes like CG43134

   fitEY1 <- stepFlexmix(jitter(y$expr) ~ y$yield, k=1, nrep=flexmix.nrep)
   icl1 <- ICL(fitEY1)

   fitEY2 <- stepFlexmix(jitter(y$expr) ~ y$yield, k=2, nrep=flexmix.nrep)
   icl2 <- ICL(fitEY2)

   if (debug) {
      print(paste0("ICL: 1-comp=",icl1,"; 2-comp=",icl2))

      print("Summary FIT1:-----------------")
      print(summary(fitEY1))

      print("Summary FIT2:-----------------")
      print(summary(fitEY2))
   }

   numComp <- 1
   bgComp <- 1 #initialize to component 1

   if( icl2 < icl1 &
      length(y$expr[clusters(fitEY2) == 1]) > 0 &
      length(y$expr[clusters(fitEY2) == 2]) > 0) {

      numComp <- 2

# Figure out which component may be the background one;
#
# * Picking the lower y-intercept doesn't always work:
#      fitYE2 <- relabel(fitEY2,by="model",which="Intercept")
#
# * Pick the one with the lower mean expression of assigned points

      exprDist1 <- y$expr[clusters(fitEY2) == 1]
      exprDist2 <- y$expr[clusters(fitEY2) == 2]

      if (debug) {
         print(paste0("gene ",geneName, " 2comp fit:"))
         print(paste0("gene ",geneName, " 2comp fit:"))
         print(paste("cluster1: ",exprDist1))
         print(paste("cluster2: ",exprDist2))

         print(paste0("component 1: n=",length(exprDist1)))
         print(paste0("component 2: n=",length(exprDist2)))
      }


# Note for some genes, eg CG12107, this throws:
#    Warning in sqrt(diag(z@vcov)[indices]) : NaNs produced

      fitEY <- fitEY2
      refitEY <- refit(fitEY2)
      out[["curves"]] <- c(
         list(c(
            yint = attr(refitEY,"components")[[1]]$Comp.1[1,1],
            slope = attr(refitEY,"components")[[1]]$Comp.1[2,1],
            pval.yint = attr(refitEY,"components")[[1]]$Comp.1[1,4],
            pval.slope = attr(refitEY,"components")[[1]]$Comp.1[2,4],
            n=length(exprDist1)
         )),
         list(c(
            yint = attr(refitEY,"components")[[1]]$Comp.2[1,1],
            slope = attr(refitEY,"components")[[1]]$Comp.2[2,1],
            pval.yint = attr(refitEY,"components")[[1]]$Comp.2[1,4],
            pval.slope = attr(refitEY,"components")[[1]]$Comp.2[2,4],
            n=length(exprDist2)
         ))
      )


      kspval <- ks.test(exprDist1, exprDist2, alternative="less")$p.value
      if (kspval < 0.01 | max(exprDist2) < max(exprDist1)) {
         bgComp <- 2
      }
      if (debug) {
         print(paste0("max exprDist1 = ",max(exprDist1)))
         print(paste0("max exprDist2 = ",max(exprDist2)))
         print(paste0("Found 2 components; component ", bgComp, " has lower expr"))
      }

   } else {

      fitEY <- fitEY1
      refitEY <- refit(fitEY1)

      if (debug) print("Found 1 component")

      out[["curves"]] <- c(
         list(c(
            yint = attr(refitEY,"components")[[1]]$Comp.1[1,1],
            slope = attr(refitEY,"components")[[1]]$Comp.1[2,1],
            pval.yint = attr(refitEY,"components")[[1]]$Comp.1[1,4],
            pval.slope = attr(refitEY,"components")[[1]]$Comp.1[2,4]
         ))
      )

   }

   print("PREDICTED VALUES FROM: ")
   print(cbind(matrix(unlist(predict(fitEY)),ncol=numComp,byrow=FALSE), y$expr))
   print("---------------")


   out[["numComp"]] <- numComp

   bgComp.num <- bgComp
   out[["bgComp"]] <- bgComp
   bgComp <- paste0("Comp.",bgComp)

# Determine if p-value of fit is significant
   pval.intercept <- attr(refitEY,"components")[[1]][[bgComp]][1,4]
   pval.slope <- attr(refitEY,"components")[[1]][[bgComp]][2,4]

   bgComp.intercept <- attr(refitEY,"components")[[1]][[bgComp]][1,1]
   bgComp.slope <- attr(refitEY,"components")[[1]][[bgComp]][2,1]

   comp1.intercept <- attr(refitEY,"components")[[1]]$Comp.1[1,1]
   comp1.slope <- attr(refitEY,"components")[[1]]$Comp.1[2,1]

   if (debug) {
      print(paste0("-> Component 1 slope: ",comp1.slope," y-int: ", comp1.intercept))
   }

   if (numComp == 2) {
      comp2.intercept <- attr(refitEY,"components")[[1]]$Comp.2[1,1]
      comp2.slope <- attr(refitEY,"components")[[1]]$Comp.2[2,1]
      if (debug) {
         print(paste0("-> Component 2 slope: ",comp2.slope," y-int: ", comp2.intercept))
      }
   }

# If slope is sig negative (< -0.3), call it bg

   if (debug) {
print(paste0("bgComp.slope=",bgComp.slope))
print(paste0("pval.intercept=",pval.intercept))
print(paste0("pval.sope=",pval.slope))
   }

   isBG <- FALSE
   if (!is.na(pval.intercept) & pval.intercept < 0.01 &
       !is.na(pval.slope) & pval.slope < 0.05 &
       !is.na(bgComp.slope) & bgComp.slope < 0) {

      isBG <- TRUE
      bgEst.yield0 <- c(10 ^ bgComp.intercept - 1)
      if (debug) {
         print(paste0("-> bgEst = ",bgEst.yield0))
      }
      out[["bgEst.yield0"]] <- bgEst.yield0


      bgEst <- c(10 ^ (bgComp.slope * y$yield + bgComp.intercept) - 1)

      decontamExpr <- c(expr - bgEst)
      decontamExpr[decontamExpr < 0] <- 0

#      if (decontamMode == "binary") {
#         decontamExpr[clusters(fitEY) == bgComp.num] <- 0
#      }


      if (debug) {
            print(paste0("on ", geneName))
            print(paste0("bgEst.yield0 = ",bgEst.yield0))
            print(paste("expr = ",expr))
            print(paste("bgEst = ",bgEst))
            print(paste("decontamExpr = ",decontamExpr))
      }


      out[["decontamExpr"]] <- decontamExpr


   } else {
      if (debug) {
         print(paste0("no bg found for ", geneName))
      }

      out[["bgComp"]] <- 0
      out[["bgEst.yield0"]] <- 0
      out[["decontamExpr"]] <- expr
   }


   if (includeBgLevel) {
      out[["decontamExpr"]] <- c(out[["decontamExpr"]], out[["bgEst.yield0"]])
   }



##   print(paste0("~~~~~~~~ BG COMP = ", out[["bgComp"]]))
   if (details) {
      return(out)
   } else {
      return(out[["decontamExpr"]])
   }

#   plot(y,col=c("red","black")[clusters(testfit)])

# constraints:
# 1. reps should be assigned same on/off value

}


plot.ERCC.calibration.curves<- function( data,
                                         figName,
                                         returnMat = FALSE,
                                         returnMols = FALSE,
                                         returnSampleInfo = FALSE) {

# Purpose: plot ERCC calibration curves

   erccSpecs <- list()
   erccSpecs$atto <- 1E-18
   erccSpecs$micro<- 1E-6          #ERCC conc in am/uL
   erccSpecs$Na   <- 6.0221409e+23 #mol / mole

   erccSpecs$vol  <- 1E-6   # 1 uL
   erccSpecs$dil  <- 1E-4   # 1:10,000


   erccInfo <- read.table(data$specs$erccConcFn,
                          sep    = "\t",
                          header = TRUE,
                          col.names = c( "Resort.ID",
                                         "ERCC.ID",
                                         "subgroup",
                                         "conc.mix1",
                                         "conc.mix2",
                                         "exp.fc.ratio",
                                         "log2.mix1.mix2" )
                          )


   allSamples <- data$specs$sampleInfo[data$specs$sampleInfo$intact.protocol %in%
                                       c("INTACT", "TAPIN") &
                                       data$specs$sampleInfo$driverType != "control",]
   tpmCols <- paste0("tpm.", allSamples$sample_name)
   print(paste0("num cols = ", length(tpmCols))) ;
   rownames(allSamples) <- tpmCols
   tpnCols <- c()


   erccTPM <- merge(data$expr$geneExpr.withERCC,
                    erccInfo[, c("ERCC.ID", "conc.mix1", "conc.mix2")],
                    by.x  = "gene_name",
                    by.y  = "ERCC.ID",
                    all.x = FALSE,
                    all.y = TRUE)


   tpmRange<-c()
   tpnRange<-c()
   erccTpmFits<-list()
   erccTpnFits<-list()
   withTpn<-c()
   allSamples$numEndogMols <- 0
   allSamples$numAllMols <- 0
   tpn.factors <- vector()
   erccCol <- "conc.mix1" #all samples used mix 1
   for (i in (1:nrow(allSamples))) {
      tpmCol <- tpmCols[i]

      erccTPM[,tpmCol] <- as.numeric(erccTPM[,tpmCol])

      tpnCol <- paste0("tpn.", allSamples$sample_name[i])

      tpmRange   <- range(c(tpmRange, erccTPM[, tpmCol]))

# calc best fit line ercc conc -- tpm table.
      erccTpmFits[[i]] <- lm(
            log10( 1 + erccTPM[erccTPM[,tpmCol] >= 1 ,tpmCol]) ~
            log10( erccTPM[erccTPM[,tpmCol] >= 1 ,erccCol]) )

# convert ercc conc to tx / nuclei if nuclei yield available and more than 0
      if (is.na(allSamples$yield.nuclei.thousands[i]) |
          allSamples$yield.nuclei.thousands[i] == "unk" |
          allSamples$yield.nuclei.thousands[i] == 0) {
         next
      }

      tpnCols<-c(tpnCols,tpnCol)
      withTpn<-c(withTpn, i)

      numNuclei <- as.numeric(allSamples$yield.nuclei.thousands[i]) * 1E3
#      print(paste0("nuclei yield of ",tpmCol," is ",numNuclei))

      tpn.factors[tpmCol] <- c( (erccSpecs$atto / erccSpecs$micro) * erccSpecs$vol *
                             erccSpecs$dil * erccSpecs$Na ) / numNuclei

      fracERCCtpm <- sum(erccTPM[,tpmCol]) /
                     sum(data$expr$geneExpr.withERCC[,tpmCol])

      allSamples[tpmCol,"numEndogMols"] <- sum(erccTPM[,erccCol]) *
         (erccSpecs$atto / erccSpecs$micro) * erccSpecs$vol *
         erccSpecs$dil * erccSpecs$Na * (1 - fracERCCtpm) / fracERCCtpm

      allSamples[tpmCol,"numAllMols"] <- sum(erccTPM[,erccCol]) *
         (erccSpecs$atto / erccSpecs$micro) * erccSpecs$vol *
         erccSpecs$dil * erccSpecs$Na / fracERCCtpm

      erccTPM[,tpnCol] <- erccTPM[,erccCol] * tpn.factors[tpmCol]

      erccTpnFits[[i]] <- lm(
            log10( 1 + erccTPM[erccTPM[,tpmCol] >= 1 ,tpmCol]) ~
            log10( erccTPM[erccTPM[,tpmCol] >= 1 , tpnCol]) )

      tpnRange <- range(c(tpnRange,
                          erccTPM[ erccTPM[,tpmCol] >= 1, tpnCol]))
   }
   if (returnMat) { return(erccTPM) }


   allSamples$shortName <- rownames(allSamples)
   allSamples$shortName <- gsub("tpm.","", allSamples$shortName)
   if (returnMols) {return(allSamples)}

   allSamples$numDetGenes <- 0
   for (tpmCol in tpmCols) {
      allSamples[tpmCol,"numDetGenes"] <- sum(exp(data$expr$geneExpr[,tpmCol]) > 1) }

   if (returnSampleInfo) {return(allSamples)}

   outPdf <- paste0(data$specs$outDir, "/",figName,"_nummols_library_stats.pdf")
   pdf(outPdf,height=15,width=10)
   par(mfrow=c(3,2))

   plot(allSamples$yield.cdna.ng_uL, allSamples$numDetGenes,
     pch=20, cex=0.5, log="x", xlab="cDNA yield",
     ylab="Number of detected genes (TPM > 0)", main="")

   text(allSamples$yield.cdna.ng_uL, allSamples$numDetGenes,
     allSamples$shortName, pos=4, cex=0.3, col="darkgray")

   points(allSamples$yield.cdna.ng_uL, allSamples$numDetGenes, pch=20, cex=0.5)


   plot(allSamples$yield.nuclei.thousands, allSamples$numDetGenes,
     pch=20, cex=0.5, log="x", xlab="nuclear yield (thousands)",
     ylab="Number of detected genes (TPM > 0)", main="")

   text(allSamples$yield.nuclei.thousands, allSamples$numDetGenes,
     allSamples$shortName, pos=4, cex=0.3, col="darkgray")

   points(allSamples$yield.nuclei.thousands, allSamples$numDetGenes, pch=20, cex=0.5)


   plot(allSamples$yield.cdna.ng_uL, allSamples$numEndogMols,
     pch=20, cex=0.5, log="xy", xlab="cDNA yield",
     ylab="Transcript molecules (endogenous)", main="")

   text(allSamples$yield.cdna.ng_uL, allSamples$numEndogMols,
     allSamples$shortName, pos=4, cex=0.3, col="darkgray")

   points(allSamples$yield.cdna.ng_uL, allSamples$numEndogMols, pch=20, cex=0.5)


   plot(allSamples$yield.nuclei.thousands, allSamples$numEndogMols,
     pch=20, cex=0.5, log="xy", xlab="nuclear yield (thousands)",
     ylab="Transcript molecules (endogenous)", main="")

   text(allSamples$yield.nuclei.thousands, allSamples$numEndogMols,
     allSamples$shortName, pos=4, cex=0.3, col="darkgray")

   points(allSamples$yield.nuclei.thousands, allSamples$numEndogMols, pch=20, cex=0.5)


   plot(allSamples$yield.cdna.ng_uL, allSamples$numAllMols,
     pch=20, cex=0.5, log="xy", xlab="cDNA yield",
     ylab="Transcript molecules (endogenous + spike)", main="")

   text(allSamples$yield.cdna.ng_uL, allSamples$numAllMols,
     allSamples$shortName, pos=4, cex=0.3, col="darkgray")

   points(allSamples$yield.cdna.ng_uL, allSamples$numAllMols, pch=20, cex=0.5)


   plot(allSamples$yield.nuclei.thousands, allSamples$numAllMols,
     pch=20, cex=0.5, log="xy", xlab="nuclear yield (thousands)",
     ylab="Transcript molecules (endogenous + spike)", main="")

   text(allSamples$yield.nuclei.thousands, allSamples$numAllMols,
     allSamples$shortName, pos=4, cex=0.3, col="darkgray")

   points(allSamples$yield.nuclei.thousands, allSamples$numAllMols, pch=20, cex=0.5)


   dev.off()


# AALL DONE ---------------------------------------------------------


# Make plot overlaying all best fit lines (alpha) onto a single plot
# plot 1. TPM vs mix conc

   outPdf<-paste0(data$specs$outDir, "/",figName,"_lineplots_ercc_calibration.pdf")
   pdf(outPdf)

   plot( log10(erccTPM[,"conc.mix1"]),
         log10(1 + erccTPM[,tpmCols[1]]),
         xlab = "ERCC original concentration (log10 aM/uL)",
         ylab = "estimated abundance (log10 TPM + 1)",
         type = "n",
         ylim = log10(1 + tpmRange) )

   for (i in c(1:nrow(allSamples))){
      abline(erccTpmFits[[i]],
             lwd = 2,
             col = rgb(0,0,0,0.25))
   }


   shortNames <- gsub("ERCC-0*","",erccTPM$gene_name)
   print(shortNames)
   for (i in c(1:nrow(allSamples))){
      plot( jitter(log10(erccTPM[,"conc.mix1"])),
            log10(1 + erccTPM[,tpmCols[i]]),
            xlab = "ERCC original concentration (aM)",
            ylab = "estimated abundance (TPM + 1)",
            main = tpmCols[i],
            ylim = log10(1 + tpmRange),
            pch=20,
            type="n")
      textxy( jitter(log10(erccTPM[,"conc.mix1"])),
              log10(1 + erccTPM[,tpmCols[i]]),
              shortNames,
              cex=0.7)
      abline(erccTpmFits[[i]],
             lwd = 2,
             col = rgb(0,0,0))
   }

   dev.off()


# plot 2. TPM vs TPN

   outPdf<-paste0(data$specs$outDir, "/",figName,"_lineplots_ercc_tpn.pdf")
   pdf(outPdf,height=2,width=2.5)
   par(mar=c(3,4.0,0.5,0.5),ps=10,xpd=FALSE,pty="s")

   plot( log10(erccTPM[,tpnCols[1]]),
         log10(1 + erccTPM[,tpmCols[1]]),
         xlab = "", ylab="",
         type = "n",
         xlim = c(-2,6),
         ylim = c(0,6),
         ann  = FALSE,
         axes = FALSE)

   mtext(side = 1, text = "true abundance\n(transcripts per nuclei)", line = 2.0)
   mtext(side = 2, text = "relative abundance\n(TPM + 1)", line = 2.8)


   axis(1,
        at     = c(-2:6),
        labels = c(10^c(-2:6)),
        cex.axis = 0.8,
        mgp=c(3,0.5,0),
        las    = 1)

   axis(2,
        at     = c(0:6),
        labels = c(10^c(0:6)),
         cex.axis = 0.8,
        las    = 1)

   for (i in withTpn){
      abline(erccTpnFits[[i]],
             lwd = 1,
             col = rgb(0, 0, 0, 0.25))
   }
   dev.off()

}


plot.readCoverage <- function( data, figName, type,
                              cdna.conc_to_ug = 0.03) {

# Purpose: Parse PICARD reports of normalized read coverage and plot individual
#          lines and smoothScatter overlay.


   allSamples <- data$specs$sampleInfo[!is.na(data$specs$sampleInfo$yield.cdna.ng_uL),]
   allSamples <- allSamples[allSamples$driverType != "control",]

   rankLib <- rank(-1 * allSamples$yield.cdna.ng_uL)
   colLib <- rgb(rankLib/length(rankLib),0,0,0.4)

   readCov<-list()
   yRange<-c()
   maxBias <- vector()
   print("Reading PICARD coverage reports")
   for (i in (1:nrow(allSamples))) {
      picardFn<-paste0(data$specs$picardDir,"/",
                       allSamples$project[i],"/",
                       allSamples$sampleID[i],"/",
                       allSamples$sampleID[i],
                       ".noERCC.picard_rnaseq_report.txt")

      readCov[[i]]<-read.table(pipe(paste0("tail -103 ",picardFn)),
                               sep       = "\t",
                               header    = TRUE,
                               col.names = c("pos","read.coverage"))
      yRange<-range(c(yRange, readCov[[i]]$read.coverage))
      maxBias[i] <- max(readCov[[i]]$read.coverage)
   }

   if (type == "txPosition") {
   print("Plotting read coverage summaries")
   outPdf<-paste0(data$specs$outDir,"/",figName,"_lineplots_picard_readcoverage.pdf")
   pdf(outPdf,height=2,width=2.5)
   par(mar=c(3,4.0,0.5,0.5),ps=10,xpd=TRUE)
   plot(readCov[[1]],
        type    = "n",
        main    = "",
        ylim    = yRange,
        xaxt = "n", yaxt="n",
        xlab = "", ylab="")

   mtext(side = 1, text = "Gene body position (%)", line = 1.5)
   mtext(side = 2, text ="Coverage bias", line = 2)
   axis(2, las=1, xpd=FALSE,cex.axis=0.8)
   axis(1, las=1, xpd=FALSE,cex.axis=0.8, mgp=c(3,0.5,0))

   for (i in (1:nrow(allSamples))) {
      lines(readCov[[i]],
            col         = rgb(0,0,0,0.25),
#            col         = colLib[i],
            lwd         = 1)
   }

#   legend("top",legend=c("sample with lowest cDNA yield",
#                         "sample with highest cDNA yield"),
#       text.col=c("red", "black"), bty="n")

   dev.off()
   }

   if (type == "maxVsYield") {

   outPdf<-paste0(data$specs$outDir,"/",
                  figName,"_genecoveragebias_vs_cdnayield.pdf")
   pdf(outPdf, height=2,width=2.5)
   par(mar=c(3,4.0,0.5,0.5),ps=10,xpd=TRUE)
   plot(allSamples$yield.cdna.ng_uL * cdna.conc_to_ug,
        maxBias,
        pch=20,
        cex=0.8,
        xaxt = "n", yaxt="n",
        xlab = "", ylab="",
        main="")

   mtext(side = 1, text = "cDNA yield (ug)", line = 1.5)
   mtext(side = 2, text ="Maximum coverage bias", line = 2)
   axis(2, las=1, xpd=FALSE,cex.axis=0.8)
   axis(1, las=1, xpd=FALSE,cex.axis=0.8, mgp=c(3,0.5,0))

   dev.off()

   }

}



plotHeatmap <- function( data,
                         tpmMat,
                         as.is = FALSE,
                         deGenes,
                         plotWidth = 7,
                         plotHeight = 10,
                         treeheight_row=NULL,
                         treeheight_col=NULL,
                         fontsize_row = 4,
                         fontsize_col = 4,
                         main,
                         legend,
                         cluster_cols = TRUE,
                         cluster_rows = TRUE,
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         outFormat = "pdf", # "png"
                         goTerms,
                         geneClass,
                         geneLists=NULL, #optional list of genes
                         filePrefix){
# Purpose: plot expression heatmaps

# If no gene_id column, make one, assume rownames are names
   if ((! "gene_id" %in% colnames(tpmMat)) & !as.is) {
      tpmMat$gene_name <- rownames(tpmMat)
      tpmMat <- merge(tpmMat,
         unique(data$expr$transcriptInfo[, c( "gene_id","gene_name")]))
   }

   library(pheatmap)

   minMaxExprThresh <- 30
   minFc <- 4

   colorScaleFull <- colorRampPalette(c("steelBlue3",
                                        "white",
                                        "darkOrange3"))(100)

   tpmCols <- colnames(tpmMat)[grep("tpm",colnames(tpmMat))]

   geneid2name <- unique(data$expr$transcriptInfo[,c("gene_id","gene_name")])
   rownames(geneid2name)<-geneid2name$gene_name


   if (missing(goTerms) & !as.is) {
      goTerms <- c()
      if (!is.null(geneLists)) {
         goTerms <- names(geneLists)
      } else {
         goTerms <- c("neuropeptide signaling pathway",
                      "G-protein coupled receptor signaling pathway",
                      "ion channel activity",
                      "sequence-specific DNA binding transcription factor activity",
                      "cell adhesion",
                      "neuron projection")
      }
   }

   if (as.is) {goTerms <- c("ASIS")}

   hmaps <- list()
   for (goTerm in goTerms) {

      print(paste0("Now on go ",goTerm))

# If gene type is a biotype, deal with that separately!

      clusterRows <- TRUE

      rowAnn <- NA
      rowAnnCol <- list()
      if (goTerm == "ASIS") {
         clusterRows <- FALSE
         clusterCols <- FALSE

      } else if (goTerm == "ntSystems") {

         allGenes <- geneid2name[geneLists$ntSystems$geneName, "gene_id"]
         clusterRows <- FALSE
         rowAnn <- geneLists$ntSystems[, c("nt", "direction")]
         rownames(rowAnn)<-geneLists$ntSystems$geneName
         rowAnn$nt<-factor(rowAnn$nt)
         rowAnn$direction<-factor(rowAnn$direction)
         rowAnnCol <- list(
            direction = c(Out = "black", In = "darkgrey")
         )

      } else if ( grepl("^IPR", goTerm) ) {

         allGenes <- data$specs$interproInfo$gene_id[
                        data$specs$interproInfo$domain == goTerm ]


      } else if ( !missing(geneLists) & (goTerm %in% names(geneLists))) {

         allGenes <- geneid2name[geneLists[[goTerm]], "gene_id"]
         print(paste0("Gene list for ",goTerm," is:"))
         print(paste0(allGenes, collapse=", "))

      } else if ( grepl("biotype", goTerm) ) {

         biotype  <- sub("biotype:", "", goTerm)
         allGenes <- unique(data$txExpr$gene_id[a$txExpr$gene_biotype == biotype])

      } else {

         allGenes <- data$specs$goInfo$gene_id[ data$specs$goInfo$go_term == goTerm ]

      }

      if (!as.is) {
      curExpr <- tpmMat[tpmMat$gene_id %in% allGenes,
                        c("gene_name", tpmCols)]
      curExpr <- tpmMat[tpmMat$gene_id %in% allGenes,
                        c("gene_name", tpmCols)]

      rownames(curExpr) <- curExpr$gene_name

      # reorder if ntSystems
      if (!is.na(rowAnn)){
         curExpr<-curExpr[geneLists$ntSystems$geneName,]
      }

      curExpr$gene_name <- NULL

      curExpr <- curExpr[ apply(curExpr, 1, max) >= minMaxExprThresh, ]
      print(paste0("   all genes (TPM > ", minMaxExprThresh, "): ", nrow(curExpr)))

      } else {
         print("SETTING curExpr to tpmMat")
         curExpr <- tpmMat
      }

      colnames(curExpr) <- gsub("tpm.", "", colnames(curExpr))

      curFn <- paste0("heatmap_tpm_allGenes_", goTerm, ".pdf")
      curFn <- gsub(" ", "_", curFn)

      if (!missing(filePrefix)) curFn <- paste0(filePrefix, ".", curFn)

      curExpr <- log2((1 + curExpr) / (1 + apply(curExpr, 1, mean)))

      curExpr[curExpr > 2.5] <- 2.5
      curExpr[curExpr < -2.5] <- -2.5
      curBreaks <- seq(-2.5, 2.5, length.out=101)

      curFn <- paste0("heatmap_relexpr_allGenes_",goTerm,".",outFormat)
      curFn <- gsub(" ", "_", curFn)

      if (!missing(filePrefix)) curFn <- paste0(filePrefix, "_", curFn)

      if (outFormat == "pdf") {
         pdf(paste0(data$specs$outDir, "/",  curFn),
             height=plotWidth, width=plotHeight,
            onefile=FALSE)
      } else {
         png(paste0(data$specs$outDir, "/",  curFn),
             width=plotWidth, height=plotHeight,
             units="in", res=300)
      }


      if (missing(cluster_rows)) { cluster_rows <- clusterRows }
      if (missing(cluster_cols)) { cluster_cols <- TRUE }

      if (is.na(rowAnn)) {

         curMain <- "relative gene expression (log2 TPM/mean)"
         curLegend <- TRUE
         if (!missing(main)) {curMain <- main}
         if (!missing(legend)) {curLegend <- legend}

         if (cluster_cols == FALSE) {treeheight_col <- NA}
         if (cluster_rows == FALSE) {treeheight_row <- NA}

         hmaps[[goTerm]] <- pheatmap( curExpr,
                scale           = "none",
               treeheight_row = treeheight_row,
               treeheight_col = treeheight_col,
                main            = curMain, legend = curLegend,
                border_color    = NA,
                col             = colorScaleFull,
                breaks          = curBreaks,
               cluster_cols = cluster_cols,
               cluster_rows = cluster_rows,
               show_rownames = show_rownames,
               show_colnames = show_colnames,
               fontsize_col     = fontsize_col,
               fontsize_row     = fontsize_row)
      } else {
         curMain <- "relative gene expression (log2 TPM/mean)"
         curLegend <- TRUE

         if (!missing(main)) {curMain <- main}
         if (!missing(legend)) {curLegend <- legend}

         print(paste0("CALLING PHEATMAP WITH cluster_cols = ",cluster_cols))
         print(paste0("CALLING PHEATMAP WITH cluster_rows = ",cluster_rows))
         if (cluster_cols == FALSE) {treeheight_col <- NA}
         if (cluster_rows == FALSE) {treeheight_row <- NA}

         hmaps[[goTerm]] <- pheatmap( curExpr,
                scale           = "none",
               treeheight_row = treeheight_row,
               treeheight_col = treeheight_col,
                main            = curMain, legend = curLegend,
                border_color    = NA,
               show_rownames = show_rownames,
               show_colnames = show_colnames,
               cluster_cols = cluster_cols,
               cluster_rows = cluster_rows,
               annotation_row = rowAnn[rownames(curExpr),],
               annotation_colors = rowAnnCol,
                col             = colorScaleFull,
                breaks          = curBreaks,
               fontsize_col     = fontsize_col,
               fontsize_row  = fontsize_row)
      }
      dev.off()

   }

   return(hmaps)

}



plotHmapPmat <- function(data, mode="gs",
                         exprMode = "p", #or raw
                         tpmOverlay = FALSE, #write TPM's in heatmap cells
                         rawNL = "mean", #or mock or Z
                         rawSamples, #if exprMode=raw, optional specify samples
                         suppMat, #non-gene features (eg, pheno) to add to matrix
                         suppColors,
                         geneClass,
                         geneList, 
                         flybaseGeneGroup,
                         colType = "cell", #or "gene"
                         cellGroups, #optional list of cell groups
                         geneLabel = FALSE,
                         lip=FALSE,
                         gaps_genes = NULL,
                         gaps_genes.group_labels,
                         gaps_cells = NULL,
                         gaps_col = NULL,
                         gaps_row = NULL,
                         minPon = 0,
                         fontsize_col=7,
                         fontsize_row=5,
                         treeheight_row=NULL,
                         treeheight_col=NULL,
                         legend,
                         plotTitle,
                         plotWidth=4,
                         plotHeight=6.5,
                         refHmap, #optional pheatmap object for row/col order
                         cluster_cols = TRUE,
                         cluster_rows = TRUE,
                         outFormat = "pdf",  # or "png"
                         figName) {

# Purpose: make heatmap of probability (or raw) expression matrix

   library(pheatmap) 

   if (exprMode == "p") {
      hmapBreaks <- seq(0,1,length.out=101)
      hmapCols <- colorRampPalette(c("steelBlue3","white","darkOrange3"))(100)

      if (mode == "gs") {
         pMat <- exp(data$pFit$logp_on_gs[,
                        setdiff(colnames(data$pFit$logp_on_gs),
                                c("tpm.min","tpm.max"))]) ;
      } else if (mode == "gd") {
         pMat <- exp(data$pFit$logp_on_gd[,
                        setdiff(colnames(data$pFit$logp_on_gd),
                                c("min","max"))]) ;
      } else if (mode == "gc") {
         pMat <- exp(data$pFit$logp_on_gc[,
                        setdiff(colnames(data$pFit$logp_on_gc),
                                c("min","max"))]) ;
      }
   
   } else if (exprMode == "raw") {

      if (missing(rawSamples)) {
         rawSamples <- sort(data$specs$sampleLists$good$allCriteria)
      }

      t.rawMat <- data$expr$geneExpr[, paste0("tpm.", rawSamples)]
      rownames(t.rawMat) <- data$expr$geneExpr$gene_name

      resMode <- "sample"
      if (mode == "gd") {
         resMode <- "driver"
      } else if (mode == "gc") {
         resMode <- "celltype"
      }

      if (resMode != "sample") {

         t.curUnits <- sort(unique(data$specs$sampleInfo[
                                       data$specs$sampleInfo$sample_name %in%
                                       rawSamples, resMode]))

         tMat <- matrix( nrow = nrow(t.rawMat),
                         ncol = length(t.curUnits))

         rownames(tMat) <- rownames(t.rawMat)
         colnames(tMat) <- t.curUnits

         for (t.curUnit in t.curUnits) {
            print(paste0("Calculating average for ",resMode," ", t.curUnit))

            t.curSamples <- data$specs$sampleInfo$sample_name[
                              data$specs$sampleInfo[[resMode]] %in%
                              c(t.curUnit)]

            t.curSamples <- intersect(t.curSamples, rawSamples)
            t.curSamples <- paste0("tpm.",t.curSamples)
            tMat[,t.curUnit] <- apply(t.rawMat[,t.curSamples],1,mean)
         }
         tMat <- as.data.frame(tMat)

         t.rawMat <- tMat
      }

      pMat <- t.rawMat
      origRawMat <- t.rawMat

      if (rawNL %in% c("mean","Z","mock")) {
         curColors <- c("steelBlue3","white","darkOrange3")
         colorScaleMin <- -2.5
         colorScaleMax <- 2.5
      } else if (rawNL == "fracmax") {
         curColors <- c("white","darkOrange3")
         colorScaleMin <- 0
         colorScaleMax <- 1
      }
      hmapBreaks <- seq(colorScaleMin,colorScaleMax,length.out=101)
      hmapCols <- colorRampPalette(curColors)(100)

   } else {

      print("ERROR: unrecognized exprMode: either raw or p")
      return(1) ;

   }

   if (!missing(suppColors)) {
      suppBreaks <- 1001:(1000 + length(suppColors))
      suppCols   <- c()

      newSuppMat <- suppMat
      for (i in 1:length(suppColors)) {
         suppVal <- names(suppColors)[i]
         suppCols <- c(suppCols, suppColors[[suppVal]])
         newSuppMat[suppMat == suppVal] <- 1000 + i
      }
      suppMat <- data.matrix(as.data.frame(newSuppMat,as.is=TRUE,stringsAsFactors=FALSE))
      colnames(suppMat) <- colnames(newSuppMat)
      rownames(suppMat) <- rownames(newSuppMat)
      print(suppMat)

      hmapBreaks <- c(hmapBreaks, suppBreaks)
      hmapCols <- c(hmapCols, suppCols)

   }


   if (!missing(suppMat)) {

      pMat <- pMat[,colnames(suppMat)]

      suppMat <- suppMat[,colnames(pMat)]

      print("adding new rows:")
      print(rownames(suppMat))

      print(paste0("old Nrows= ", nrow(pMat)))
      pMat <- rbind(pMat, suppMat)
      print(paste0("new Nrows= ", nrow(pMat)))

   }

   outPre <- ""
   if (!missing(figName)) { outPre <- paste0(figName,".") }


   if (!missing(flybaseGeneGroup)) {
      tx <- read.table(dat$specs$flyBaseGroupsFn[[flybaseGeneGroup]],
                       header=FALSE, sep="\t")
      colnames(tx) <- c("geneID", "currentGeneName")

      geneList <- dat$specs$transcriptInfo$gene_name[
         match(tx$geneID, dat$specs$transcriptInfo$gene_id)]
   }



   if (!missing(geneList)) {

      if (length(setdiff(geneList, rownames(pMat))) > 0) {
         print("MISSING GENES: ")
         print(paste0(sort(setdiff(geneList, rownames(pMat))),collapse=", "))
      }

      pMat <- pMat[geneList,]
      print(paste0("HEYO NUMBER OF GENE: ",nrow(pMat)))

   } else if (!missing(geneClass)) {

      if (geneClass == "TF") {
         tfList <- loadTFs(data)
         tfList <- unique(data$specs$transcriptInfo$gene_name[
                           data$specs$transcriptInfo$gene_id %in% tfList])
         pMat <- pMat[rownames(pMat) %in% tfList,]
         if (nrow(pMat) < length(tfList)) {
            t1 <- setdiff(tfList,rownames(pMat))
            print(paste0("WARNING: missing ",length(t1)," TF genes in pMat: "))
            print(paste(sort(t1), sep=", "))
         }
      }
   }


   geneList <- rownames(pMat)

   if (exprMode == "p" & missing(refHmap)) { #skip if reference heatmap provided
      pMat <- pMat[apply(pMat,1,max) >= minPon,]
   }

   print(paste0("GENELIST SIZE: ",length(geneList)))
   print(paste0("-> vs NUMBER OF GENES IN HEATMAP: ",nrow(pMat)))

   if (outFormat == "png") {
      png(paste0(data$specs$outDir,"/",outPre,"pmat_hmap.png"),
         width=plotWidth, height=plotHeight, units="in", res=300)
   } else {
      pdf(paste0(data$specs$outDir,"/",outPre,"pmat_hmap.pdf"),
         width=plotWidth,height=plotHeight, onefile=FALSE)
   }
   if (!missing(cellGroups)) {
      if (colType == "gene") {
         if (missing(cluster_rows)) { cluster_rows <- FALSE }
      } else {
         if (missing(cluster_cols)) { cluster_cols <- FALSE }
      }
      gaps_cells <- c()
      cellOrder <- c()
      for (i in 1:length(cellGroups)){
         cellOrder <- unlist(c(cellOrder, cellGroups[i]))

         if (i < length(cellGroups)) {
            gaps_cells <- c(gaps_cells, length(cellOrder)) }
      }
      if (length(setdiff(cellOrder, colnames(pMat))) > 0) {
         print("CANT FIND:")
         print(setdiff(cellOrder, colnames(pMat)))
      }
      
      pMat <- pMat[, cellOrder]
   }


   if (exprMode == "raw") {
      if (rawNL == "fracmax") {
         pMat <- pMat / apply(pMat,1,max)
      } else if (rawNL == "mean") {


         pMat <- log2(1 + pMat)
         pMat <- pMat - apply(pMat,1,mean)

      } else if (rawNL == "Z") {


         pMat <- (pMat - apply(pMat,1,mean)) / apply(pMat,1,sd)

      } else if (rawNL == "mock") {


         pMat <- log2(1 + pMat)
         mockSamples <-  c("mock_201611_rep1", "mock_201611_rep2")

         mockProfile <- apply(data$expr$geneExpr[,paste0("tpm.",mockSamples)], 1, mean)
         mockProfile <- log2(1 + mockProfile)

         pMat <- pMat - mockProfile

      }

      pMat <- pMat
      pMat[pMat < colorScaleMin] <- colorScaleMin
      pMat[pMat > colorScaleMax] <- colorScaleMax

   }



   if (colType == "gene") {
      gaps_row <- gaps_cells
      pMat <- t(pMat)

      if (!missing(gaps_genes)) {gaps_col <- gaps_genes}

   } else {
      gaps_col <- gaps_cells
      if (!missing(gaps_genes)) {gaps_row <- gaps_genes}
   }


   if (!missing(refHmap)) { # expectes gene-x-sample matrix
      print("Reordering pMat to match reference heatmap")
      cluster_cols <- FALSE
      cluster_rows <- FALSE

      if (length(setdiff(refHmap$tree_col$labels,colnames(pMat))) > 0) {
         print("DIFFERENT COLUMN NAMES")
         if (any(grepl("^tpm.", refHmap$tree_col$labels))) {
            print("deletitng tpm. from refHmap col labels")
            refHmap$tree_col$labels <- gsub("^tpm.","",refHmap$tree_col$labels)
         } else {
            print("adding tpm. to refHmap col labels")
            refHmap$tree_col$labels <- paste0("tpm.",refHmap$tree_col$labels)
         }
         print(refHmap$tree_col$labels)
      }

      if (length(setdiff(refHmap$tree_row$labels,rownames(pMat))) > 0) {
         print("DIFFERENT ROW NAMES!")
         print(setdiff(refHmap$tree_row$labels,rownames(pMat)))
      }

      pMat <- pMat[refHmap$tree_row$labels[refHmap$tree_row$order],
                   refHmap$tree_col$labels[refHmap$tree_col$order]]
   }

   if (missing(plotTitle)) {
      if (exprMode == "p") {
         plotTitle <- paste0("P(on_",mode,")")
      } else if (exprMode == "raw") {
         plotTitle <- paste0("log2(TPM + 1) / mean")
      }
   }


   show_rownames <- TRUE
   show_colnames <- TRUE

   if (geneLabel == FALSE) {
      if (colType == "gene") {
         show_colnames <- FALSE
      } else {
         show_rownames <- FALSE
      }
   }

   if (fontsize_row == 0) {show_rownames <- FALSE}
   if (fontsize_col == 0) {show_colnames <- FALSE}

   curLegend <- TRUE
   if (!missing(legend)) {curLegend <- legend}

# HERENOW 180402_1121: get genes and sample/driver/cell of final matrix and
#                      if tpmOverlay is set, make a matching matrix and
#                      give to pheatmap as display_numbers


   displayNumbers <- FALSE
   if (tpmOverlay & missing(suppMat)) { #doesn't behave well with suppMat
      print("GOT HERE making tpmOverlay matrix")
      mat.geneNames <- colnames(pMat)
      mat.colNames <- rownames(pMat)
      print(paste0("gene Names: ", mat.geneNames))
      print(paste0("col Names: ", mat.colNames))

      if (mode == "gc") {
            rawponSet <- "cells"
      } else if (mode == "gd") {
            rawponSet <- "drivers"
      } else if (mode == "gs") {
            rawponSet <- "samples"
            mat.colNames <- paste0("tpm.",mat.colNames)
      }
      if (exprMode == "raw") {
         print("t.rawMat structure:")
         print(str(origRawMat))
         print("request mat.colNames not in the matrix:")
         print(paste0(setdiff(mat.colNames,colnames(origRawMat))))
         displayNumbers <- as.matrix(round(origRawMat[mat.colNames, mat.geneNames], digits=0))
      } else {
         displayNumbers <- as.matrix(t(dat$expr$rawpon[[rawponSet]]$raw[
                                          mat.geneNames, mat.colNames]))
      }

      displayNumbers <- round(displayNumbers,digits=0)
   }


   labels_col <- colnames(pMat)
   labels_row <- rownames(pMat)

   if (!show_colnames) {labels_col <- NULL}
   if (!show_rownames) {labels_row <- NULL}

   if (!missing(gaps_genes.group_labels)) {
      print(paste0("TRYING TO PRINT: ",gaps_genes.group_labels))
      new.labels <- rep("", ncol(pMat))
      new.labels[1] <- gaps_genes.group_labels[1]
      for (i in 1:length(gaps_genes)) {
         print(paste0("-> NOW ON: ", gaps_genes.group_labels[i + 1]))
         new.labels[1 + gaps_genes[i]] <- gaps_genes.group_labels[i + 1] }

      if (colType == "gene")  {
         show_colnames <- TRUE
         labels_col <- new.labels
      } else {
         show_rownames <- TRUE
         labels_row <- new.labels
      }
   }



   print("displayNumbres structure:")
   print(str(displayNumbers))
   print("pMat structure:")
   print(str(pMat))
   print("GOT HERE FINE!")
   hmap <- pheatmap(pMat,
            breaks = hmapBreaks,
            col = hmapCols,
            cluster_cols = cluster_cols,
            cluster_rows = cluster_rows,
            gaps_col = gaps_col,
            gaps_row = gaps_row,
            display_numbers = displayNumbers,
            number_color = "black",
            fontsize_number=c(fontsize_row - 1),
            fontsize_row=fontsize_row,
            fontsize_col=fontsize_col,
            treeheight_row = 0,
            treeheight_col = 0,
            border_color = NA,
            show_rownames = show_rownames,
            show_colnames = show_colnames,
            labels_row = labels_row,
            labels_col = labels_col,
            legend = curLegend,
            main = plotTitle)

   dev.off()

   return(hmap)

}

plotHmapPmat.nt_and_rec <- function( dat,
                                     figName1  = "msFig6A",
                                     figName2  = "msFig7A") {

# Purpose: make heatmap of neurotransmitter output and input genes
   
   outGenes <- c("brp", "nSyb", "Syt1", "Snap25","cpx",
                 "Hdc", "t", "CarT", "e","CG3790", 
                 "Gad1", "VGAT",
                 "VAChT", "ChAT",
                 "VGlut", "Eaat1",
                 "ple", "Vmat", "DAT",
                 "CG8468", "Eaat2", "prt", "Nos")
   outGenes.gaps <- c(5,10,12,14,16,19)

   inGenes <- dat$specs$ntSystems[
               dat$specs$ntSystems$geneName %in% rownames(dat$pFit$logp_on_gc) &
               dat$specs$ntSystems$direction == "In" &
               dat$specs$ntSystems$nt %in% c("ACh", "DA", "GABA", "Glu", "Hist", "Oct", "5HT"),]

   inGenes$nt[inGenes$nt == 'Hist'] <- "1.Hist"
   inGenes$nt[inGenes$nt == 'GABA'] <- "2.GABA"
   inGenes$nt[inGenes$nt == 'ACh'] <- "3.ACh"
   inGenes$nt[inGenes$nt == 'Glu'] <- "4.Glu"
   inGenes$nt[inGenes$nt == 'DA'] <- "5.DA"
   inGenes$nt[inGenes$nt == 'Oct'] <- "6.Oct"
   inGenes$nt[inGenes$nt == '5HT'] <- "7.5HT"

   inGenes <- inGenes[order(inGenes$nt),]
   gaps_genes <- c()
   for (i in 1:nrow(inGenes)) {
      if (i > 1) {
         if (inGenes[i,"nt"] != inGenes[(i - 1), "nt"]) {
            gaps_genes<- c(gaps_genes, (i - 1)) } } }

   inGenes <- inGenes$geneName
   print(gaps_genes)

   print("WARNING MISSING GENES BECAUSE OF NAME / CG / ID bullshit!")
   print(inGenes)
  
# Add special column denoting NT call
#   - color code by adding values outside 0-1 range
#   eg, 2 = red, 3 = blue, 4 = XXX, etc.

   ntCalls <- callNTout(dat)
   orig.outGenes.gaps <- outGenes.gaps
   orig.outGenes <- outGenes

   outGenes.gaps <- c(outGenes.gaps, length(outGenes))
   outGenes <- c(outGenes, colnames(ntCalls$ntMat))

   plotHmapPmat(dat, figName = figName1, mode = "gc",
      plotWidth         = 5.5,
      plotHeight        = 6.5,
      plotTitle = "", legend=FALSE,
      geneLabel         = TRUE,
      cluster_cols      = FALSE,
      cluster_rows      = FALSE,
      geneList          = outGenes,
      gaps_genes        = outGenes.gaps,
      colType           = "gene",
      cellGroups        = dat$specs$cellGroups,
      suppMat           = t(ntCalls$ntMat),
      suppColors        = list( "unk"  = "#CCCCCC",
                                "Hist" = "#F26522",
                                "GABA" = "#21409A",
                                "ACh"  = "#7F3F98",
                                "Glu"  = "#00A14B",
                                "DA"   = "#915B41")
   )
   
   plotHmapPmat(dat,
      figName = paste0(figName1,"_relTPM"),
      mode = "gc",
      exprMode = "raw", rawSamples = dat$specs$sampleLists$all,
      rawNL = "mean",
      plotWidth = 5.5,
      plotHeight = 10,
      plotTitle = "", legend=FALSE,
      geneLabel = TRUE,
      cluster_cols= FALSE,
      cluster_rows= FALSE,
      geneList = orig.outGenes,
      gaps_genes = orig.outGenes.gaps,
      colType = "gene",
      cellGroups = dat$specs$cellGroups.all
   )

   plotHmapPmat(dat, figName = figName2, mode = "gc",
      geneLabel = TRUE,
      cluster_cols= FALSE,
      cluster_rows= FALSE,
      fontsize_col=5,
      fontsize_row=5,
      plotTitle = "", legend=FALSE,
      plotWidth = 6,
      plotHeight=7,
      geneList = inGenes,
      colType = "gene",
      gaps_genes = gaps_genes,
      cellGroups = dat$specs$cellGroups
   )
   
   plotHmapPmat(dat,
      figName = paste0(figName2,"_relTPM"),
      mode = "gc",
      exprMode = "raw", rawSamples = dat$specs$sampleLists$all,
      rawNL = "mean",
      geneLabel = TRUE,
      cluster_cols= FALSE,
      cluster_rows= FALSE,
      fontsize_col=5,
      fontsize_row=5,
      plotTitle = "", legend=FALSE,
      plotWidth = 6,
      plotHeight=10,
      geneList = inGenes,
      colType = "gene",
      gaps_genes = gaps_genes,
      cellGroups = dat$specs$cellGroups.all
   )
   
}


plotHmapPmat.cam <- function(dat, figName,
                             mode = "sparse" #"all"
                             ) {

# Purpose:  make heatmap of cell adhesion molecules

   geneIDs <- dat$specs$goInfo$gene_id[ dat$specs$goInfo$go_term == "cell adhesion"]
   geneList <- dat$specs$transcriptInfo$gene_name[dat$specs$transcriptInfo$gene_id %in% geneIDs]
   
   geneList <- intersect(geneList, rownames(dat$pFit$logp_on_gc))
   
   # only select those with sparse patterns = < 50% of cells
   
   tMat <- dat$pFit$logp_on_gc[geneList,
            setdiff(colnames(dat$pFit$logp_on_gc), c("min","max"))]

   if (mode=="sparse") {
      geneSparsity <- apply(tMat,1, function(x) {sum(exp(x) < 0.2)})
      names(geneSparsity) <- rownames(tMat)
      geneList <- names(geneSparsity[geneSparsity > 0.3 * ncol(tMat)])

      tMat <- tMat[geneList,]

   }
   
   plotHmapPmat(dat, figName = figName, mode = "gc",
      fontsize_col=3,
      plotWidth=5,
      plotHeight=6.85,
      plotTitle="",legend=FALSE,
       geneLabel = TRUE,
       cluster_cols = TRUE,
       cluster_rows = FALSE,
       geneList = geneList,
       colType = "gene",
       cellGroups = dat$specs$cellGroups
   )
   
}


plotDetectionVsYield <- function(data,
                                 figName,
                                 cdna.conc_to_ug = 0.03,
                                 yieldType = "yield.cdna.ng_uL" # yield.cdna.ng_uL or yield.nuclei.thousands
                                ) {

# Purpose: Plot number of genes detected vs cDNA yield

   curSamples <- sort(intersect(data$specs$sampleLists$all,
      data$specs$sampleInfo$sample_name[!is.na(data$specs$sampleInfo[,yieldType]) &
                                        data$specs$sampleInfo[,yieldType] != "unk"]))
   outPre <- ""
   if (!missing(figName)) { outPre <- paste0(figName,".") }

   tpmCols <- paste0("tpm.",curSamples)
   plotMat <- data$expr$geneExpr[,tpmCols]

   numDet <- vector()
   for (tpmCol in tpmCols) {
       numDet[tpmCol] <- sum(plotMat[,tpmCol] > 0) }

   curYield <- as.numeric(data$spec$sampleInfo[
                           match(curSamples, data$spec$sampleInfo$sample_name),
                           yieldType])

   if (yieldType == "yield.cdna.ng_uL") {
      curYield <- curYield * cdna.conc_to_ug
   }

   pdf(paste0(data$specs$outDir,"/",outPre,"detection_vs_",yieldType,".pdf"),
       height=2,width=2.5)
   par(mar=c(3,4.0,0.5,0.5),ps=10)

   if (yieldType == "yield.nuclei.thousands") {
      xLab <- "nuclear yield (thousands)"
   } else {
      xLab <- "cDNA yield (ug)"
   }

   plot(curYield, numDet,
           pch=20,
           cex=0.75,
           log="x",
           xlab="",
           ylab="",
           las = 1,
           cex.axis=0.8,
           xaxt="n",
           main="")

   axis(1,
           cex.axis=0.8,
           mgp=c(3,0.5,0))
   mtext(xLab, 1, line=1.5)
   mtext("genes detected", 2, line=3.2)
   dev.off()

}


plotTxMolVsYield <- function(data,
                             figName,
                             yieldType = "yield.nuclei.thousands" # yield.cdna.ng_uL or yield.nuclei.thousands
                            ) {

# Purpose: plot ERCC-estimated number of transcript molecules vs yield

   curSamples <- sort(intersect(data$specs$sampleLists$all,
      data$specs$sampleInfo$sample_name[!is.na(data$specs$sampleInfo[,yieldType]) &
                                        data$specs$sampleInfo[,yieldType] != "unk"]))
   outPre <- ""
   if (!missing(figName)) { outPre <- paste0(figName,".") }

   tpmCols <- paste0("tpm.",curSamples)
   plotMat <- data$expr$geneExpr[,tpmCols]

   numMols <- plot.ERCC.calibration.curves(data, returnMols = TRUE)
   numMols <- numMols[numMols$sample_name %in% curSamples,]
   numMols <- numMols[match(curSamples, numMols$sample_name),]

   if (yieldType == "yield.cdna.ng_uL") {
      xLab <- "cDNA yield (ng/uL)"
   } else {
      xLab <- "nuclear yield (thousands)"
   }

   pdf(paste0(data$specs$outDir,"/",outPre,"transcripts_vs_",yieldType,".pdf"),
       height=2,width=2.5)
   par(mar=c(3,4.0,0.5,0.5),ps=10)
   plot(as.numeric(numMols[,yieldType]),
        numMols$numEndogMols,
        pch=20,
        cex=0.75,
        cex.axis=0.8,
        log="xy",
        xlab="",
        ylab="",
        las = 1,
        xaxt="n",
        main="")

   axis(1, cex.axis=0.8, mgp=c(3,0.5,0))
   mtext(xLab, 1, line=1.5)
   mtext("transcript molecules", 2, line=3.2)
   dev.off()

}



plotRepCorVsYield <- function(data,
                              figName,
                              cdna.conc_to_ug = 0.03,
                              yieldType = "yield.cdna.ng_uL" # yield.cdna.ng_uL or yield.nuclei.thousands
                            ) {

# Purpose: Plot replicate correlation vs yield

   curSamples <- sort(intersect(data$specs$sampleLists$all,
      data$specs$sampleInfo$sample_name[!is.na(data$specs$sampleInfo[,yieldType]) &
                                        data$specs$sampleInfo[,yieldType] != "unk"]))
   outPre <- ""
   if (!missing(figName)) { outPre <- paste0(figName,".") }

   tx1 <- compareReplicates(data, returnDatOnly=TRUE, mode="bioReps",
                            plotRepScatter = FALSE,
                            selectFrom = curSamples)

   tx2 <- as.data.frame(rbind(as.matrix(tx1[,c("rep1","cor")]),
                              as.matrix(tx1[,c("rep2","cor")])),
                        stringsAsFactors=FALSE) %>%
          group_by(rep1) %>%
          summarise_all(funs(max))
   colnames(tx2) <- c("sample_name", "maxCor")
   tx2 <- as.data.frame(tx2)

   tx1 <- merge(data$specs$sampleInfo, tx2)

   if (yieldType == "yield.cdna.ng_uL") {
      tx1[,yieldType] <- tx1[,yieldType] * cdna.conc_to_ug
      xLab <- "cDNA yield (ug)"
#      xLab <- "cDNA yield (ng/uL)"
   } else {
      xLab <- "nuclear yield (thousands)"
   }

   pdf(paste0(data$specs$outDir,"/",outPre,"replicateCor_vs_",yieldType,".pdf"),
       height = 2, width  = 2.5)
   par(mar=c(3,4.0,0.5,0.5),ps=10,xpd=TRUE)

   plot(as.numeric(tx1[,yieldType]),
        tx1$maxCor,
        pch=20,
        cex=0.75,
        las = 1,
        log="x",
        xlab="",
        ylab="",
        yaxt="n",
        xaxt="n",
        main="")

   mtext(side = 1, text = xLab, line = 1.5)
   mtext(side = 2, text = "Replicate correlation\n(max pearson R)", line = 2)

   axis(1, las=1, xpd=FALSE,cex.axis=0.8, mgp=c(3,0.5,0))
   axis(2, las=1, xpd=FALSE,cex.axis=0.8)

   dev.off()

}



evaluatePcalls <- function(data,
                           figName = "msFigS3",
                           summaryOnly = FALSE,
                           labelGenes = TRUE
                           ) {
# Evaluate p(on) calls against a benchmark of known + and - genes

   library(png)

   pMat <- data$pFit$logp_on_gc

   numMatches <- 0
   numMisMatches <- 0
   bm.pons <- list(on = c(), off = c())
   bm.ponsGenes <- list(on = c(), off = c())
   bm.ponsCells <- list(on = c(), off = c())
   usedBM.cells <- c()
   usedBM.genes <- c()

   bmSet  <- data.frame("cell" = character(),
                        "gene" = character(),
                        "state" = character(),
                        "ref" = character(),
                        stringsAsFactors = FALSE)
   mismatches <- list()
   for (i in 1:nrow(data$specs$bmExpr)) {
      curGene <- data$specs$bmExpr$geneName[i]
      curCell <- data$specs$bmExpr$celltype[i]
      curRef  <- data$specs$bmExpr$references[i]
      expState <- data$specs$bmExpr$expression[i]

      if ( (! grepl("immuno", data$specs$bmExpr$evidence[i])) &
           (! grepl("good", data$specs$bmExpr$evidence[i])) &
           (! grepl("protein", data$specs$bmExpr$evidence[i])) &
           (! grepl("western", data$specs$bmExpr$evidence[i]))) next;

      if (curCell %in% c("opticlobe", "medulla", "lamina")) {next;}

      if (! curGene  %in% rownames(pMat)) {
         print(paste0("dont have benchmark gene ",curGene," in pMat"))
         next;
      }

      curCells <- curCell

      if (curCell == "R" | curCell == "R1-8") {
         curCells <- c("R1-6", "R7", "R8_Rh5", "R8_Rh6")
      } else if (curCell == "R1-6") {
         curCells <- c("R1-6")
      } else if (curCell == "R7") {
         curCells <- c("R7")
      } else if (curCell == "R8") {
         curCells <- c("R8_Rh5", "R8_Rh6")
      } else if (curCell == "R7-8") {
         curCells <- c("R7", "R8_Rh5", "R8_Rh6")
      } else if (curCell == "Rh5") {
         curCells <- c("R8_Rh5")
      } else if (curCell == "Rh6") {
         curCells <- c("R8_Rh6")
      }

      for (j in 1:length(curCells)) {
         if (! curCells[j] %in% colnames(pMat)) {
            next;
         }

         bmSet[nrow(bmSet) + 1,] <- c(curCells[j], curGene, expState, curRef)

         usedBM.cells <- unique(c(usedBM.cells, curCells[j]))
         usedBM.genes <- unique(c(usedBM.genes, curGene))

         obsState <- 'ambiguous (0.2 < p < 0.8)'
         if (exp(pMat[curGene, curCells[j]]) < 0.2)  {
            obsState <- '-'
         } else if (exp(pMat[curGene, curCells[j]]) > 0.8)  {
            obsState <- '+'
         }

         if (obsState == expState) {
            numMatches <- numMatches + 1
         } else {
            numMisMatches <- numMisMatches + 1
            print(paste0("*******MISMATCH: ",curGene," in ",curCells[j],
                         " obs=",obsState,
                         " vs exp=",expState,
                         " in ",curCell,
                         " REF: ",curRef))
            mismatches[[curGene]] <- c(mismatches[[curGene]], curCells[j])
         }

         if (expState == '+') {
            bm.pons$on <- c(bm.pons$on, exp(pMat[curGene, curCells[j]]))
            bm.ponsGenes$on <- c(bm.ponsGenes$on, curGene)
            bm.ponsCells$on <- c(bm.ponsCells$on, curCells[j])
         } else {
            bm.pons$off <- c(bm.pons$off, exp(pMat[curGene, curCells[j]]))
            bm.ponsCells$off <- c(bm.ponsCells$off, curCells[j])
         }
      }
   }

# print bmSet to file


   print(paste0("numMatches = ",numMatches,
                ". numMisMatches = ",numMisMatches,
                " agreement=",
                sprintf("%.2f",100 * numMatches / (numMatches + numMisMatches))))

   print(paste0("usedBM.cells: ",length(usedBM.cells)))
   print(paste0("usedBM.genes: ",length(usedBM.genes)))

   if (summaryOnly) {
      print("BENCHMARK OVERVIEW: ")
      print(paste0(" total # bm points n=", nrow(bmSet)))
      print(paste0(" -> positive n=", sum(bmSet$state == "+")))
      print(paste0(" -> negative n=", sum(bmSet$state == "-")))
      print(paste0(" -> unique cells: ",length(unique(bmSet$cell))))
      print(paste0(" -> unique genes: ",length(unique(bmSet$gene))))

      return(1);
   }

   write.table(bmSet,
               paste0(dat$specs$outDir,"/",figName,"_benchmark_entries.tsv"),
               quote=FALSE,
               row.names=TRUE,
               col.names=NA,
               sep="\t")

   don <- density(bm.pons$on)
   doff <- density(bm.pons$off)

   tmppngfn <- tempfile()
   png(file = tmppngfn, height=3.1, width=3.1, units="in", res=300,
       family = "ArialMT")
   par(mar = c(0,0,0,0))

   plot(ecdf(bm.pons$on),
        main="",
        lwd=2,
        las=1,
        cex=2,
        xlab="",ylab="",xaxt="n",yaxt="n",
        verticals=T, col.01line = NULL, pch=20)

   plot(ecdf(as.vector(exp(pMat))),
        main="",
        lwd=2,
        cex=2,
        col="darkgray",
        las=1,
        verticals=T, col.01line = NULL, pch=20, add=TRUE)

   dev.off()
   pngbg <- readPNG(tmppngfn)
   pngbg <- as.raster(pngbg)

   pdf(paste0(data$specs$outDir, "/",figName,"_benchmark_curves.pdf"),
       width=2.5, height=2)
   par(mar=c(3,3.5,0.5,1),ps=10)

   plot(ecdf(bm.pons$on),
        main="",
        lwd=2,
        las=1,
        xlab="",ylab="",xaxt="n",yaxt="n", do.points=FALSE,
        verticals=T, col.01line = NULL, pch=20)

   lim <- par()

   rasterImage(pngbg, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])

   mtext("P(expression)", 1, line=1.5)
   mtext("Cumulative fraction of\npositive benchmark genes", 2, line=1.8)
   axis(1, cex.axis=0.8, mgp=c(3,0.5,0), las=1)
   axis(2, cex.axis=0.8, las=1)

   mismatchX <- bm.pons$on[bm.pons$on < 0.8]
   mismatchLab <- bm.ponsGenes$on[bm.pons$on < 0.8]
   text(mismatchX, 0.2, mismatchLab, srt=60,cex=0.5,adj=c(0,0))

   lines(c(0,1,1), c(0,0,1),
         col="darkOrange3",
         lwd=2)

   legend("topleft",
          legend=c("Observed", "Perfect", "Random"),
          text.col=c("black","darkOrange3","darkgray"), bty="n",
          y.intersp = 0.8, x.intersp=0)

   dev.off()

   return(list(usedBM.cells = usedBM.cells, usedBM.genes = usedBM.genes,
               bm.pons = bm.pons,
               bmSet = bmSet,
               mismatches = mismatches))
}


loadBenchmark <- function(specs) {
# Purpose: Load benchmark expression patterns.

   print("Load sample information")
   bmExpr <- read.table( specs$benchmarkLitFn,
                         header       = TRUE,
                         sep          = "\t",
                         as.is        = TRUE,
                         comment.char = '#')


}


buildToggleTree <- function( data,
                             figName,
                             treeType   = "fastme_bal", #hca, pvclust
                             nCores     = 1, # for pvclust
                             alignLabels = TRUE,
                             makePlots  = TRUE,
                             plotHeight = 6,
                             plotWidth = 4,
                             hiSNR.only = FALSE, 
                             geneList   = c("gl"),
                             treeGeneType = "TF", #TF or GO term or "connectome" #Genes used to build tree
                             sampleType = "cell", #or driver or sample
                             geneType,            #
                             queryCell,
                             nodeLabels = FALSE) {

# Purpose: build tree based on number of TF state flips

   library(ape)
   library(phytools)
   library(pheatmap)
   library(dendextend) #for nice plotting; esp pvclust trees

   cellList <- setdiff(colnames(data$pFit$logp_on_gc), c("min","max"))
   cellList <- cellList[!cellList %in% c("VGlut", "Gad1", "ChAT")]

   eMat <- exp(data$pFit$logp_on_gc[,cellList])

   if (hiSNR.only) {
      muons <- vector()
      muoffs <- vector()
      for (gene in names(data$pFit$geneEfit)) {
         if (!"bimPars" %in% names(data$pFit$geneEfit[[gene]])){next;}
         muons[gene] <- data$pFit$geneEfit[[gene]]$bimPars$mu_on
         muoffs[gene] <- data$pFit$geneEfit[[gene]]$bimPars$mu_off
      }
      dmu <- muons - muoffs

      genes.hiSNR <- intersect(names(muons)[muons >= 3], names(dmu)[dmu >= 1])
      print(paste0("ORIG: ",nrow(eMat)," genes"))
      eMat <- eMat[genes.hiSNR,]
      print(paste0("-> hiSNR filter: ",nrow(eMat)," genes"))
   }

   fullMat <- eMat
   cellList <- colnames(eMat)
   nCells <- length(cellList)

# update cellnaming to avoid ':' characters
   newNames <- gsub("\\:",".",cellList)
   cellList <- newNames
   colnames(eMat) <- cellList
   colnames(fullMat) <- cellList


   tfGenes <- loadTFs.flybase(data)
   tfGenes <- sort(intersect(tfGenes, rownames(eMat)))

   tfMat   <- fullMat[tfGenes,]


   fullMat[fullMat > 0.8] <- 1
   fullMat[fullMat < 0.2] <- 0
   fullMat[fullMat > 0.2 & fullMat < 0.8] <- 0.5


# For plotting number of detected genes vs tree
   onMat <- eMat
   onMat[onMat > 0.8] <- 1
   onMat[onMat < 1] <- 0
   nOnGenes <- colSums(onMat)
   names(nOnGenes) <- colnames(onMat)


   if (missing(geneType)) {geneType <- treeGeneType;}

   if (treeGeneType == "all") {
      tfMat   <- fullMat
   } else if (treeGeneType == "TF") {
      tfMat   <- fullMat[tfGenes,]
   } else if (treeGeneType == "ozkan") {
      ozkanList <- loadOzkanPPI(dat$specs)
      ozkanList <- intersect(ozkanList$proteinInfo$Symbol,
                             rownames(dat$pFit$logp_on_gc))
      tfMat   <- fullMat[ozkanList,]
   } else if (treeGeneType == "connectome") {
      tfMat   <- t(cMat.binary)
   } else if (treeGeneType != "all") { #assume GO term
      curGenes <- data$specs$goInfo$gene_id[ data$specs$goInfo$go_term == treeGeneType]
      curGenes <- unique(data$specs$transcriptInfo$gene_name[data$specs$transcriptInfo$gene_id %in% curGenes])
      curGenes <- sort(intersect(curGenes, rownames(eMat)))
      tfMat <- fullMat[curGenes,]
      treeGeneType  <- gsub(" ", "_", treeGeneType)
      treeGeneType  <- gsub("-", "_", treeGeneType)
      print("SELECTED ONLY GOTERM GENES")
      print(rownames(tfMat))
   }
   nGenes <- nrow(tfMat)


   if (treeType == "dna.fastme") {

      dnaMat <- tfMat
      dnaMat[dnaMat == 0.5] <- "N"
      dnaMat[dnaMat == 0] <- "A"
      dnaMat[dnaMat == 1] <- "C"
      dnaMat <- t(dnaMat)
      dnaX <- as.DNAbin(dnaMat)
      rownames(dnaX) <- colnames(tfMat)

      dnaDmat <- dist.dna(dnaX,model="raw")
      nBootstrap <- 1000

      treeFun <- function(x) fastme.bal(dist.dna(x,model="raw"))
      tr <- treeFun(dnaX)

      outGroupCells <- c( "Muscles_App" )
      print(paste0("ROOTING TREE ******** using ", paste0(outGroupCells,collapse=", ")))
      tr <- root(tr, outGroupCells, resolve.root=TRUE)

      bp <- boot.phylo(tr,dnaX,treeFun,quiet=FALSE,
                       B=nBootstrap, mc.cores=nCores,trees=TRUE)
      edgeSupp <- prop.clades(tr,bp$trees)
      edgeSupp <- round( edgeSupp / nBootstrap, digits=2)
      print(paste0("edgeSupp has ",length(edgeSupp)," items"))

      if (makePlots) {

      pdf(paste0(dat$specs$outDir,"/",figName,
                  "_tree_",treeType,"_",treeGeneType,".pdf"),
          height=plotHeight,
          width=plotWidth)


      if (!alignLabels) {
         plotTree(tr, fsize=0.4, ftype="reg",offset=0.3)
         nodelabels(pie=cbind(100 * edgeSupp,100-(100 * edgeSupp)),
                    piecol=c("red","white"),cex=0.3)
      } else {

         par(fg="transparent")
         plotTree(tr, fsize=0.4, ftype="reg",offset=0.3)
         nodelabels(pie=cbind(rep(100,length(edgeSupp)),
                              rep(0,length(edgeSupp))),
                    piecol=c("black","white"),cex=0.5, pch=1)
         nodelabels(pie=cbind(100 * edgeSupp,100-(100 * edgeSupp)),
                    col="black", piecol=c("red","white"),cex=0.3)
         obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
         par(fg="black")
         text(rep(max(obj$xx[1:Ntip(tr)]),Ntip(tr)),obj$yy[1:Ntip(tr)],
               labels=tr$tip.label,pos=4,cex=0.4)
         for(i in 1:Ntip(tr)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(tr)])),
                               rep(obj$yy[i],2),lty="dotted")
      }
      dev.off()

      }

      return(list(tr = tr, bp=bp,
                  nBootstrap = nBootstrap,
                  edgeSupp = edgeSupp))

   } else if (treeType == "pvclust") {
      tfMat[tfMat == 0.5] <- NA
      cl <- parallel::makePSOCKcluster(nCores)
      tx <- pvclust(data=tfMat, method.dist="manhattan",nboot=100, method.hclust="complete", parallel=cl)
      parallel::stopCluster(cl)


      dend <- as.dendrogram(tx)
      dend <- dend %>% set("labels_cex", 0.4)

      if (makePlots) {
      pdf(paste0(dat$specs$outDir,"/",figName,
                  "_tree_",treeType,"_",treeGeneType,".pdf"))

      dend %>% pvclust_show_signif_gradient(tx,signif_type="au") %>%
         plot(main = "")
      dev.off()
      }

      return(tx)
   }



   dMat <- matrix(nrow=nCells, ncol=nCells, data = NA,
                  dimnames = list(cellList, cellList))

   for (i in 1:nCells) {
      for (j in 1:nCells) {
         dMat[cellList[i],cellList[j]] <-
            sum(abs(tfMat[,i] - tfMat[,j])) / nGenes
      }
   }

   dMat2 <- as.dist(dMat,diag=FALSE,upper=FALSE)


   if (treeType == "fastme_bal") {
      tree <- fastme.bal(dMat2)
   } else if (treeType == "fastme_ols") {
      tree <- fastme.ols(dMat2)
   } else if (treeType == "hca") {
      hca <- hclust(dMat2, method="complete")
      tree <- as.phylo(hca)
   }
   print("-> DONE")

   if (makePlots) {
   pdf(paste0(dat$specs$outDir,"/",figName,
              "_TFflip_tree_",treeType,"_",treeGeneType,".pdf"))
   plotTree(tree, fsize=0.4, ftype="i", offset=0.3)
   dev.off()
   }

# types, eg ion channels, the goruping wasn't monophyletic

   outGroupCells <- c( "Muscles_App", "Muscles_Head")

   print(paste0("ROOTING TREE ******** using ", paste0(outGroupCells,collapse=", ")))
   if (treeGeneType != "connectome") {
      tree <- root(tree, outGroupCells, resolve.root=TRUE)
      print(is.rooted(tree))
      tree$edge.length[tree$edge.length==0]<-max(nodeHeights(tree))*1e-6
   } else {
      tree <- root(tree, c( "C2"), resolve.root=TRUE)
      print(is.rooted(tree))
      tree$edge.length[tree$edge.length==0]<-max(nodeHeights(tree))*1e-6
   }

   if (makePlots) {
   pdf(paste0(dat$specs$outDir,"/",figName,
              "_TFflip_tree_",treeType,"_nOnGenes.pdf"))
   contMap(tree, nOnGenes, type="phylogram",fsize=0.4,ftype="i",ooffset=0.3)
   dev.off()
   }


# reconstruct ancestral state of geneList genes, and plot as color on tree
   for (gene in geneList) {
      curVec <- as.character(fullMat[gene,])
      names(curVec) <- colnames(fullMat)
      x <- curVec

      x <- x[tree$tip.label]

      print(paste0("Ancestral reconsturction of ", gene))
      fitER <- ace(x ,phy=tree,model="ER",type="discrete")
      print("-> DONE")


      cols<-setNames(palette()[1:length(unique(x))],sort(unique(x)))

      if (makePlots) {
      pdf(paste0(dat$specs$outDir,"/",figName,
                 "_TFflip_tree_",treeType,"_",treeGeneType,"_asr_",gene,".pdf"))

      plotTree(tree,fsize=0.4,ftype="i", offset=0.3, main=gene)

      tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)

      if (nodeLabels) {
      nodelabels(node=1:tree$Nnode+Ntip(tree),
          pie=fitER$lik.anc,piecol=cols,cex=0.3)
      nodelabels(bg="white",cex=0.3,)
      }

      add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                        y=-max(nodeHeights(tree)),fsize=0.5,
                        vertical=FALSE)

      dev.off()
      }
   }

   tfMat <- tfMat[,tree$tip.label]
   fullMat <- fullMat[,tree$tip.label]

   if (geneType == "TF") {
      tfMat   <- fullMat[tfGenes,]
   } else if (geneType == "ozkan") {
      tfMat   <- fullMat[ozkanList,]
   } else if (geneType != "all") { #assume GO term
      curGenes <- data$specs$goInfo$gene_id[ data$specs$goInfo$go_term == geneType]
      curGenes <- unique(data$specs$transcriptInfo$gene_name[data$specs$transcriptInfo$gene_id %in% curGenes])
      curGenes <- sort(intersect(curGenes, rownames(eMat)))
      tfMat <- fullMat[curGenes,]
      geneType  <- gsub(" ", "_", geneType)
      geneType  <- gsub("-", "_", geneType)
      print("SELECTED ONLY GOTERM GENES")
      print(rownames(tfMat))
   } else {
      tfMat <- fullMat
   }


   hits <- c()
# Search for group-v-rest pattern match
   for (node in 1:(nCells+tree$Nnode)) {

      descCells <- tree$tip.label[
         intersect(getDescendants(tree,node), 1:nCells)]

      queryPattern <- vector(mode="numeric",length=nCells)
      names(queryPattern) <- tree$tip.label

      queryPattern[descCells] <- 1
      curPattString <- paste(sort(descCells), collapse=",")

      for (gene in rownames(tfMat)) {
         curScore <- sum(abs(queryPattern - tfMat[gene,]))
         if (curScore > 0.2 * sum(queryPattern)) {next;}
         if (min(tfMat[gene,]) > 0) {next;}
         hits <- unique(c(hits, gene))

         if (!missing(queryCell)) {
         if (any(grep(queryCell, descCells)) & length(descCells) < 10) {
            print(paste0("****************************************",
                         curPattString," vs rest: ", gene))
         }
         }

      }
   }


# Search for L-vs-R descendant pattern match
   nodeDesc <- list()
   for (node in (nCells + 1):(nCells+tree$Nnode)) {
      curDesc <- tree$edge[tree$edge[,1] == node,2]
      desc.R <- tree$tip.label[intersect(getDescendants(tree,curDesc[1]), 1:nCells)]
      desc.L <- tree$tip.label[intersect(getDescendants(tree,curDesc[2]), 1:nCells)]
      allDesc <- c(desc.R, desc.L)
      print(paste0("Looking for genes that mark ", paste0(desc.R,collapse=", "),
                   " vs ", paste0(desc.L,collapse=", ")))
      curPattString <- paste0(paste(sort(desc.R), collapse=","), " vs ",
                              paste(sort(desc.L), collapse=","))


      if (length(desc.R) < 2 & length(desc.L) < 2 ) { next; }

      curTfMat <- tfMat[,allDesc]

      queryPattern <- list()
      queryPattern$r <- vector(mode="numeric",length=length(allDesc))
      queryPattern$l <- vector(mode="numeric",length=length(allDesc))
      names(queryPattern$r) <- allDesc
      names(queryPattern$l) <- allDesc

      queryPattern$r[desc.R] <- 1
      queryPattern$l[desc.L] <- 1

      for (gene in rownames(curTfMat)) {

         curScore1 <- sum(abs(queryPattern$r - curTfMat[gene,]))
         curScore2 <- sum(abs(queryPattern$l - curTfMat[gene,]))

         if (curScore1 > 0.1 * sum(queryPattern$r) &
             curScore2 > 0.1 * sum(queryPattern$l)) {next;}

         if (min(tfMat[gene,]) > 0) {next;}

         if (!missing(queryCell)) {
         if (any(grep(queryCell, allDesc)) & length(allDesc) < 30) {
            print(paste0("-> ",gene))
         }
         }

         hits <- unique(c(hits, gene))
      }
   }
   print(paste0("hits = ",hits))

   pdf(paste0(dat$specs$outDir,"/",figName,
              "_TFflip_tree_",treeType,"_",treeGeneType,"_hmap_allpatt.pdf"),
       onefile=FALSE, height=5,width=8)

   pheatmap(fullMat[hits,],
            breaks = seq(0,1,length.out=101),
            cluster_cols=as.hclust(chronos(tree)),
            fontsize=8,
            treeheight_row = 0,
            treeheight_col = 0,
            fontsize_row=3,
            fontsize_col=4,
            border_color = NA,
            col=colorRampPalette(c("steelBlue3","white","darkOrange3"))(100))
   dev.off()

   return(list(
      dMat = dMat,
      dMat2 = dMat2,
      tree = tree
   ))

}


writeDataTables <- function(dat, tabName = "dataTable") {

# Purpose: Write Data tables for GEO submission

   timeStamp <- format(Sys.time(), "%Y%m%d")

   fh.readme<-file(paste0(dat$specs$outDir,"/README_",timeStamp,".md"), "w")
   cat("README: supplemental data tables for Davis et al., 2018\n",file=fh.readme)
   cat("-------------------------------------------------------\n",file=fh.readme)
   cat("\n",file=fh.readme)
   cat("\n",file=fh.readme)
   cat("These data tables report our expression measurements at several levels of\n", file=fh.readme)
   cat("processing, from sample-level transcript abundance to cell type-level binary\n", file=fh.readme)
   cat("expression calls.\n", file=fh.readme)
   cat("",file=fh.readme)
   cat("",file=fh.readme)
   cat("1. transcript x sample. transcript-level abundance for all genes in all samples\n",file=fh.readme)
   cat("   (includes coding, noncoding, INTACT reporter, and ERCC spike-ins).\n",file=fh.readme)

   print("- Data table 1")
   {
      outFn <- paste0(tabName,"1.transcripts_x_samples_TPM.all_genes.txt.gz")
      outTab <- dat$expr$txExpr.withERCC[,
         c("transcript_id","gene_id","gene_name",
           paste0("tpm.", dat$specs$sampleInfo$sample_name))]
      tpmCols <- colnames(outTab)[grep("^tpm", colnames(outTab))]
      outTab[,tpmCols] <- round(outTab[,tpmCols], digits=2)
      colnames(outTab) <- gsub("^tpm.","",colnames(outTab))

      gz1<-gzfile( paste0(dat$specs$outDir,"/",outFn), "w") ;
      write.table(outTab, file=gz1, col.names=TRUE,
                  row.names=FALSE, sep="\t", quote=FALSE);
      close(gz1) ;
   }

   print("- Data table 2")
   cat("\n",file=fh.readme)
   cat("2. gene x sample. gene-level abundance for all genes in all samples\n",file=fh.readme)
   cat("   (includes coding, noncoding, INTACT reporter, and ERCC spike-ins).\n",file=fh.readme)
   {
      outFn <- paste0(tabName,"2.genes_x_samples_TPM.all_genes.txt.gz")
      outTab <- dat$expr$geneExpr.withERCC[,
         c("gene_id","gene_name",
           paste0("tpm.", dat$specs$sampleInfo$sample_name))]
      tpmCols <- colnames(outTab)[grep("^tpm", colnames(outTab))]
      outTab[,tpmCols] <- round(outTab[,tpmCols], digits=2)
      colnames(outTab) <- gsub("^tpm.","",colnames(outTab))

      gz1<-gzfile( paste0(dat$specs$outDir,"/",outFn), "w") ;
      write.table(outTab, file=gz1, col.names=TRUE,
                  row.names=FALSE, sep="\t", quote=FALSE);
      close(gz1) ;
   }

   print("- Data table 3")
   cat("\n",file=fh.readme)
   cat("3. gene x sample. abundance of protein coding genes in all samples.\n",file=fh.readme)
   {
      outFn <- paste0(tabName,"3.genes_x_samples_TPM.coding_genes.txt.gz")
      outTab <- dat$expr$geneExpr[,
         c("gene_id","gene_name",
           paste0("tpm.", dat$specs$sampleInfo$sample_name))]
      tpmCols <- colnames(outTab)[grep("^tpm", colnames(outTab))]
      outTab[,tpmCols] <- round(outTab[,tpmCols], digits=2)
      colnames(outTab) <- gsub("^tpm.","",colnames(outTab))

      gz1<-gzfile( paste0(dat$specs$outDir,"/",outFn), "w") ;
      write.table(outTab, file=gz1, col.names=TRUE,
                  row.names=FALSE, sep="\t", quote=FALSE);
      close(gz1) ;
   }


# NEW per-cell average (good samples only!)
   print("- Data table 4")
   cat("\n",file=fh.readme)
   cat("4. gene x cell. cell type-mean abundance of protein coding genes in QC-passed samples.\n",file=fh.readme)
   {
      outFn <- paste0(tabName,"4.genes_x_cells_TPM.coding_genes_QCpass.txt.gz")
      t.rawMat <- dat$expr$geneExpr[,paste0("tpm.",
                                 dat$specs$sampleLists$good$allCriteria)]
      rownames(t.rawMat) <- dat$expr$geneExpr$gene_name

      t.curCells <- sort(unique(dat$specs$sampleInfo$celltype[
                        dat$specs$sampleInfo$sample_name %in%
                        dat$specs$sampleLists$good$allCriteria]))

      t.gcMat <- matrix(nrow=nrow(t.rawMat),ncol=length(t.curCells),
                        dimnames = list(rownames(t.rawMat), t.curCells))

      for (t.curCell in t.curCells) {
         t.curSamples <- paste0("tpm.",
                                intersect(dat$specs$sampleLists$good$allCriteria,
                                   dat$specs$sampleInfo$sample_name[
                                      dat$specs$sampleInfo$celltype %in%
                                      c(t.curCell)]))
         t.gcMat[,t.curCell] <- apply(t.rawMat[,t.curSamples],1,mean)
      }
      t.gcMat <- round(t.gcMat , digits=2)
      t.gcMat <- as.data.frame(t.gcMat)
      gz1<-gzfile(paste0(dat$specs$outDir,"/",outFn), "w")
      write.table(t.gcMat, file=gz1, col.names=NA, row.names=TRUE,
               sep="\t", quote=FALSE)
      close(gz1)
   }

   cat("\n", file=fh.readme)
   cat("\n", file=fh.readme)
   cat("\n", file=fh.readme)
   cat("Tables 5 - 7. raw levels and inferred expression probabilities for coding genes\n", file=fh.readme)
   cat("              with at least 10TPM abundance in at least one sample.\n", file=fh.readme)
   cat("              NOTE: only considers samples that pass QC.\n", file=fh.readme)

   addCols <- as.data.frame(matrix(NA,ncol=8,nrow=nrow(dat$pFit$logp_on_gc))) ;
   colnames(addCols) <- c( "diff.elpd", "diff.elpd.se",
                           "bimodal.on_level", "bimodal.off_level",
                           "bimodal.sd", "bimodal.pi", "unimodal.level", "unimodal.sd") ;
   rownames(addCols) <- rownames(dat$pFit$logp_on_gc)
   exprType <- vector()
   for (t.gene in rownames(dat$pFit$logp_on_gc)){
      exprType[t.gene] <- dat$pFit$geneEfit[[t.gene]]$exprType
      addCols[t.gene,] <- c(dat$pFit$geneEfit[[t.gene]]$diff.elpd,
                            dat$pFit$geneEfit[[t.gene]]$diff.elpd.se,
                            dat$pFit$geneEfit[[t.gene]]$bimPars$mu_on,
                            dat$pFit$geneEfit[[t.gene]]$bimPars$mu_off,
                            dat$pFit$geneEfit[[t.gene]]$bimPars$sd_on,
                            dat$pFit$geneEfit[[t.gene]]$bimPars$pi_on,
                            dat$pFit$geneEfit[[t.gene]]$uniPars$mu1,
                            dat$pFit$geneEfit[[t.gene]]$uniPars$sd1)
   }
   addCols <- round(addCols,digits=2)
   addCols$exprType[match(names(exprType), rownames(addCols))] <- exprType

   curTabNum <- 5
   for (mode in c("samples", "drivers", "cells")) {
      print(paste0("- Data table ",curTabNum))
   
      nUnits <- ncol(dat$expr$rawpon[[mode]]$raw)

   
      {
         outFn <- paste0(tabName,curTabNum,"a.",
                     "genes_x_",mode,"_TPM.modeled_genes.txt.gz")
         outTab <- round(dat$expr$rawpon[[mode]]$raw, digits=2)
         colnames(outTab) <- gsub("^tpm.","",colnames(outTab))
         rownames(outTab) <- rownames(dat$expr$rawpon[[mode]]$raw)
                 
         cat("\n",file=fh.readme)
         cat(paste0(curTabNum,"a. genes x ",mode," matrix of transcript abundance (TPM).\n"),file=fh.readme)
         gz1<-gzfile( paste0(dat$specs$outDir,"/",outFn), "w") ;
         write.table(outTab, file=gz1, col.names=NA,
                     row.names=TRUE, sep="\t", quote=FALSE);
         close(gz1) ;
      }
   
      {
         outFn <- paste0(tabName,curTabNum,"b.",
                     "genes_x_",mode,"_p_expression.modeled_genes.txt.gz")
         outTab <- round(exp(dat$expr$rawpon[[mode]]$p),digits=2)
         colnames(outTab) <- gsub("^tpm.","",colnames(outTab))
   
         outTab <- merge(outTab, addCols, by="row.names")
         rownames(outTab) <- outTab$Row.names
         outTab$Row.names <- NULL
   
         cat(paste0(curTabNum,"b. genes x ",mode," matrix of expression probabilities.\n"),file=fh.readme)
      cat(paste0("   - columns 1-",nUnits,". ",mode,"\n"),file=fh.readme)
      cat(paste0("   - columns ",(1 + nUnits),"-",ncol(outTab),". model details\n"),file=fh.readme)
      cat("    - diff.elpd: difference in bimodal - unimodal model fit: positive = bimodal better\n",file=fh.readme)
      cat("    - bimodal.on_level: mean log-scale expression level of inferred on component\n",file=fh.readme)
      cat("    - bimodal.off_level: mean log-scale expression level of inferred off component\n",file=fh.readme)
      cat("    - bimodal.sd: mean log-scale standard deviation of bimodal components\n",file=fh.readme)
      cat("    - bimodal.pi_on: inferred fraction of samples expressing the gene\n",file=fh.readme)
      cat("    - unimodal.level: mean log-scale expression level of inferred unimodal fit\n",file=fh.readme)
      cat("    - unimodal.sd: mean log-scale standard deviation of inferred unimodal fit\n",file=fh.readme)
      cat("    - expr.type: inferred expression type: bimodal or unimodal\n", file=fh.readme)
   
         gz1<-gzfile( paste0(dat$specs$outDir,"/",outFn), "w") ;
         write.table(outTab, file=gz1, col.names=NA,
                     row.names=TRUE, sep="\t", quote=FALSE);
         close(gz1) ;
      }

      curTabNum <- curTabNum + 1
   }


   close(fh.readme)

   return(1)

}

vizPfit <- function(dat,
                    figName.tpmHisto,
                    figName.tpm_pon,
                    geneList    = c("ninaE"),
                    pointCol    = "black",
                    pointPch    = 20,
                    pointCex = 1,
                    labelCol    = "darkOrange3",
                    mode        = "drivers", # "sample" or "cell"
                    plotUni     = FALSE,
                    autoAxes    = FALSE,
                    labelPch    = 20,
                    labelCells = c(),         #= c("Rh1"),
                    labelDirectly = FALSE, #vs TRUE will label beside point
                    labelLegend = "*",         #= "R1-6 photoreceptor"
                    labelPtCex  = 1,
                    singlePlot  = FALSE
                    ) {

# Purpose: visualize inferState model fit for individual genes

# if sample, tpm is actual sample-level; if driver or cell, then average over
# the specified sampleType

# NOTE: only uses the high quality samples that went into building the model

   rawMat <- dat$expr$rawpon[[mode]]$raw[geneList,,drop=FALSE]
   pMat <- dat$expr$rawpon[[mode]]$p[geneList,,drop=FALSE]

   for (gene in geneList) {
      print(paste0("on gene ",gene))
      rawTPM <- unlist(log(1 + rawMat[gene,]))
      curPon <- pMat[gene,]

      curXrange <- range(rawTPM)

      labelTpmCols <- labelCells
      labelTPM <-  unlist(rawTPM[labelTpmCols])
      labelPon <-  unlist(curPon[labelTpmCols])

      raw.hist <- hist(rawTPM, plot=FALSE,breaks=50)

# Plot 1: histogram of raw TPM on log scale; overlaid with bi-modal and unimodal component densities
      pdf(paste0(dat$specs$outDir,"/",figName.tpmHisto,"_",gene,"_",mode,"_TPM_histogram.pdf"),
          height=2,width=2.5)
      if (!singlePlot) {
         par(mar=c(3,3.5,0.5,1),ps=10)
      } else {
         par(mar=c(3,2.25,0.5,2.25),ps=10)
      }

      plot(raw.hist,main="",
           xlab="",
           ylab="",
           xlim=curXrange,
           col="grey60",
           cex.axis=0.8,
           border=NA,
           las = 1,
           yaxt="n",
           xaxt="n")

      if (!autoAxes) {
         tpmLabs <- c(1,2,5,10,20,50,
                      100,200,500,
                      1000, 2000,5000,10000,50000)
         tpmLabs <- sort(tpmLabs)
         tpmLabs <- tpmLabs[tpmLabs > min(rawTPM)]

         if (!singlePlot) {
            axis(1,at=log(tpmLabs),labels=tpmLabs,
                 cex.axis=0.8, mgp=c(3,0.5,0))
            axis(2,cex.axis=0.8,mgp=c(3,0.5,0),las=1)
         } else {

            axis(1,at=log(tpmLabs), labels=tpmLabs,
                 cex.axis=0.8, mgp=c(3,0.5,0))
            axis(2,cex.axis=0.8,col="grey60",
                 col.axis="grey60",col.ticks="grey60",mgp=c(3,0.5,0),las=1)
         }
      } else {
         if (! singlePlot) {
            axis(1,at=axTicks(1),cex.axis=0.8, mgp=c(3,0.5,0))
         } else {
            axis(1,at=axTicks(1),cex.axis=0.8, mgp=c(3,0.2,0))
         }
      }

      mtext(paste0(gene," abundance (TPM+1)"), 1, line=1.5)

      if (!singlePlot){
         mtext(paste0("# ",mode),2, line=2.0)
      } else {
         mtext(paste0("# ",mode),2, line=1,col="grey60")
      }


      bim.pion <- dat$pFit$geneEfit[[gene]]$bimPars$pi_on

      curX <- seq(curXrange[1],curXrange[2],0.05)
      fitY.on <- bim.pion * diff(raw.hist$mids[1:2]) * length(rawTPM) *
                 dnorm(curX, dat$pFit$geneEfit[[gene]]$bimPars$mu_on,
                       dat$pFit$geneEfit[[gene]]$bimPars$sd_on)
      lines(curX,fitY.on,lwd=2,col="darkOrange3",main="",type="l",
           xlab="",ylab="",
           xlim = curXrange)

      fitY.off <- (1 - bim.pion) *  diff(raw.hist$mids[1:2]) * length(rawTPM) *
                  dnorm(curX,
                        dat$pFit$geneEfit[[gene]]$bimPars$mu_off,
                        dat$pFit$geneEfit[[gene]]$bimPars$sd_off)
      lines(curX,fitY.off,lwd=2,col="steelBlue3",main="",type="l",
           xlab="",ylab="",
           xlim = curXrange)

      if (plotUni) {
         fitY.uni <- diff(raw.hist$mids[1:2]) * length(rawTPM) * dnorm(curX,
                        dat$pFit$geneEfit[[gene]]$uniPars$mu1,
                        dat$pFit$geneEfit[[gene]]$uniPars$sd1)
         lines(curX,fitY.uni,lwd=2,col="darkgray",main="",type="l",
           xlab="",ylab="",
           xlim = curXrange)

         legend("topright", legend=c("bimodal on", "bimodal off", "unimodal"),
                         text.col=c("darkOrange3","steelBlue3","darkgray"),bty="n")
      } else {
         legend("topright", legend=c("on", "off"),
                         text.col=c("darkOrange3","steelBlue3"),bty="n")
      }


      if (!singlePlot) {
         dev.off()

# Plot 2: estimated p(on) vs TPM for observed data.

         pdf(paste0(dat$specs$outDir,"/",figName.tpm_pon,"_",gene,"_",mode,"_TPM_vs_pon.pdf"),
          height=2,width=2.5)
         par(mar=c(3,3.5,0.5,1),ps=10)
      } else {
         par(mar=c(3,2.25,0.5,2.25),ps=10,new=TRUE)
      }

      plot(exp(rawTPM), exp(curPon),
           xlim=exp(curXrange),
           log="x",
           las = 1,
           main="",
           xlab="",
           ylab = "",
           col = pointCol,
           xaxt="n",
           yaxt="n",
           cex=pointCex,
           cex.axis=0.8,
           pch=pointPch)


      if (!singlePlot) {
         mtext(paste0(gene," abundance (TPM+1)"), 1, line=1.5)
         axis(2,cex.axis=0.8,las=1)
         if (!autoAxes) {
            axis(1,at=tpmLabs,labels=tpmLabs, cex.axis=0.8, mgp=c(3,0.5,0),las=1)
         } else {
            print(axTicks(1))
            axis(1, at=axTicks(1,log=TRUE),
                     labels=exp(axTicks(1,log=TRUE)),
                     cex.axis=0.8, mgp=c(3,0.5,0),las=1)
         }
   
         mtext("P(expression)",2, line=2.0)
      } else {
         axis(4,cex.axis=0.8,las=1,mgp=c(3,0.5,0))
         mtext("P(expression)",4, line=1.25)
      }

      if (length(labelTpmCols) > 0) {
         otherLabel <- "other samples"
         if (labelLegend == "") {
            otherLabel <- ""
         }

         points(exp(labelTPM), exp(labelPon),
                cex = labelPtCex,
                pch=labelPch, col = labelCol)

         if (!labelDirectly) {
            legend( "top", legend=c(labelLegend,
                                 otherLabel),
                  text.col=c(labelCol,"black"),
                  bty="n")
         } else {

            # splitting for visual clarity -- although still need adjustment
            curXlab <- exp(labelTPM)
            curYlab <- exp(labelPon)

            offLabs <- which(curYlab < 0.5)
            onLabs <- which(curYlab >= 0.5)

            if (length(offLabs) > 0) {
            text(curXlab[offLabs], curYlab[offLabs] + 0.1,
                 labels=labelTpmCols[offLabs],
                 col=labelCol,
                 adj=c(0,0.5), srt=90,
                 cex = 0.8)
            }

            if (length(onLabs) > 0) {
            text(curXlab[onLabs], curYlab[onLabs] - 0.1,
                 labels=labelTpmCols[onLabs],
                 col=labelCol,
                 adj=c(1,0.5), srt=90,
                 cex = 0.8)
            }
         }
      }

      dev.off()
   }

}


findReg.NT <- function(dat, figName) {
# Purpose: predict TFs that regulate NT phenotype

   tfList <- loadTFs.flybase(dat)
   tfList <- intersect(tfList, rownames(dat$pFit$logp_on_gc))

   cur.thresh <- list()
   cur.thresh$discOn <- 0.8
   cur.thresh$discOff <- 0.2
   cur.thresh$marker.minFracOverlap <- 0.2  #TF should express in at least this fraction of NT+ cells

   cellGroups <- c(dat$specs$cellGroups, list(c("ChAT", "Gad1","VGlut")))
   discMat <- exp(dat$pFit$logp_on_gc[,c(unlist(cellGroups))])

   discMat[discMat < cur.thresh$discOff] <- 0
   discMat[discMat >= cur.thresh$discOff & discMat < cur.thresh$discOn] <- 0.5
   discMat[discMat >= cur.thresh$discOn] <- 1

   markers <- c("Hdc", "Gad1", "VAChT", "VGlut", "Vmat", "prt")
   marker.cors <- list()
   hmapGenes <-  c()
   geneGaps <- c()
   for (marker in markers) {

      marker.numOn <- sum(discMat[marker,] == 1)

      hmapGenes <- c(hmapGenes, marker)
      print(marker)
      marker.cors[[marker]] <- apply(discMat,
                         1,
                         function(x,y) {
                           if (sum(x == 1) > 1) {
                              score <- (sum(x == 1 & y == 1) /
                                        sum(x == 1))
                           } else {
                              score <- 0
                           }
                           return(score)
                         },
                         y = discMat[marker,])
      
      marker.cors[[marker]] <- data.frame(
          score = marker.cors[[marker]],
          numOn = apply(discMat,1, function(x) {return(sum(x == 1))})
      )
      marker.cors[[marker]]$numOverlap <- marker.cors[[marker]]$numOn * 
                                          marker.cors[[marker]]$score
      rownames(marker.cors[[marker]]) <- rownames(discMat)
      marker.cors[[marker]] <- marker.cors[[marker]][tfList,]

# cutoff: top-10 of score-ranking, at least 20% of marker+ cells
      candGenes <- marker.cors[[marker]][marker.cors[[marker]]$numOverlap >=
                                         cur.thresh$marker.minFracOverlap * marker.numOn,]

      candGenes <- head(candGenes[order(candGenes$score,candGenes$numOn,decreasing=TRUE),],n=10)
      hmapGenes <- c(hmapGenes, rownames(candGenes))

      if (marker != markers[length(markers)]) {
         geneGaps <- c(geneGaps, length(hmapGenes)) }

      #print(t(t(head(sort(marker.cors[[marker]][names(marker.cors[[marker]]) %in% tfList],decreasing=TRUE),n=30))))
   }

   print(hmapGenes)

# Make heatmap of candidate TFs, cellgrouped as before, in blocks with NT marker as first column

   plotHmapPmat(dat, figName = figName, mode = "gc",
       geneLabel = TRUE,
       cluster_cols= FALSE,
       cluster_rows= FALSE,
       geneList = hmapGenes,
       gaps_genes = geneGaps,
       plotTitle = "", legend=FALSE,
       plotWidth = 8,
       plotHeight = 7,
       colType = "gene",
       cellGroups = cellGroups 
   )

   return(marker.cors)

}



findCellTypeDEG.limmavoom <- function(dat,
                                      sampleType="all", # "good"
                                      mode = "lax" # 'deg'
                                      ) {
# Purpose: use limma/voom to define cell type-specific genes on all samples
# as well as high-quality subset used for inferState modeling

# useful if want to explore data from # sub-optimal libraries -- eg NPF should
#   be wildly high in NPF line regardless of quality.

   print(paste0("Finding markers in ", mode," mode"))

   library(limma)
   library(pheatmap)

   allSamples <- dat$specs$sampleInfo
   allSamples <- allSamples[allSamples$sample_name %in%
                            dat$specs$sampleLists[[sampleType]],]

   # using same geneExpr matrix as inferState does.
   curMat.counts <- dat$expr$geneExpr[, paste0("est_counts.", allSamples$sample_name)]
   rownames(curMat.counts) <- dat$expr$geneExpr$gene_name

# only genes detected in at least half the samples
   detGenes <- rownames(curMat.counts)[which(apply(curMat.counts,1,function(x) { sum(x > 1)}) > ncol(curMat.counts) / 2)]
   curMat.counts <- curMat.counts[detGenes,]
   print(paste0("Evaluating ", length(detGenes)," genes detected in at least half of the samples"))

   tpmMat <- dat$expr$geneExpr[, paste0("tpm.", allSamples$sample_name)]
   tpmMat <- round(tpmMat,digits=3)
   tpmMat$gene_name <- dat$expr$geneExpr$gene_name
   rownames(tpmMat) <- tpmMat$gene_name
   tpmMat$gene_name <- NULL


# cellType markers; later cellGroup markers (photoreceptor, glia, muscle, neuron)
# - photoreceptor vs all or photoreceptor vs other neurons?
# - glia vs neuron or glia vs all?
# - muscle vs neuron, glia or muscle vs all?

# HERENOW:170910_1620  need to add cellGroup.XX columns to the sampleInfo table
# to have the right comparators for photoreceptor-, glia-, neuron-, and muscle markers

   allSamples$cellGroup.glia <- "not.glia"
   allSamples$cellGroup.glia[grepl("glia", allSamples$celltype)] <- "glia"

   allSamples$cellGroup.pr <- "not.pr"
   allSamples$cellGroup.pr[grepl("^R1", allSamples$celltype)] <- "pr"
   allSamples$cellGroup.pr[grepl("^Rh", allSamples$celltype)] <- "pr"
   allSamples$cellGroup.pr[grepl("R7", allSamples$celltype)] <- "pr"

   allSamples$cellGroup.muscle <- "not.muscle"
   allSamples$cellGroup.muscle[grepl("uscle", allSamples$celltype)] <- "muscle"

   cellGroups <- colnames(allSamples)[grep("^cellGroup", colnames(allSamples))]

   laxMarkers <- list()
   diffGenes <- list()
   for (groupType in c("celltype", cellGroups)) {
      diffGenes[[groupType]] <- list()
      print(paste0("Finding ",groupType,"-specific genes"))
      ct <- factor(allSamples[,groupType])
      design <- model.matrix(~0 + ct)
      colnames(design) <- levels(ct)

      print("Voom")
      v <- voom(counts = curMat.counts, design = design, normalize="quantile")

      if (mode == "deg") {
      print("->lmFit")
      fit <- lmFit(v, design)
      }

      allGroups <- unique(allSamples[,groupType])
      print(paste0("Comparing across ",length(allGroups)," ",groupType))

      for (group in sort(allGroups)) {
         print(paste0("Finding ",groupType," ",group,"-specific genes"))
         otherGroups <- setdiff(allGroups,group)

         if (mode == "deg") {
            ct <- relevel(ct, ref=group)
            design <- model.matrix(~ct)
            fit <- lmFit(v, design)
            fit$genes <- rownames(v)
            fit2 <- fit[,c(paste0("ct",otherGroups))]
            fit2 <- eBayes(fit2, trend=TRUE)
   
            results <- decideTests(fit2, lfc=1, adjust.method="BH", p.value=0.05)
   #         print(summary(results))
   
            posSig <- rowSums(results<0)==length(otherGroups)
#         print(paste0("-> positive signature: ", paste0(rownames(fit2[posSig,]),collapse=", ")))
            posGenes <- rownames(fit2[posSig,])
   
         } else {

            t.fgcols <- paste0("est_counts.",allSamples$sample_name[
                                 allSamples[,groupType] == group])
            t.bgcols <- setdiff(colnames(v$E), t.fgcols)
            posGenes <- rownames(v$E[apply(v$E[,t.fgcols],1,mean) > 5 &
                                   apply(v$E[,t.fgcols],1,mean) + 0.5 >
                                   5 * apply(v$E[,t.bgcols],1,mean) + 0.5,])
      
         }
         print(paste0("-> positive signature: ", paste0(posGenes,collapse=", ")))

#         negSig <- rowSums(results>0)==length(otherGroups)
#         print(paste0("-> negative signature: ", paste0(rownames(fit2[negSig,]),collapse=", ")))

         diffGenes[[groupType]][[group]] <- list(
            posSig = posGenes)

         if (1 & length(posGenes) > 0) {
            print("---> drawing heatmap")
            curMat <- tpmMat[posGenes,]

            curMat <- curMat + 1
            curMat <- log2(curMat / apply(curMat,1,mean))

            curMat[curMat < -2] <- -2
            curMat[curMat > 2] <- 2
            curBreaks <- seq(-2, 2, length.out=101)

            curHeight <- 7.5
            curWidth <- 20

            rowFont <- 6

            if (grepl("cellGroup", groupType)) {
               curHeight <- 20 }

            if (nrow(curMat) > 50) { rowFont <- 3.5 }
            if (nrow(curMat) > 80) { rowFont <- 3 }

            pdf(paste0(dat$specs$outDir,"/voom_onemat_sigGenes_heatmap.",mode,".",
                       groupType,"_",group,".pdf"), onefile=FALSE,family="ArialMT",
                height=curHeight, width=curWidth)
            pheatmap(curMat,
                     main="relative expression log2(TPM/mean)",
                     border_color = NA,
                     show_rownames = TRUE,
                     show_colnames = TRUE,
#                     gaps_col = tpmCol.gaps,
#                     annotation_col = tpmCol.annot,
#                     annotation_colors = tpmCol.annotColors,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
#                     annotation_legend = TRUE,
#                     annotation_names_col = TRUE,
                     breaks = curBreaks,
                     color = colorRampPalette(c("steelBlue2","white","darkOrange2"))(100),
                     fontsize      = 8,
                     fontsize_row  = rowFont,
                     font.main = 1
            )
            dev.off()
         }

      }

   }

   return(list(diffGenes = diffGenes,
               laxMarkers = laxMarkers))

}


plotHmapPmat.trueMarkers <- function(dat, figName="msFig4F") {

# Purpose: show heatmap of genes that only express in single cell types
# - exception: T4T5/female/male

   pMat <- dat$pFit$logp_on_gc[,setdiff(colnames(dat$pFit$logp_on_gc),
                                        c("min","max"))]

   markerGenes <- list()
   markerGenes <- rownames(pMat[apply(exp(pMat), 1,
            function(x) {return(sum(x >= 0.9))}) == 1, ])
   markerGenes.T4T5 <- rownames(pMat[
      exp(pMat[,"T4.T5"]) >= 0.9 &
      apply(exp(pMat), 1, function(x) {return(sum(x >= 0.9))}) == 1, ])

   markerGenes.C2C3 <- rownames(pMat[
      exp(pMat[,"C2.C3"]) >= 0.9 &
      ( exp(pMat[,"C2"]) >= 0.9 |
        exp(pMat[,"C3"]) < 0.9) &
      apply(exp(pMat), 1, function(x) {return(sum(x >= 0.9))}) == 2, ])

   markerGenes.L1L2 <- rownames(pMat[
      exp(pMat[,"L1.L2"]) >= 0.9 &
      ( exp(pMat[,"L1"]) >= 0.9 |
        exp(pMat[,"L2"]) < 0.9) &
      apply(exp(pMat), 1, function(x) {return(sum(x >= 0.9))}) == 2, ])

   markerGenes <- sort(unique(c(markerGenes,
                                markerGenes.T4T5,
                                markerGenes.C2C3,
                                markerGenes.L1L2)))

   print("Ignoring CR genes")
   markerGenes <- markerGenes[!grepl("^CR",markerGenes)]
   print(paste0(length(markerGenes)," true marker genes"))

   muons <- vector()
   muoffs <- vector()
   for (gene in markerGenes) {
      if (!"bimPars" %in% names(dat$pFit$geneEfit[[gene]])){next;}
      muons[gene] <- dat$pFit$geneEfit[[gene]]$bimPars$mu_on
      muoffs[gene] <- dat$pFit$geneEfit[[gene]]$bimPars$mu_off
   }
   dmu <- muons - muoffs
   markerGenes <- markerGenes[dmu >= 1]
   print(paste0(length(markerGenes)," with dmu >= 1"))

# STRATEGY 1: pick top-n most clearly on genes
   print(paste0(length(markerGenes)," true marker genes"))
   sortMarkers <- markerGenes
   if (0) {
   sortMarkers <- markerGenes[order(dmu[markerGenes],decreasing=TRUE)]
   sortMarkers <- head(sortMarkers,n=40)
   print(sortMarkers)
   }

# Strategy 2: pick best marker for each cell. order genes by cell group so
# nice stair-step

   eMat <- exp(pMat[markerGenes,unlist(dat$specs$cellGroups)])
   print("Emat columns:")
   print(colnames(eMat))
   sortMarkers <- c()
   for (cell in unlist(dat$specs$cellGroups)) {
      print(paste0("ON ", cell))
      curMarkers <- rownames(eMat)[eMat[,cell] >= 0.9]
      print(paste0("->", curMarkers))

      curMarkers <- head(curMarkers[order(dmu[curMarkers],decreasing=TRUE)],n=2)

      sortMarkers <- c(sortMarkers, curMarkers)
   }
   sortMarkers <- unique(sortMarkers)

   plotHmapPmat(dat, figName = figName, mode = "gc",
      fontsize_col=3,
      plotWidth=6,
      plotHeight=7,
      geneLabel = TRUE,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      geneList = sortMarkers,
      colType = "gene",
      cellGroups = dat$specs$cellGroups
   )

   return(1)

}


plotHmapPmat.markers <- function(dat,
                                 cellGroups.marker,
                                 cellGroups.heatmap,
                                 plotHeight = 6,
                                 plotWidth = 6,
                                 useSuboptimal = FALSE,
                                 figName.qc_pass ="msFig4",
                                 tabName = "msTabS4",
                                 figName.withSuboptimal ="msFig4.lax") {

# Purpose: make heatmap of genes marking cell types or groups

   pMat <- exp(dat$expr$rawpon$cells$p)
   rawMat <- dat$expr$rawpon$cells$raw

   if (missing(cellGroups.heatmap)) {
      cellGroups.heatmap <- dat$specs$cellGroups
   }

   if (missing(cellGroups.marker)) {
   cellGroups.marker <- list(
      pr = c("R1-6", "R7","R8_Rh5","R8_Rh6"),
      pr.outer = c("R1-6"),
      pr.inner = c("R7","R8_Rh5","R8_Rh6"),
      glia = c("Glia_Eg", "Glia_Psg", "Glia_Mg"),
      muscle = c("Muscles_App", "Muscles_Head"),
      pigment = c("Tm20", "Dm1"),
      kenyoncell =  c("KC_ab_c", "KC_ab_c.p.s", "KC_ab_p", "KC_ab_s", "KC_gd")
   )
   }

   maxOtherCells <- list(
      pr = 2,
      pr.outer = 2,
      pr.inner = 2,
      glia = 2,
      muscle = 2,
      pigment = 2,
      kenyoncell = 5,
      centralcomplex = 2,
      peptidergic = 2
   )

   pOnThresh <- 0.9

   markerGenes <- list()
# Set up "normal" marker gene search
   for (cellGroup in names(cellGroups.marker)) {
      print("-----------------")
      print(cellGroup)
      curCells <- cellGroups.marker[[cellGroup]]
      otherCells <- setdiff(colnames(dat$pFit$logp_on_gc),curCells)
      otherCells <- setdiff(otherCells, c("min","max"))

      curMaxOn <- length(curCells) + maxOtherCells[[cellGroup]]

      curSet <- rownames(dat$pFit$logp_on_gc)[
         exp(dat$pFit$logp_on_gc[,curCells[1]]) >= pOnThresh]

      if (length(curCells) > 1) {
      for (i in 2:length(curCells)) {
         curSet <- intersect(curSet,
            rownames(dat$pFit$logp_on_gc)[
                exp(dat$pFit$logp_on_gc[,curCells[i]]) >= pOnThresh])
      }
      }


# minTPM across curCells should be higher than maxTPM across all other cells

      print(paste0("curCells:",paste0(curCells,collapse=", ")))

      tpmSet <- rownames(dat$expr$rawpon$cells$raw)[
                           apply(dat$expr$rawpon$cells$raw[,curCells,drop=FALSE], 1, min) >
                           apply(dat$expr$rawpon$cells$raw[,otherCells,drop=FALSE], 1, max)]


      curSet <- intersect(curSet,
         rownames(dat$pFit$logp_on_gc)[
            apply(exp(dat$pFit$logp_on_gc), 1,
               function(x, pOnThresh) {return(sum(x >= pOnThresh))},
               pOnThresh = pOnThresh) <= curMaxOn])

      markerGenes[[cellGroup]] <- intersect(curSet,tpmSet)

      print(paste0("-> found markerGenes for ", cellGroup," n=",
                   length(markerGenes[[cellGroup]])))
   }

   marker.df <- data.frame(
      cellGroup = rep(names(markerGenes), unlist(lapply(markerGenes,length))),
      markerGene = unlist(markerGenes),
      stringsAsFactors = FALSE
   )
   outFn <- paste0(dat$specs$outDir,"/",tabName,"_markerGenes.txt")
   write.table(marker.df,
               file     = outFn,
               col.names= TRUE,
               row.names= FALSE,
               quote    = FALSE,
               sep      = "\t")
   return(1)

   muons <- vector()
   for (gene in names(dat$pFit$geneEfit)) {
      if (!"bimPars" %in% names(dat$pFit$geneEfit[[gene]])){next;}
      muons[gene] <- dat$pFit$geneEfit[[gene]]$bimPars$mu_on
   }

   geneList <- c()
   geneGaps <- c()
   geneGroups <- names(markerGenes)
   for (i in 1:length(geneGroups)) {
      geneList <- c(geneList, markerGenes[[geneGroups[i]]])
      if (i < length(geneGroups)) {
         geneGaps <- c(geneGaps, length(geneList)) }

      print("*******************************************************")
      print(geneGroups[i])
      print(paste0("-> ",length(markerGenes[[geneGroups[i]]])," genes"))
      curGenes <- markerGenes[[geneGroups[i]]]

      curGenes <- curGenes[!grepl("^mt:",curGenes)]
      print(paste0("-> without mt: genes:", length(curGenes)))

      if (! geneGroups[i] %in% names(cellGroups.marker)) {
         curGenes <- curGenes[order(muons[curGenes])]
      } else {
         if (length(cellGroups.marker[[geneGroups[i]]]) > 1) {
         curGenes <- curGenes[order(apply(rawMat[
            match(curGenes,rownames(rawMat)),
            cellGroups.marker[[geneGroups[i]]]],
            1, mean), decreasing=TRUE)]
         } else {
         curGenes <- curGenes[order(rawMat[
                                       match(curGenes,rownames(rawMat)),
                                       cellGroups.marker[[geneGroups[i]]]],
                                    decreasing=TRUE)]
         }
      }
      print(paste0(curGenes, collapse=", "))
      print("*******************************************************")

# Top-10 in expression level; just use pFit mu.on to order.
   }


   plotHmapPmat(dat, figName = figName.qc_pass, mode = "gc",
      fontsize_col=2,
      plotWidth=plotWidth,
      plotHeight=plotHeight,
      legend=FALSE, plotTitle="",
      geneLabel = FALSE,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      gaps_genes = geneGaps,
      geneList = geneList,
      colType = "gene",
      cellGroups = cellGroups.heatmap
   )

   if (useSuboptimal) {
      laxMarkers <- findCellTypeDEG.limmavoom(dat)
      suboptimalCells <- c("Lat", "Pdf", "NPY", "LargeLNv", "kdm2")
      markerGenes2 <- c()
      rowGaps2 <- c()
      for (i in 1:length(suboptimalCells)) {
          t.cell <- suboptimalCells[i] ;
          print(t.cell)
          print("ORIG MARKERS:" )
   
          print(laxMarkers$diffGenes$celltype[[t.cell]]$posSig)
   
          curMarkers <- setdiff(laxMarkers$diffGenes$celltype[[t.cell]]$posSig,
           union(laxMarkers$diffGenes$celltype$Rh1$posSig,
                 laxMarkers$diffGenes$cellGroup.pr$pr$posSig)) ;
   
          print("AFTER PR SUBTRACTION:" )
          print(curMarkers)
          markerGenes2 <- c(markerGenes2, curMarkers) ;
          if (i < length(suboptimalCells)) {
              rowGaps2 <- c(rowGaps2, length(markerGenes2)) ;
          }
      }
   
      hMat2.suboptimal <- dat$expr$geneExpr[
          match(markerGenes2, dat$expr$geneExpr$gene_name),
          paste0("tpm.", dat$specs$sampleLists$all)]
      hMat2.suboptimal <- (hMat2.suboptimal + 1) / apply((hMat2.suboptimal + 1), 1, max)
   
   
      png(paste0(dat$specs$outDir,"/",figName.withsuboptimal,"_markers_suboptimal.png"),
          width=7, height=3, units="in",res=300)
      pheatmap(hMat2.suboptimal,
            scale = "none",
            border_color = NA,
            cluster_cols = TRUE,
            treeheight_col = 0,
            cluster_rows = FALSE,
            show_rownames = FALSE,
            show_colnames = FALSE,
            main = "TPM / max",
            breaks = seq(0,1,length.out=101),
            col=colorRampPalette(c("steelBlue3","white","darkOrange3"))(100),
            gaps_row = rowGaps2)
      dev.off()
   }

   return(list(markerGenes = markerGenes))

}


makeSamplesTable <- function(dat, tabName = "msTableS1") {
# Purpose: make supplemental table of sample info

   allLines  <- dat$specs$sampleInfo

   allLines$QC <- "suboptimal"
   allLines$QC[allLines$sample_name %in% dat$specs$sampleLists$good$all] <- "pass"

   outCols <- c("sampleID",
                "sample_name",
                "driverType",
                "QC",
                "celltype",
                "driver",
                "intact.reporter",
                "intact.protocol",
                "yield.nuclei.thousands",
                "yield.cdna.ng_uL")

   outTab <- allLines[,outCols]

   colnames(outTab) <- c(
      "sampleID.internal",
      "sampleName",
      "driverType",
      "qualityControl",
      "simpleCellTypes",
      "driverID",
      "intact.reporter",
      "intact.protocol",
      "yield.nuclei.thousands",
      "yield.cdna.ng_uL")


   outFn1 <- paste0(dat$specs$outDir, "/",tabName,"A_all_samples.tsv")
   write.table(outTab,file=outFn1,col.names=TRUE,
               row.names=FALSE,quote=FALSE, sep="\t")

   outTab <- unique(allLines[,c(
      "driver",
      "driverType",
      "Identifier..Stable.split.number.or.similar.",
      "simpleCellTypes",
      "Main.cell.type.s.",
      "Additional.expression",
      "How.identified",
      "SPLIT.AD.or.Genotype",
      "SPLIT.DBD",
      "Reference")])

   outTab <- outTab[order(outTab$driver),]

   colnames(outTab) <- c(
      "driverID",
      "driverType",
      "driverID.external",
      "celltype",
      "celltype.details",
      "celltype.additional.expression",
      "how.identified",
      "split.AD.or.Genotype",
      "split.DBD",
      "reference")


   outFn2 <- paste0(dat$specs$outDir, "/",tabName,"B_all_drivers.tsv")
   write.table(outTab,file=outFn2,col.names=TRUE,
               row.names=FALSE,quote=FALSE, sep="\t")


   geoTab.samples <- data.frame(sampleID = allLines$sampleID,
                                sampleName = allLines$sample_name,
                                simpleCellTypes = as.character(allLines$simpleCellTypes),
                                organism = "D. melanogaster",
                                driver = allLines$driver,
                                driverType = allLines$driverType,
                                driver.external = as.character(allLines$Identifier..Stable.split.number.or.similar),
                                intact.reporter = allLines$intact.reporter,
                                intact.protocol = allLines$intact.protocol,
                                molecule.type = "nuclear RNA",
                                description = paste0("QC: ", allLines$QC),
                                raw.file = paste0(allLines$sampleID,".fastq.gz"),
                                stringsAsFactors = FALSE)

   geoTab.samples <- geoTab.samples[order(geoTab.samples$sampleID),]
   geoTab.samples$molecule.type[geoTab.samples$driverType == "dissected"] <- "total RNA"

   geoTab.samples$simpleCellTypes[grepl("^KC",geoTab.samples$simpleCellTypes)] <-
      paste0("Kenyon Cell ",geoTab.samples$simpleCellTypes[grepl("^KC",geoTab.samples$simpleCellTypes)])

   colnames(geoTab.samples) <- 
      c("Sample Name",
        "title",
        "source name",
        "organism",
        "characteristics: driver",
        "characteristics: driver type",
        "characteristics: external driver ID",
        "characteristics: INTACT reporter",
        "characteristics: INTACT protocol",
        "molecule",
        "description",
        "raw file")

   outFn3 <- paste0(dat$specs$outDir, "/",tabName,"_geo_submission_samples.tsv")
   write.table(geoTab.samples, file=outFn3,col.names=TRUE,
               row.names=FALSE,quote=FALSE, sep="\t")

   return(geoTab.samples)
}


plotScatterTPM <- function(data, figName = "msFigS2",
                           xSamples,
                           ySamples, 
                           celltype1,
                           celltype2,
                           driver1,
                           driver2,
                           labelGenes,
                           labelGenes.right,
                           labelGenes.left,
                           colorGenes,
                           colorGenes.color = "darkOrange3",
                           colorGenes.legend,
                           mode = "coding",
                           xlab, ylab) {
# Purpose: Draw pairwise scatterplot of expression levels (units of TPM)

   library(png)

   if (!missing(celltype1)) {
      xSamples <- data$specs$sampleInfo$sample_name[
         data$specs$sampleInfo$celltype == celltype1]
   } else if (!missing(driver1)) {
      xSamples <- data$specs$sampleInfo$sample_name[
         data$specs$sampleInfo$driver == driver1]
   }

   if (!missing(celltype2)) {
      ySamples <- data$specs$sampleInfo$sample_name[
         data$specs$sampleInfo$celltype == celltype2]
   } else if (!missing(driver1)) {
      ySamples <- data$specs$sampleInfo$sample_name[
         data$specs$sampleInfo$driver == driver2]
   }


   if (mode == "all.genes") { #remove ERCC
      curMat <- dat$expr$geneExpr.withERCC[
         !grepl("^ERCC[0-9]", dat$expr$geneExpr.withERCC$gene_name),
         c("gene_name",paste0("tpm.",c(xSamples,ySamples)))]
      rownames(curMat) <- curMat$gene_name
      curMat$gene_name <- NULL
      curMat <- sweep(curMat,2,colSums(curMat),`/`) * 1E6

      xval <- apply(1 + curMat[,paste0("tpm.",xSamples)],1,mean)
      yval <- apply(1 + curMat[,paste0("tpm.",ySamples)],1,mean)
      geneNames <- rownames(curMat)

   } else if (mode == "coding") {

      xval <- apply(1 + dat$expr$geneExpr[,paste0("tpm.",xSamples)],1,mean)
      yval <- apply(1 + dat$expr$geneExpr[,paste0("tpm.",ySamples)],1,mean)
      geneNames <- dat$expr$geneExpr$gene_name
   }

# GOAL: PDF plot, with PNG datapoints.

   exprRange <- range(xval, yval)

   tmppngfn <- tempfile()
   png(file = tmppngfn, height=3.1, width=3.1, units="in", res=300,
       family = "ArialMT")
   par(mar = c(0,0,0,0),pty="s")
   plot(xval, yval,
        xlab=xlab, ylab=ylab,
        main="", log="xy",
        xlim = exprRange,
        ylim = exprRange,
        pch=20, cex=0.5)

   if (!missing(colorGenes)) {
      colorX <- xval[geneNames %in% colorGenes]
      colorY <- yval[geneNames %in% colorGenes]
      points(colorX, colorY,
             pch = 20, cex=0.5,
             col = colorGenes.color)
   }

   dev.off()

   pngbg <- readPNG(tmppngfn)
   pngbg <- as.raster(pngbg)

   pdf(paste0(dat$specs$outDir,"/",figName,"_scatter.pdf"),
       height = 2, width  = 2.5)

   par(mar=c(3,4.5,0.5,0),ps=10, xpd=TRUE, pty="s")

   plot(xval, yval,
        xlab="", ylab="",
        main="", log="xy",
        xlim = exprRange,
        ylim = exprRange,
        xaxt = "n", yaxt = "n",
        type = "n",
        las = 1,
        pch=20, cex=0.5)
   lim <- par()

   mtext(side = 1, text = xlab, line = 1.5)
   mtext(side = 2, text = ylab, line = 2.5)
   axis(1, las=1, xpd=FALSE,cex.axis=0.8, mgp=c(3,0.5,0))
   axis(2, las=1, xpd=FALSE,cex.axis=0.8)

   corVal <- cor(x=log(xval), y=log(yval))

   rasterImage(pngbg, 10^lim$usr[1], 10^lim$usr[3],
                      10^lim$usr[2], 10^lim$usr[4])
   legendText <- paste0("Pearson r = ", sprintf("%.2f", corVal))
   legendCol <- c("black")
   if (!missing(colorGenes.legend)) {
      legendText <- c(legendText, colorGenes.legend)
      legendCol <- c(legendCol, colorGenes.color)
   }

   legend("topleft", legendText, text.col = legendCol,
          bty="n", cex=0.75, y.intersp=0.75, x.intersp=0)

# Select top-5 up and down
   allFC <- c(log(xval) - log(yval))

   df <- data.frame(xval, yval, allFC)
   colnames(df) <- c("xval","yval","allFC")
   rownames(df) <- geneNames
   df$maxE <- apply(df[,c("xval","yval")], 1, max)
   df$minE <- apply(df[,c("xval","yval")], 1, min)


   df <- df[order(df$allFC,decreasing=TRUE),]

   df <- df[! grepl("^CR", rownames(df)),]
   df <- df[! grepl(":CR", rownames(df)),]

   if (missing(labelGenes) & missing(labelGenes.right) & missing(labelGenes.left)){

      df <- df[df$maxE >= 10 & df$minE > 1,]
      outliers.up <- head(df,n=5)
      outliers.down  <- tail(df,n=5)

   } else if (!missing(labelGenes)) {

      df <- df[rownames(df) %in% labelGenes,]
      outliers.up <- df
      outliers.down <- data.frame()

   } else if (!missing(labelGenes.right)) {

      outliers.up <- df[rownames(df) %in% labelGenes.right,]
      outliers.down <- df[rownames(df) %in% labelGenes.left,]
   }

   print("LABELED POINT POSITIONS!")
   print(outliers.up)
   print(outliers.down)

   lab_x <- c(outliers.down$xval, outliers.up$xval)
   lab_y <- c(outliers.down$yval, outliers.up$yval)

   lab_text <- c(rownames(outliers.down), rownames(outliers.up))
   lab_textadj <- c(rep(1, nrow(outliers.down)), rep(0, nrow(outliers.up)))
   lab_textx <- c(rep(10, nrow(outliers.down)), rep(5000, nrow(outliers.up)))
   laborder <- order(lab_y)

   texty_logo <- ( 0.15 * log(max(exprRange), base = 10))
   texty_loginc <- ((0.75 * log(max(exprRange), base = 10)) / length(lab_y))

   for (j in 1:length(laborder)) {
      k <- laborder[j]
      lab_texty <- (10**( (j - 1) * texty_loginc + texty_logo))
      text(lab_textx[k], lab_texty, lab_text[k], cex = 0.5, adj = lab_textadj[k])
      segments(lab_textx[k], lab_texty, lab_x[k], lab_y[k], lwd = 1, col = "grey")
   }

   dev.off()
   return(data.frame(
   gene = geneNames,
   xval = xval,
   yval = yval))


}

makeTAPINyieldFig <- function(dat,
                              figName = "msFig1E",
                              cdna.conc_to_ug = 0.03) {

# Purpose: plot yield at each TAPIN step

   plotSamples <- list()
   plotSamples$input <- dat$specs$sampleInfo$sample_name[dat$specs$sampleInfo$driver == "Tm1_d1_input"]
   plotSamples$cap1 <- dat$specs$sampleInfo$sample_name[dat$specs$sampleInfo$driver == "Tm1_d1_capture1"]
   plotSamples$cap2 <- dat$specs$sampleInfo$sample_name[dat$specs$sampleInfo$driver == "Tm1_d1"]
   plotSamples$mock1<- dat$specs$sampleInfo$sample_name[dat$specs$sampleInfo$driver == "Tm1_d1_mock1"]
   plotSamples$mock2<- dat$specs$sampleInfo$sample_name[dat$specs$sampleInfo$driver == "Tm1_d1_mock2"]

   allYields <- list()
   avgYields <- list()
   for (type in names(plotSamples)) {
      avgYields[[type]] <- mean(dat$specs$sampleInfo$yield.cdna.ng_uL[
                                   dat$specs$sampleInfo$sample_name %in%
                                   plotSamples[[type]]]) * cdna.conc_to_ug

      allYields[[type]] <- dat$specs$sampleInfo$yield.cdna.ng_uL[
                                   dat$specs$sampleInfo$sample_name %in%
                                   plotSamples[[type]]] * cdna.conc_to_ug
   }

   heights <- c(avgYields$input, avgYields$cap1, avgYields$mock1,
                avgYields$cap2, avgYields$mock2)
   pdf(paste0(dat$specs$outDir,"/",figName,"_TAPIN_yields.pdf"),
       height=2,width=2)
   par(mar=c(3,3.5,0,1),xpd=TRUE)
   mp <- barplot(rev(heights), main="",
           names.arg=rev(c("Input", "Capture 1", "(mock)",
                           "Capture 2", "(mock)")),
           col=rev(c("darkgreen","darkgreen","grey60","darkgreen","grey60")),
           xlab="",
           xlim=c(0, max(unlist(allYields))),
           space = c(0.4, 0, 0.4, 0, 0.4),
           cex.names=0.6,
           cex.axis=0.6,
           cex.lab=0.6,
           axes=FALSE,
           horiz=TRUE,
           border=NA,
           las=1)
   axis(1,cex.axis=0.6, mgp=c(3, .5, 0))

   points(allYields$mock2,rep(mp[1],2), pch=20,cex=0.5)
   points(allYields$cap2,rep(mp[2],2), pch=20,cex=0.5)
   points(allYields$mock1,rep(mp[3],2), pch=20,cex=0.5)
   points(allYields$cap1,rep(mp[4],2), pch=20,cex=0.5)
   points(allYields$input,rep(mp[5],2), pch=20,cex=0.5)

   mtext("cDNA yield (ug)", 1,line=1.5,cex=0.6)
   dev.off()

}

plotReplicateCorrHistos <- function(data, figName="msFig3I") {

# Purpose: make histogram of replicate correlation values

   repcomp.bioreps <- compareReplicates(data, mode="bioReps",
                                        plot=FALSE,returnDatOnly=TRUE)

   repcomp.drivers <- compareReplicates(data, mode="drivers",
                                        plot=FALSE,returnDatOnly=TRUE)

   repcomp.protocols <- compareReplicates(data, mode="protocols",
                                          plot=FALSE,returnDatOnly=TRUE)
   

   pdf(paste0(dat$specs$outDir,"/",figName,"_repCorR_histogram.pdf"),
          height=2,width=2.5)
   par(mar=c(3,3.5,0.5,1),ps=10)

   xRange <- range(c(repcomp.bioreps$cor,
                     repcomp.drivers$cor,
                     repcomp.protocols$cor))

   d1 <- density(repcomp.bioreps$cor)
   d2 <- density(repcomp.drivers$cor)
   d3 <- density(repcomp.protocols$cor)

   yRange <- range(c(d1$y, d2$y,d3$y))

   plot(density(repcomp.bioreps$cor),
        xlab = "", ylab = "", main="",
        yaxt = "n", xaxt="n",
        xlim=xRange, ylim=yRange,
        lwd=2)

   axis(1, las=1, xpd=FALSE,cex.axis=0.8, mgp=c(3,0.5,0))
   axis(2, las=1, xpd=FALSE,cex.axis=0.8)
   mtext(side = 1, text = "Correlation (pearson R, logTPM)", line = 1.5)
   mtext(side = 2, text = "Density", line = 2)

   lines(density(repcomp.drivers$cor),
        lwd=2, col="darkOrange3")

   lines(density(repcomp.protocols$cor),
        lwd=2, col="steelBlue3")

   legend("topleft",
          legend=c("Replicates", "Alt Drivers", "Alt Protocols"),
          text.col=c("black","darkOrange3","steelBlue3"), bty="n",
          y.intersp = 0.8, x.intersp=0)

   dev.off()

}


plotExprBreadthHisto <- function(data,
                                 figName.histo  = "msFig4G",
                                 figName.genegroups = "msFig4C",
                                 figName.cdf    = "msFig4H"
                                 ) {
# Purpose: draw histogram/CDF of expression breadth for all genes.
#          also stratifies by functional gene group.

   allGenes <- rownames(data$pFit$logp_on_gc)
   print(paste0(length(allGenes)," genes"))

   cell2type <- unique(data$specs$sampleInfo[,c("driverType","celltype")])
   fullMat <- exp(data$pFit$logp_on_gc[allGenes,])
   fullMat <- fullMat[,colnames(fullMat) %in%
                      cell2type$celltype[cell2type$driverType == "pureCellType"]]
   print(paste0("NOTE!: only analyzing across ",ncol(fullMat)," pure cell types"))

# Now count num markerGenes on in each cell type. (p(on) >= 0.9)
   subMat <- fullMat
   numCells <- apply(subMat >= 0.9,1,sum)

   geneID2name<- unique(data$specs$transcriptInfo[,c("gene_id", "gene_name")])
   numCells.geneID <- numCells
   names(numCells.geneID) <- geneID2name[match(names(numCells), geneID2name$gene_name),"gene_id"]

   outFn <- paste0(data$specs$outDir,"/",figName.histo,"_exprBreadthHisto.pdf")

   pdf(outFn, height=2,width=2.5)
   par(mar=c(3,4.0,0.5,0.5),ps=10,xpd=FALSE,pty="s")
   raw.hist <- hist(numCells, plot=FALSE,breaks=50)
   plot(raw.hist,main="",xlab="",ylab="",
        col="grey60", cex.axis=0.8,border=NA,las=1,yaxt="n",xaxt="n")

   axis(1, cex.axis=0.8, mgp=c(3,0.5,0))
   axis(2, cex.axis=0.8, mgp=c(3,0.5,0),las=1)

   mtext("# genes", 2, line=2.0)
   mtext("Number of expressing cells", 1, line=1.5)

   if (0) {
   for (gene in c("fkh", "Ets65A")) {
      abline(v=numCells[[gene]],lwd=2,col="black")
      text(numCells[[gene]], max(raw.hist$counts) / 2,
           paste0(gene,": ", numCells[[gene]]),
           col="black",
           adj=c(0,0), offset=0.2)
   }
   }

   dev.off()

# Draw CDF for all genes, overlay TF subset, and homeodomains
#   either (SUPERFAMILY IPR009057 or PFAM-based IPR001356)

   goInfo <- loadGO(data$specs)

   geneSubsets <- list(
#      TF = data$specs$transcriptInfo$gene_name[
#                     data$specs$transcriptInfo$gene_id %in% loadTFs(dat)[["tf"]]],
      TF = loadTFs.flybase(dat),
      homeodomain = unique(data$specs$transcriptInfo$gene_name[
                     data$specs$transcriptInfo$gene_id %in%
                     data$specs$interproInfo$gene_id[
                        data$specs$interproInfo$domain == "IPR001356"]]),
      veryspecific = unique(data$specs$transcriptInfo$gene_name[
                        data$specs$transcriptInfo$gene_id %in%
                        goInfo$gene_id[goInfo$go_term ==
                        "neuropeptide hormone activity"]]),
      cns.morpho = unique(data$specs$transcriptInfo$gene_name[
                        data$specs$transcriptInfo$gene_id %in%
                        goInfo$gene_id[goInfo$go_term ==
                        "central nervous system morphogenesis"]]),
      verybroad = unique(data$specs$transcriptInfo$gene_name[
                              data$specs$transcriptInfo$gene_id %in%
                              goInfo$gene_id[goInfo$go_term ==
                              "synaptic vesicle endocytosis"]])
   )

   outFn <- paste0(data$specs$outDir,"/",figName.cdf,"_exprbreadth_cdf_tf.pdf")
   pdf(outFn, height=2, width=2.5)
   par(mar=c(3,4.0,0.5,0.5),ps=10,xpd=FALSE,pty="s")
   plot(ecdf(numCells[geneSubsets$TF]),pch=20,cex=0.2,col="black",
         col.01line = NA, verticals=T,
         main="",
         xlab="",ylab="",
         xaxt="n", yaxt="n"
         )
   plot(ecdf(numCells[geneSubsets$homeodomain]),pch=20,cex=0.2,col="orange",
        add=TRUE, col.01line = NA, verticals=T)
   plot(ecdf(numCells),pch=20,cex=0.2,col="gray60",add=TRUE, col.01line = NA,
        verticals=T)
   plot(ecdf(numCells[geneSubsets$verybroad]),
        pch=20,cex=0.2,col="steelBlue4",add=TRUE, col.01line = NA,
        verticals=T)
   plot(ecdf(numCells[geneSubsets$veryspecific]),
        pch=20,cex=0.2,col="red",add=TRUE, col.01line = NA,
        verticals=T)

   mtext("Number of expressing cells",1, line=1.5)
   mtext("Cumulative fraction", 2, line=1.8)
   axis(1, cex.axis=0.8, mgp=c(3,0.5,0), las=1)
   axis(2, cex.axis=0.8, las=1)

   legend("top",
          legend=c("all","TF","homeodomain", "s. vesicle endocytosis", "neuropeptide"),
          text.col=c("gray60","black","orange","steelBlue4","red"),
          cex=0.5,
          bty="n")

   dev.off()



# NEW SECTION: FlyBase Gene Groups 180428_1440

   print("LOADING FLYBASE GROUP ASSIGNMENTS") 
   fbggInfo <- loadFlyBaseGeneGroups(fn = data$specs$fullFlyBaseGroupsFn,
      terminalGroupsOnly = TRUE)

# RESTRICT to genes that exist in the numCells.geneID
   print(paste0("-> groups with at least 10 members: ",
                sum(table(fbggInfo$FB_group_id) >= 10)))

   print(paste0("-> Starting with group assignemnts for n=",
         length(unique(fbggInfo$Group_member_FB_gene_id))," genes"))

   fbggInfo <- fbggInfo[fbggInfo$Group_member_FB_gene_id %in%
                        names(numCells.geneID),]
   print(paste0("-> Restricting to those with expr pMat eentries, n=",
         length(unique(fbggInfo$Group_member_FB_gene_id))," genes"))

   print(paste0("-> now, groups with at least 10 membres: ",
                sum(table(fbggInfo$FB_group_id) >= 10)))


   fbgg2medianExprCells <- vector()
   fbgg2nExprCells <- list()
   for (fbgg in sort(unique(fbggInfo$FB_group_id))) {
      fbggName <- head(fbggInfo$FB_group_name[fbggInfo$FB_group_id == fbgg], n=1)

      if (nrow(fbggInfo[fbggInfo$FB_group_id == fbgg,]) < 10) {next;}

      fbgg2medianExprCells[fbggName] <- median( na.omit(numCells.geneID[
         fbggInfo$Group_member_FB_gene_id[fbggInfo$FB_group_id == fbgg]]))

      fbgg2nExprCells[[fbggName]] <- na.omit(numCells.geneID[
         fbggInfo$Group_member_FB_gene_id[ fbggInfo$FB_group_id == fbgg]])

      if (fbgg == "FBgg0000581") {print("YUP GOT HERE ON ION CHANNELS!")}
   }
   print(paste0("PROCESSED n=",length(fbgg2medianExprCells)," groups with at least 10 members"))

   fbgg2medianExprCells <- sort(fbgg2medianExprCells)
   print(paste0("-> post-sort n= ",length(fbgg2medianExprCells)))

   outFn <- paste0(data$specs$outDir,"/",figName.genegroups,
                   "_exprbreadth_medians.pdf")
   pdf(outFn, height=4, width=2.0)
   par(mar=c(4,2.5,0.5,0.5),ps=10,xpd=FALSE)
   xRange <- c(0, max(numCells.geneID))
   yRange <- c(0, length(fbgg2medianExprCells))
   plot(xRange, yRange,
        type="n",ylab="",axes=FALSE,xlab="Number of expressing cells")
   axis(1)

   nRows <- length(fbgg2medianExprCells)
   ggList <- names(fbgg2medianExprCells)

   mtext(side=2,text=paste0("FlyBase gene groups (n=",nRows,")"),line=0.5)

   hiGroups <- list(
      "DPR-INTERACTING PROTEINS" = "#1f78b4",
      "BEAT FAMILY" = "#7b3294",
      "DEFECTIVE PROBOSCIS EXTENSION RESPONSE" = "#e66101"
   )

   for (i in 1:nRows) {
      curX <- fbgg2nExprCells[[ggList[i]]]
      curY <- rep(nRows - i, length(curX))
      points(curX,curY,cex=0.5, col="gray70",pch=20)
      points(median(curX), curY[1],cex=1,col="black",pch=20)
   }

   for (i in 1:nRows) {
      if (! ggList[i] %in% names(hiGroups)) {next;}
      curColor <- hiGroups[[ggList[i]]]
      curX <- fbgg2nExprCells[[ggList[i]]]
      curY <- rep(nRows - i, length(curX))
      points(curX,curY,cex=0.5, col=curColor,pch=20)
      points(median(curX), curY[1],cex=1,col=curColor,bg="black",pch=21)
   }

   dev.off()
   return(list(genegroups = fbgg2medianExprCells, fbgg2nExprCells = fbgg2nExprCells))


   uniqGOterms <- sort(unique(goInfo$go_term))
   fSpecificGO <- vector(length=length(uniqGOterms));
   names(fSpecificGO) <- uniqGOterms
   numGOgenes <- fSpecificGO
   for (goTerm in names(fSpecificGO)) {
      curGeneID <- unique(goInfo$gene_id[goInfo$go_term == goTerm])
      curGeneName <- unique(data$specs$transcriptInfo$gene_name[data$specs$transcriptInfo$gene_id %in% curGeneID])
      curGeneName <- intersect(curGeneName,names(numCells))
      nCurGene <- length(numCells[curGeneName])
      fSpecificGO[goTerm] <- sum(numCells[curGeneName] <= 30) / nCurGene
      numGOgenes[goTerm] <- nCurGene
   }


   uniqIPR <- sort(unique(data$specs$interproInfo$domain))
   fSpecificIP <- vector(length=length(uniqIPR));
   names(fSpecificIP) <- uniqIPR
   numIPgenes <- fSpecificIP
   for(ipr in names(fSpecificIP)) {
      curGeneID <- unique(data$specs$interproInfo$gene_id[
                           data$specs$interproInfo$domain == ipr])
      curGeneName <- unique(data$specs$transcriptInfo$gene_name[data$specs$transcriptInfo$gene_id %in% curGeneID])
      if (ipr == "IPR001356") {print("GOT IN on IPR001356!!"); print(curGeneName)}

      curGeneName <- intersect(curGeneName,names(numCells))
      if (ipr == "IPR001356") {print("POST-name check"); print(curGeneName)}

      nCurGene <- length(numCells[curGeneName])
      fSpecificIP[ipr] <- sum(numCells[curGeneName] <= 30) / nCurGene
      numIPgenes[ipr] <- nCurGene
   }

   return(list(genegroups = fbgg2medianExprCells,
               go = cbind(numGOgenes, fSpecificGO),
               interpro = cbind(numIPgenes, fSpecificIP)))

}


prepRawPonMats <- function(data, mode="cells") {
# Purpose: prepare matched matrices of sample-, driver-, or cell- level TPM
#          averages and probability calls

   geneList <- rownames(data$pFit$logp_on_gs)

   rawMat <- dat$expr$geneExpr[
               match(geneList, dat$expr$geneExpr$gene_name),
               setdiff(colnames(dat$pFit$inDat$logC), c("tpm.min","tpm.max")),
               drop=FALSE]

   rownames(rawMat) <- geneList
   sampleNames <- gsub("^tpm.", "", colnames(rawMat))

   if ( mode == "samples" ) {

      pMat <- data$pFit$logp_on_gs

   } else if ( mode == "drivers" ) {

      t.curDrivers <- setdiff(sort(unique(dat$pFit$inDat$driver)), c("min","max"))
      tpmCols <- t.curDrivers
      pMat <- data$pFit$logp_on_gd[, tpmCols]

      t.gdMat <- matrix(nrow=nrow(rawMat),
                        ncol=length(t.curDrivers),
                        dimnames = list( rownames(rawMat), t.curDrivers ))

      for (t.curDriver in t.curDrivers) {
         t.curSamples <- sampleNames[ dat$pFit$inDat$driver == t.curDriver]
         t.gdMat[,t.curDriver] <- apply(rawMat[, paste0("tpm.",t.curSamples) ],
                                        1, mean)
      }
      rawMat <- as.data.frame(t.gdMat)

   } else if ( mode == "cells" ) {

      pMat <- data$pFit$logp_on_gc

      t.curCells <- setdiff(sort(unique(dat$pFit$inDat$cell)), c("min","max"))
      tpmCols <- t.curCells
      pMat <- data$pFit$logp_on_gc[, tpmCols]

      t.gcMat <- matrix(nrow=nrow(rawMat),
                        ncol=length(t.curCells),
                        dimnames = list( rownames(rawMat), t.curCells ))

      for (t.curCell in t.curCells) {
         t.curSamples <- sampleNames[ dat$pFit$inDat$cell == t.curCell]
         t.gcMat[,t.curCell] <- apply(rawMat[, paste0("tpm.",t.curSamples) ],
                                      1, mean)
      }
      rawMat <- as.data.frame(t.gcMat)

   }

   return(list( raw = rawMat, p = pMat))

}


plotOnRangeHisto <- function(data,
                             figName="msFig5F",
                             labelGenes) {
# Purpose: make histogram of dynamic range of on expression state

   pMat <- exp(data$expr$rawpon$cells$p)
   rawMat <- data$expr$rawpon$cells$raw

   onRange <- vector(length=nrow(pMat))
   names(onRange) <- rownames(pMat)
   onRange <- NA
   for (gene in rownames(pMat)) {
      ttx <- pMat[gene,] >= 0.9
      if (sum(ttx) < 2) next;
      t.range <- range(rawMat[gene, ttx])
      if (t.range[2] < 10) next
      onRange[gene] <- (t.range[2] + 1) / (t.range[1] + 1)
   }

# for each gene, calculate min-max range for cells with pMat >= 0.9

   onRange <- na.omit(onRange)
   rawHist <- hist( log(onRange), plot=FALSE,breaks=50 )

   outFn <- paste0(data$specs$outDir,"/",figName,"_onRangeHisto.pdf")

   pdf(outFn, height=2,width=2.5)
   par(mar=c(3,3.5,0.5,1),ps=10)
   plot(rawHist, main="",xlab="",ylab="",
        col="grey60", cex.axis=0.8,border=NA,las=1,yaxt="n",xaxt="n")

   fcLabs <- c(1,2,5,10,20,50,100,200,500,1000)

   axis(1,at=log(fcLabs),labels=fcLabs,
                 cex.axis=0.8, mgp=c(3,0.5,0))
   axis(2,cex.axis=0.8,mgp=c(3,0.5,0),las=1)

   if (!missing(labelGenes)) {
   for (gene in labelGenes) {
      abline(v=log(onRange[[gene]]),lwd=2,col="black")
      text(log(onRange[[gene]]), max(rawHist$counts) / 2,
           paste0(gene,": ", round(onRange[[gene]],digits=0)),
           col="black",
           adj=c(0,0), offset=0.2)
   }
   }

   mtext("# genes", 2, line=2.0)
   mtext("On-state dynamic range\n(TPMmax/min)", 1, line=2)

   dev.off()

}


plotMuOnOffHistos <- function(data,
                              figName.muOnHisto = "msFig5F.new1",
                              figName.muOffHisto = "msFig5F.new2",
                              figName.dmuOnOffHisto ="msFig5F.new3",
                              labelGenes
                              ) {
# Purpose: make histograms of expression state mean levels (on, off, and don-off)

   allGenes <- rownames(data$pFit$logp_on_gc)
   print("Ignoring CR genes")
   allGenes <- allGenes[!grepl("^CR",allGenes)]
   print(paste0(length(allGenes)," coding genes"))

# Only makes sense for genes that are bimodal -- quick list of these?

   isBimodal <- apply(dat$pFit$logp_on_gc[allGenes,], 1,
                      function(x) {min(x) != max(x)})
   allGenes <- names(isBimodal[isBimodal == TRUE])

# Good-sample Cell type

   muons <- vector()
   muoffs  <- vector()
   for (gene in allGenes) {
      muons[gene] <- dat$pFit$geneEfit[[gene]]$bimPars$mu_on
      muoffs[gene] <- dat$pFit$geneEfit[[gene]]$bimPars$mu_off
   }
   dmu <- exp(muons) / exp(muoffs)

   outFn <- paste0(data$specs$outDir,"/",figName.muOnHisto,"_muOnHisto.pdf")


   rawHist <- hist(muons, plot=FALSE, breaks=50)
   pdf(outFn, height=2,width=2.5)
   par(mar=c(3,3.5,0.5,1),ps=10)
   plot(rawHist, main="",xlab="",ylab="",
        col="grey60", cex.axis=0.8,border=NA,las=1,yaxt="n",xaxt="n")

   tpmLabs <- c(1,10,100,1000,10000)

   axis(1,at=log(tpmLabs + 1),labels=tpmLabs,
                 cex.axis=0.8, mgp=c(3,0.5,0))
   axis(2,cex.axis=0.8,mgp=c(3,0.5,0),las=1)

   if (!missing(labelGenes)) {
   for (gene in labelGenes) {
      abline(v=muons[[gene]],lwd=2,col="black")
      text(muons[[gene]], max(rawHist$counts) / 2,
           paste0(gene,": ", round(exp(muons[[gene]]) - 1, digits=0)),
           col="black",
        adj=c(0,0), offset=0.2)
   }
   }

   mtext("# genes", 2, line=2.0)
   mtext("On-state abundance\n(TPM + 1)", 1, line=2)

   dev.off()



   outFn <- paste0(data$specs$outDir,"/",figName.muOffHisto,"_muOffHisto.pdf")
   rawHist <- hist(muoffs, plot=FALSE, breaks=50)
   pdf(outFn, height=2,width=2.5)
   par(mar=c(3,3.5,0.5,1),ps=10)
   plot(rawHist, main="",xlab="",ylab="",
        col="grey60", cex.axis=0.8,border=NA,las=1,yaxt="n",xaxt="n")

   tpmLabs <- c(1,10,100,1000,10000)

   axis(1,at=log(tpmLabs + 1),labels=tpmLabs,
                 cex.axis=0.8, mgp=c(3,0.5,0))
   axis(2,cex.axis=0.8,mgp=c(3,0.5,0),las=1)

   if (!missing(labelGenes)) {
   for (gene in labelGenes) {
      abline(v=muoffs[[gene]],lwd=2,col="black")
      text(muoffs[[gene]], max(rawHist$counts) / 2,
        paste0(gene,": ", round(exp(muoffs[[gene]]) - 1, digits=0)),
        col="black",
        adj=c(0,0), offset=0.2)
   }
   }

   mtext("# genes", 2, line=2.0)
   mtext("Off-state abundance\n(TPM + 1)", 1, line=2)

   dev.off()


   outFn <- paste0(data$specs$outDir,"/",figName.dmuOnOffHisto,"_muDiffOnOffHisto.pdf")
   rawHist <- hist(log(dmu), plot=FALSE, breaks=50)
   pdf(outFn, height=2,width=2.5)
   par(mar=c(3,3.5,0.5,1),ps=10)
   plot(rawHist, main="",xlab="",ylab="",
        col="grey60", cex.axis=0.8,border=NA,las=1,yaxt="n",xaxt="n")

   fcLabs <- c(1,2,5,10,20,50,100,200,500,1000)

   axis(1,at=log(fcLabs),labels=fcLabs, cex.axis=0.8, mgp=c(3,0.5,0))
   axis(2,cex.axis=0.8,mgp=c(3,0.5,0),las=1)

   if (!missing(labelGenes)) {
   for (gene in labelGenes) {
      abline(v=log(dmu[[gene]]),lwd=2,col="black")
      text(log(dmu[[gene]]), max(rawHist$counts) / 2,
         paste0(gene, ": ", round(dmu[[gene]], digits=0)),
         col="black",
         adj=c(0,0), offset=0.2)
   }
   }

   mtext("# genes", 2, line=2.0)
   mtext("On/Off dynamic range\n(TPM/TPM)", 1, line=2)

   dev.off()

}



callNTout <- function(dat) {
# Purpose: For each modelled celltype, call NT output using simple rules

   ponThresh <- 0.9

   ntSig <- list(
                 'Hist' = c("Hdc"),
                 'GABA' = c("Gad1"),
                 'ACh'  = c("VAChT"),
                 'Glu'  = c("VGlut"),
                 'DA'   = c("ple", "Vmat")
                 )

   pMat <- dat$pFit$logp_on_gc[,
      setdiff(colnames(dat$pFit$logp_on_gc), c("min","max"))]

   ntCalls <- vector(mode="character",length=ncol(pMat))
   names(ntCalls) <- colnames(pMat)

   maxNumNT <- 0
   for (i in 1:length(ntCalls)) {
      curCell <- names(ntCalls)[i]
      curCall <- c()
      for (nt in names(ntSig)) {
         if (all(exp(pMat[ntSig[[nt]],curCell]) >= ponThresh)){
            curCall <- c(curCall, nt)
         }
      }
      if (length(curCall) == 0) { curCall <- "unk"}
      if (length(curCall) > maxNumNT) {maxNumNT <- length(curCall)}
      ntCalls[curCell] <- paste(curCall,collapse=".")
   }

   ntMat <- matrix(nrow=length(ntCalls), ncol=maxNumNT, data = "unk")
   rownames(ntMat) <- names(ntCalls)
   print(paste0("Max # NT calls: ",maxNumNT))

   for (i in 1:length(ntCalls)) {
      curNT <- unlist(strsplit(ntCalls[i], "\\."))
      for (j in 1:length(curNT)) {
         ntMat[names(ntCalls)[i], j] <- curNT[j]
      }
   }
   colnames(ntMat) <- paste0("NTcall",1:ncol(ntMat))

   return(list( ntString        = ntCalls,
                ntMat           = ntMat))

}


makeCoTree <- function(data,
                       figName.cotree = "cotree_allgenes_tfgenes",
                       figName.allgenes = "tree_allgenes",
                       figName.TFgenes = "tree_TFgenes") {
# Purpose: make 'phylogenetic' tree of cell types based on expression pattern
#   of either all genes or TF genes

   tree.all <- buildToggleTree(data,
                               figName = figName.allgenes,
                               treeType = "dna.fastme",
                               treeGeneType     ="all",
                               hiSNR.only = TRUE,
                               makePlots = TRUE,
                               nCores   = 4)

   tree.all$tr$node.label <- tree.all$edgeSupp

   tree.tf <- buildToggleTree(data,
                              figName = figName.TFgenes,
                              treeType  ="dna.fastme",
                              makePlots = TRUE,
                              treeGeneType="TF",
                              nCores    = 4)
   tree.tf$tr$node.label <- tree.tf$edgeSupp

   coTree <- cophylo(tree.all$tr,tree.tf$tr)

   pdf(paste0(dat$specs$outDir,"/",figName.cotree,"_cophylo_allgene_TFgenes.pdf"))
   plot(coTree,
        fsize=0.4, ftype="reg",
        link.type='curved',
        link.lwd=1,link.lty='solid')

# match edges in new rotated tree to old one.

   nodelabels.cophylo(pie=cbind(100 * as.numeric(coTree$trees[[1]]$node.label),
                                100-(100 * as.numeric(coTree$trees[[1]]$node.label))),
                      piecol=c("red","white"),cex=0.25)

   nodelabels.cophylo(pie=cbind(100 * as.numeric(coTree$trees[[2]]$node.label),
                                100-(100 * as.numeric(coTree$trees[[2]]$node.label))),
                      piecol=c("red","white"),cex=0.25, which="right")

   dev.off()

}


makeCoHeatmap <- function(dat, figName) {
# Purpose: Draw paired heatmaps for TPM and p(z_gd=on)

   t.gdMat <- dat$expr$rawpon$drivers$raw
   colnames(t.gdMat) <- paste0("tpm.",colnames(t.gdMat))

   print("Heatmapping p(on) matrix")
   hmap.pon <- plotHmapPmat( dat,
                             figName            = paste0(figName,"_bottom"),
                             mode               = "gd",
                             legend             = FALSE,
                             plotTitle          = "",
                             treeheight_row     = 0,
                             treeheight_col     = 0,
                             fontsize_row       = 0,
                             fontsize_col       = 0,
                             outFormat          = "png",
                             plotWidth          = 2.5,
                             plotHeight         = 2.5 )

# reordering t.gdMat to match p(on) heatmap
   t.gdMat <- t.gdMat[hmap.pon$tree_row$labels[hmap.pon$tree_row$order],
                      paste0("tpm.",hmap.pon$tree_col$labels[hmap.pon$tree_col$order])]

   print("Heatmapping relTPM matrix")
   hmap.tpm <- plotHeatmap( dat,
                            tpmMat              = t.gdMat,
                            legend              = FALSE, main="",
                            as.is               = TRUE,
                            cluster_rows        = FALSE,
                            cluster_cols        = FALSE,
                            geneList            = list(all = rownames(t.gdMat)),
                            treeheight_row      = 0,
                            treeheight_col      = 0,
                            fontsize_col        = 2,
                            show_rownames       = FALSE,
                            fontsize_row        = 0,
                            show_colnames       = FALSE,
                            outFormat           = "png",
                            plotWidth           = 2.5,
                            plotHeight          = 2.5,
                            filePrefix          = paste0(figName,"_top") )

}



loadRiveraalba11 <- function(dat) {
# Purpose: Load lamina connectome from Rivera-Alba et al., 2011

   library(readxl)

   print(paste0("Loading connectome: RiveraAlba 11 :",
                dat$specs$riveraalba11$fn))
   riveraalba11 <- read.table(dat$specs$riveraalba11$conMatrixFn,
                              header    = TRUE,
                              sep       = "\t",
                              comment   = "#",
                              row.names = 1,
                              as.is     = TRUE)

# RiveraAlba: column=pre-synaptic cell, row=post-synaptic cell
   riveraalba11 <- as.matrix(t(riveraalba11))

# 1. evaluate sums (certain + uncertain synapses)
   riveraalba11 <- apply( riveraalba11,
                          c(1,2),
                          function(x) {eval(parse(text=x))})

# synonyms: "La wf" -> Lawf1, Lawf2 = Lawf
   newcol <- riveraalba11[, "La wf"]
   riveraalba11 <- cbind(riveraalba11, Lawf1=newcol)
   riveraalba11 <- cbind(riveraalba11, Lawf2=newcol)

   newcol <- riveraalba11[, "Am"]
   riveraalba11 <- cbind(riveraalba11, Lai=newcol)

   newcol <- riveraalba11[, "Ep.glia"]
   riveraalba11 <- cbind(riveraalba11, Glia_Eg=newcol)

   newcol <- riveraalba11[, "S.glia"]
   riveraalba11 <- cbind(riveraalba11, Glia_Psg=newcol)

   newcol <- riveraalba11[, "M.glia"]
   riveraalba11 <- cbind(riveraalba11, Glia_Mg=newcol)

   riveraalba11 <- as.data.frame(riveraalba11)
   delCells  <- c( "La wf", "Am", "Ep.glia", "S.glia", "M.glia",
                   "L4+x", "L4-y", "Orph.", "Total" )

   for (delCell in delCells) {
      riveraalba11[, delCell] <- NULL }

   riveraalba11 <- as.matrix(riveraalba11)

   riveraalba11 <- rbind(riveraalba11, Lawf1=riveraalba11["La.wf", ])
   riveraalba11 <- rbind(riveraalba11, Lawf2=riveraalba11["La.wf", ])
   riveraalba11 <- riveraalba11[!rownames(riveraalba11) %in%
                                c("La.wf", "L4.x", "L4.y", "Orph", "Total"), ]

   riveraalba11 <- rbind(riveraalba11, Lai=riveraalba11["Am", ])
   riveraalba11 <- riveraalba11[!rownames(riveraalba11) %in% c("Am"), ]


   curRowNames <- rownames(riveraalba11)
   curRowNames[curRowNames == "Ep.glia"] <- "Glia_Eg"
   curRowNames[curRowNames == "S.glia"] <- "Glia_Psg"
   curRowNames[curRowNames == "M.glia"] <- "Glia_Mg"
   rownames(riveraalba11) <- curRowNames

   riveraalba11 <- rbind("R1-6"=colSums(riveraalba11[paste0("R",1:6),]), riveraalba11)
   riveraalba11 <- riveraalba11[!rownames(riveraalba11) %in% paste0("R",1:6),]

   riveraalba11 <- cbind("R1-6" = rowSums(riveraalba11[,paste0("R",1:6)]), riveraalba11)
   riveraalba11 <- riveraalba11[,!colnames(riveraalba11) %in% paste0("R",1:6)]

   return(riveraalba11)
}

loadTakemura13 <- function(dat) {
# Purpose: Load connectome tables from Takemura et al., 2013

   library(readxl)
   print(dat$specs$takemura13$synapsesFn)
   tx <- data.frame(read_xls(dat$specs$takemura13$synapsesFn,
                             col_names=FALSE))
   synPairs <- tx[- c(1:3,29:36),]
   colnames(synPairs) <- c("Postsynaptic.Sites",
                           "Presynaptic.Name",
                           "Presynaptic.x",
                           "Presynaptic.y",
                           "Presynaptic.z",
                           "Presynaptic.comments",
                           "Postsynaptic.Name",
                           "Postsynaptic.x",
                           "Postsynaptic.y",
                           "Postsynaptic.z",
                           "Postsynaptic.comments",
                           "Proofreading.details")

   synPairs$Presynaptic.Name  <- gsub("LaWF", "Lawf", synPairs$Presynaptic.Name)
   synPairs$Postsynaptic.Name  <- gsub("LaWF", "Lawf", synPairs$Postsynaptic.Name)

   synPairs$Presynaptic.id <- apply(synPairs[,c("Presynaptic.Name",
                                                "Presynaptic.x",
                                                "Presynaptic.y",
                                                "Presynaptic.z")], 1,
                                    paste,collapse=".")
   return(list(
      synPairs = synPairs
   ))

}

getTakemura13Partners <- function(dat,
                                  cellTypes,
                                  minSyn = 5) {

   takemura13 <- loadTakemura13(dat)

   curDat <- takemura13$synPairs
   curDat$Postsynaptic.celltype <- gsub(" .*", "", curDat$Postsynaptic.Name)
   curDat$Postsynaptic.celltype[curDat$Postsynaptic.celltype == ""] <- NA

   partners <- list()
   for (cellType in cellTypes) {
      curPart <- table(curDat$Postsynaptic.celltype[curDat$Presynaptic.Name == cellType])
      partners[[cellType]] <- sort(curPart[curPart >= minSyn])
   }
   return(partners)

}

plotSynapseReceptors <- function(dat,
                                 presynaptic.name="R8 111",
                                 genes = c("ort", "HisCl1"),
                                 figName = "postSynExpr",
                                 plot.width = 1.5,
                                 plot.height = 3.5,
                                 plot.ptCex = 0.5,
                                 plot.labCex = 0.8
                                 ) {

# Purpose: plots active zones for a query neuron and colors by gene expression
#  Uses data in Takemura et al., 2013

   library(RColorBrewer)
   colorPal <- c(brewer.pal(2 ** length(genes), name="Set1"), "darkgray")

# Load synapse informatino
   takemura13 <- loadTakemura13(dat)

   curDat <- takemura13$synPairs[takemura13$synPairs$Presynaptic.Name == presynaptic.name,]
   curDat$Postsynaptic.celltype <- gsub(" .*", "", curDat$Postsynaptic.Name)
   curDat$Postsynaptic.celltype[curDat$Postsynaptic.celltype == ""] <- NA

   uniqCellTypes <- na.omit(unique(curDat$Postsynaptic.celltype))

# Get individual gene p(expr)
   geneCalls <- matrix(nrow=length(genes), ncol=length(uniqCellTypes),
                       dimnames=list(genes,uniqCellTypes),data=NA)
   for( cell in uniqCellTypes) {
      if (cell == "R8") {
         subCells <- c("R8_Rh5", "R8_Rh6")
      } else if (!cell %in% colnames(dat$expr$rawpon$cells$p)){
         print(paste0("DONT HAVE EXPR FOR ",cell))
         next;
      } else {
         subCells <- c(cell)
      }
      print(paste0("subCells == ",subCells))

      if(length(subCells) == 1) {
         geneCalls[genes,cell] <- exp(dat$expr$rawpon$cells$p[genes,subCells])
      } else {
         geneCalls[genes,cell] <- apply(exp(dat$expr$rawpon$cells$p[genes,subCells]),1,max)
      }
   }

# Discretize
   discCalls <- geneCalls
   discCalls[geneCalls >= 0.8] <- 1
   discCalls[geneCalls < 0.8] <- 0
   geneCalls <- discCalls

# Combine into gene pattern for each cell
   cell2expr <- list()
   for (cell in colnames(geneCalls)) {
      if (any(is.na(geneCalls[,cell]))) {
         cell2expr[[cell]] <- NA
      } else {
         cell2expr[[cell]] <- paste0(geneCalls[,cell],collapse="")
      }
   }

# Name each category; positive/negative naems
   typeNames <- c()
   names(typeNames)
  

# Output: matrix where row = presynaptic identifier, col = gene on/off combo,
#   value is number of post-synapses with that gene signature

   synTypeMat <- matrix(nrow = length(unique(curDat$Presynaptic.id)),
                    ncol = 1 + 2 ** length(genes),
                    data = 0)
   {
      nGenes <- length(genes)
      tx <- list()
      for (i in 1:nGenes) { tx <- c(tx, list(c(0,1))) }
      geneCombo <- apply(expand.grid(tx),1,function(x){paste0(x, collapse="")})
   }

   rownames(synTypeMat) <- unique(curDat$Presynaptic.id)
   colnames(synTypeMat) <- c(geneCombo,"unk")
   names(colorPal) <- colnames(synTypeMat)

   {
      tx <- list()
      for (i in 1:nGenes) { tx <- c(tx, list(c(paste0(genes[i],"-"),
                                               paste0(genes[i],"+")))) }
      geneComboNames <- c(apply(expand.grid(tx),1,function(x){paste0(x, collapse=" ")}),"unknown")
   }
   geneComboNames <- unlist(geneComboNames)
   if (length(genes) == 2) {
      geneComboNames[!grepl("\\+",geneComboNames) &
                      grepl("\\-",geneComboNames) ] <- "neither"
   } else {
      geneComboNames[!grepl("\\+",geneComboNames) &
                      grepl("\\-",geneComboNames) ] <- "none"
   }
   print(geneComboNames)

# Iterate over pre/postsynapse pairs and tally up receptor types
   for (i in 1:nrow(curDat)) {
      if ( is.na(curDat$Postsynaptic.celltype[i]) |
           is.na(cell2expr[curDat$Postsynaptic.celltype[i]]) ) {

         synTypeMat[curDat$Presynaptic.id[i], "unk"] <- 
            1 + synTypeMat[curDat$Presynaptic.id[i], "unk"]

      } else {

         curPostExpr <- cell2expr[[curDat$Postsynaptic.celltype[i]]]

         synTypeMat[curDat$Presynaptic.id[i], curPostExpr] <- 
            1 + synTypeMat[curDat$Presynaptic.id[i], curPostExpr]

      }
   }


# figure out plot width
   width.numSyn <- sum(apply(synTypeMat,2,max))
   height.numSyn <- nrow(synTypeMat)

# Multiplier for x / y position
   xIncr <- 1/width.numSyn
   yIncr <- -1 / height.numSyn

# Baseline X position for each expression category; for labels and data points
   typeLabel.x <- cumsum(apply(synTypeMat,2,max))
   typeLabel.x <- c(0,typeLabel.x[1:(length(typeLabel.x) - 1)])
   typeLabel.x <- typeLabel.x * xIncr

   outFn <- paste0(dat$specs$outDir,"/",figName,"_",
                   presynaptic.name,".",
                   paste(genes,collapse="_"),
                   ".pdf")
   outFn <- gsub(" ", "_",outFn)
   nSyn <- nrow(synTypeMat)
   pdf(outFn, width=plot.width, height=plot.height)
   par(mar=c(3,3,4,0.5))
   plot(c(0,1), c(-1,0), type="n", axes=FALSE)
   axis(side=2,labels=FALSE,line=0.25,tcl=0)

   mtext(side=1,line=1, text="post-synaptic\nexpression", cex=plot.labCex)
   mtext(side=2,line=0.7,text=paste0(presynaptic.name," active zones"), cex=plot.labCex)

   nonZeroTypes.ind <- which(apply(synTypeMat,2,max) > 0)
   nonZeroTypes <- colnames(synTypeMat)[nonZeroTypes.ind]
   print(paste0("nonzerotypes = ",nonZeroTypes))

   text(x=typeLabel.x[nonZeroTypes.ind],
        y=par("usr")[4],
        srt=60,xpd=TRUE,
        labels = geneComboNames[nonZeroTypes.ind],
        adj=c(0,0),
        col = colorPal[nonZeroTypes],
        cex=plot.labCex)

   mtext(side=2,at=0,text="lateral", line=0.5,las=1, cex=plot.labCex)
   mtext(side=2,at=-1,text="medial", line=0.5,las=1, cex=plot.labCex)

   for (i in 1:nrow(synTypeMat)) {
      yPos <- yIncr * i
      for (j in 1:ncol(synTypeMat)){
         curPattern <- colnames(synTypeMat)[j]
         curVal <- synTypeMat[i,curPattern]
         if (curVal  == 0) next;
         xVal <- typeLabel.x[j] + (c(1:curVal) - 1) * xIncr
         points(xVal, rep(yPos, length(xVal)),
                pch=20,
                col=colorPal[curPattern],
                cex=plot.ptCex)
      }
   }

   dev.off()

   return(list(
      takemura13 = takemura13,
      geneCalls =geneCalls,
      cell2expr = cell2expr,
      synTypeMat = synTypeMat
   ))

}


plotExprOzkanPairs <- function(dat,
                               figName = "numExprOzkanPairs",
                               cells = c("R1-6","L1","L2","L3","L4","L5", "Lai",
                                         "T1", "C2", "C3", "Lawf1", "Lawf2",
                                         "Glia_Eg", "Glia_Mg", "Glia_Psg")) {

# Purpose: count number of expressed interacting protein pairs between pairs of
#          Laminar neurons

   nCells <-  length(cells)

   ozkan13 <- loadOzkanPPI(dat$specs)
   intPairs <- unique(data.frame(symbol1 = c(ozkan13$genePairs$symbol1,
                                             ozkan13$genePairs$symbol2),
                                 symbol2 = c(ozkan13$genePairs$symbol2,
                                             ozkan13$genePairs$symbol1),
                                 stringsAsFactors=FALSE))

# Draw heatmap depicting number of protein pairs
   for (cell in cells) {
      curExpr1 <- exp(dat$expr$rawpon$cells$p[match(intPairs$symbol1,
                                           rownames(dat$expr$rawpon$cells$p)),
                                     cell])

      curExpr2 <- exp(dat$expr$rawpon$cells$p[match(intPairs$symbol2,
                                           rownames(dat$expr$rawpon$cells$p)),
                                     cell])

      curExpr1[curExpr1 >= 0.8] <- 1
      curExpr1[curExpr1 < 0.8] <- 0

      curExpr2[curExpr2 >= 0.8] <- 1
      curExpr2[curExpr2 < 0.8] <- 0

      intPairs[,paste0("expr1.",cell)] <- curExpr1
      intPairs[,paste0("expr2.",cell)] <- curExpr2
   }
   intPairs <-na.omit(intPairs)

   ppMat <- matrix(nrow=nCells, ncol=nCells, data=0, dimnames=list(cells,cells))
   for (cell1 in cells) {
      for (cell2 in cells) {
         ppMat[cell1, cell2] <- sum(intPairs[,paste0("expr1.",cell1)] == 1 &
                                    intPairs[,paste0("expr2.",cell2)] == 1)
      }
   }

   outFn <- paste0(dat$specs$outDir,"/",figName,"_heatmap.pdf")
   pdf(outFn,onefile=FALSE,height=5,width=6)

   library(ComplexHeatmap)
   library(circlize)
   library(RColorBrewer)
   ht1 <- ComplexHeatmap::Heatmap(ppMat,name="Extracellular\ninteracting\nprotein pairs",
                                  cluster_columns=FALSE,
                                  cluster_rows=FALSE,
                                  column_title = "Cell type 1",
                                  row_title = "Cell type 2",
      col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
   draw(ht1)
   dev.off()


   ccMat <- loadRiveraalba11(dat)
   ccMat <- ccMat[cells,cells]
   outFn <- paste0(dat$specs$outDir,"/",figName,"_connectome_heatmap.pdf")
   pdf(outFn,onefile=FALSE,height=5,width=6)
   ht1 <- ComplexHeatmap::Heatmap(log10(1+ccMat),name="Laminar\nSynapses\n(log10 1+n)",
                                  cluster_columns=FALSE,
                                  cluster_rows=FALSE,
                                  column_title = "Post-synaptic cell",
                                  row_title = "Pre-synaptic cell",
      col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
   draw(ht1)
   dev.off()
   


   return(ppMat)

}


compareToKonstantinides18 <- function( dat,
                                       figName = "figSX_scClusterMarkers",
                                       facsExpr) {
# Purpose: compare our data to optic lobe FACS-seq and scRNA-seq reported in
#   Konstantinides et al., 2018

# Read in cluster markers;

   if (0) {
   print(" Expression of single cell cluster makers in our dataset")
   scMarkers <- read.table(dat$specs$konstantinides18$markersFn,
                           header=TRUE,sep="\t",
                           stringsAsFactors=FALSE,
                           as.is=TRUE)
   t.1 <- scMarkers$gene_name
   print(paste0("Original markers table, n=",nrow(scMarkers)," genes"))
   scMarkers <- scMarkers[scMarkers$gene_name %in% dat$expr$geneExpr$gene_name,]
   t.2 <- scMarkers$gene_name
   print(paste0("-> in our matrix: n=",nrow(scMarkers)," genes"))
   print(paste0("Missing genes!:"))
   print(paste0(sort(setdiff(t.1,t.2)),collapse=", "))

   curCellGroups <- dat$specs$cellGroups[c(2,3,4)]

   geneList <- c()
   geneGaps <- c()
   geneGaps.labels <- c()
   for (curCluster in sort(unique(scMarkers$cluster_id))) {
      print(paste0("cluster = ",curCluster))
      curMarkers <- scMarkers$gene_name[scMarkers$cluster_id == curCluster]
      print(paste0(" -> ",length(curMarkers)," markers"))
      if (curMarkers == 0) {
         print(paste0("WARNING*************************** ",
                      "NO MARKERS FOR CLustER ",curCluster))}

      geneList <- c(geneList, curMarkers)
      geneGaps <- c(geneGaps, length(geneList))
      geneGaps.labels <- c(geneGaps.labels, curCluster)

   }
   geneGaps <- geneGaps[1:(length(geneGaps) - 1)]

   if (0) {
   plotHmapPmat(dat,
      figName = paste0(figName,"_pexpr"),
      mode = "gc",
      fontsize_col=2,
      plotWidth=12,
      plotHeight=7,
      legend=FALSE, plotTitle="",
      geneLabel = FALSE,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      gaps_genes = geneGaps,
      geneList = geneList,
      colType = "gene",
      cellGroups = curCellGroups
   )
   }

   plotHmapPmat(dat,
      figName = paste0(figName,"_relTPM"),
      mode = "gc",
      exprMode = "raw",
      rawNL = "mean",
      fontsize_col=4,
      plotWidth=8,
      plotHeight=5,
      legend=FALSE, plotTitle="",
      geneLabel = FALSE,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      gaps_genes = geneGaps,
      gaps_genes.group_labels = geneGaps.labels,
      geneList = geneList,
      colType = "gene",
      cellGroups = curCellGroups
   )
   }

# Compare FACS datasets
# ultimately; heatmap of our data on one axis, their samples on the other.

   if (1) {

   library(ComplexHeatmap)
   library(circlize)
   library(RColorBrewer)

   if (missing(facsExpr)) {
      facsExpr <- loadExpr.Konstantinides18(dat)
   }

# HERENOW  180616_1353
   our.tpmCols <- colnames(dat$expr$geneExpr.bycell)[grep("^tpm",
                     colnames(dat$expr$geneExpr.bycell))]
   facs.tpmCols <- colnames(facsExpr$geneExpr.bycell)[grep("^tpm",
                     colnames(facsExpr$geneExpr.bycell))]

   ourMat <- dat$expr$geneExpr.bycell[,our.tpmCols]
   rownames(ourMat) <- dat$expr$geneExpr.bycell$gene_name
   facsMat <- facsExpr$geneExpr.bycell[,facs.tpmCols]
   rownames(facsMat) <- facsExpr$geneExpr.bycell$gene_name

      if (0) {

         print(" GOT HERE")
         corMat <- cor(log1p(facsMat), log1p(ourMat))
         print(" GOT HERE")
         pdf(paste0(dat$specs$outDir,"/",figName,"_corMat_FACS_TAPIN.pdf"),
             height=5, width=20)
         print(" GOT HERE")
         ht1 <- ComplexHeatmap::Heatmap(corMat,
                  cluster_rows=FALSE,
                  cluster_columns=FALSE,
                  row_title = "FACS-seq",
                  col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                  column_title = "INTACT/TAPIN-seq")
         draw(ht1)
         dev.off()
         print(" GOT HERE")
      }

      if (1) {

         sharedCells <- sort(intersect(unique(facsExpr$sampleInfo$celltype),
                                       unique(dat$specs$sampleInfo$celltype)))


# Smart-Sort  to get Dm12 after Dm8
         {
            smartSplit.x <- gsub("[0-9]","",sharedCells)
            smartSplit.y <- as.numeric(gsub("[A-Za-z]","",sharedCells))
            sharedCells <- sharedCells[order(smartSplit.x, smartSplit.y)]
            print("NEW ORDER!:")
            print(paste(sharedCells, collapse=", "))
         }

         print("Shared cells:")
         print(paste(sharedCells, sep=", "))

         sharedCols <- paste0("tpm.",sharedCells)

         print("sharedCols")
         print(sharedCols)
# HERENOW 180620_0037

         ourMat <- ourMat[,sharedCols]
         facsMat <- facsMat[,sharedCols]

         ourMat  <- ourMat[ rownames(dat$expr$rawpon$cells$p),]
         facsMat <- facsMat[ rownames(dat$expr$rawpon$cells$p),]

      # MARKER STRAEGY: HIGHEST FC v MEAN
      # ALT: highest correlaton to bit vector pattern -- eg, 00001000
         fcMat <- facsMat[apply(facsMat,1,max) >= 50,]
         fcMat <- log2(1 + fcMat)
         fcMat <- fcMat - apply(fcMat,1,mean)
         facsMarkers <- list()
         geneGaps <- c()
         for (facsCol in colnames(facsMat)) {
            print(paste0("MARKER FOR ",facsCol))
            curMat <- fcMat[,facsCol,drop=FALSE]
            curMat <- curMat[order(curMat[,facsCol],decreasing=TRUE),,drop=FALSE]
            curMat <- curMat[curMat[,facsCol] > 2,,drop=FALSE]

#            curMat <- facsMat[facsMat[,facsCol] >= 50,]
#            curMat <- curMat[curMat[,facsCol] > 2 * 
#                              apply(curMat[,setdiff(colnames(curMat),facsCol)],1,max),]
#            curMat <- curMat[order(curMat[,facsCol],decreasing=TRUE),,drop=FALSE]

            facsMarkers[[facsCol]] <- head(rownames(curMat),n=10)

            if (length(geneGaps) == 0) {
               geneGaps <- length(facsMarkers[[facsCol]])
            } else {
               geneGaps <- c(geneGaps, geneGaps[length(geneGaps)] +
                                       length(facsMarkers[[facsCol]]))
            }
            print(paste0(facsMarkers[[facsCol]], collapse=", "))
         }
         geneGaps.labels <- gsub("tpm.","",colnames(facsMat))


   if (0) {
   plotHmapPmat(dat,
      figName = paste0(figName,"_FACSmarkers_pexpr"),
      mode = "gc",
      fontsize_col=2,
      plotWidth=10,
      plotHeight=5,
      legend=FALSE, plotTitle="",
      geneLabel = FALSE,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      gaps_genes = geneGaps,
      geneList = unlist(facsMarkers),
      colType = "gene",
      cellGroups = list(sharedCells)
#      cellGroups = dat$specs$cellGroups
   )
   }

   plotHmapPmat(dat,
      figName = paste0(figName,"_FACSmarkers_relTPM"),
      mode = "gc",
      exprMode = "raw",
      rawNL = "mean",
      fontsize_col=5,
      plotWidth=5,
      plotHeight=2.5,
      legend=FALSE, plotTitle="",
      geneLabel = FALSE,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      gaps_genes = geneGaps,
      gaps_genes.group_labels = geneGaps.labels,
      geneList = unlist(facsMarkers),
      colType = "gene",
      cellGroups = list(sharedCells)
#      cellGroups = dat$specs$cellGroups
   )
      }

      return(1)
   }


# Cluster size vs true cell type abundance
   if (1) {

# Read in cluster label
      clusterLabels <- read.table(dat$specs$konstantinides18$scClusterLabelsFn,
                        sep="\t",header=FALSE)
      colnames(clusterLabels) <- c("clusterID","clusterName")

# Read in cluster size

      clusterAss <- read.csv(dat$specs$konstantinides18$scClusterAssFn,
                             header=TRUE)
      colnames(clusterAss) <- c("cellID", "clusterID")
      clusterAss$clusterID <- gsub("[^0-9]","",clusterAss$clusterID,perl=TRUE)
      clusterSize <- table(clusterAss$clusterID)
      print(clusterSize)

      trueAbund.perBrain <- c(
#         "perineurial glia" = ,
#         "subperineurial glia" = ,
#         "ensheathing glia" = ,
#         "neuropile glia" = ,
#         "astrocyte-like glia" = ,
#         "chiasm glia" = ,
#         "chortex glia" = ,
         "T1" = 1600,
         "T4/T5" = 12800,
         "Dm8/Tm5c" = 1200,
         "Mt1" = NA,
         "Tm5ab" = NA,
         "T2/T3" = 3200,
         "Lawf1/Lawf2" = 400,
         "C2/C3" = 3200,
         "Dm12" = 250,
         "Tm9" = 1600,
         "TmY14" = NA,
         "Dm2" = NA,
         "Tm1/TmY8" = NA,
         "Mi1/Tm2/Tm3" = 4800,
         "Pm3" = 50,
         "Pm1/Pm2" = NA
      )
      trueAbund.perBrain <- na.omit(trueAbund.perBrain)
      print("Number of celltype/clusters:")
      print(length(trueAbund.perBrain))


      plotMat <- matrix(data = rep(NA, 2 * length(trueAbund.perBrain)), ncol=2)
      rownames(plotMat) <- names(trueAbund.perBrain)
      colnames(plotMat) <- c("clusterSize", "trueAbund.perBrain")

      for (celltype in names(trueAbund.perBrain)) {
         curClusterID <- clusterLabels$clusterID[clusterLabels$clusterName == celltype]
         print(paste0(celltype," is cluster ID ", curClusterID))
         plotMat[celltype,"clusterSize"] <- clusterSize[names(clusterSize) == curClusterID]
         plotMat[celltype,"trueAbund.perBrain"] <- trueAbund.perBrain[celltype]
      }
      print(plotMat)

      if (0) {
      pdf(paste0(dat$specs$outDir,"/",figName,"_scClusterSize_vs_trueAbund.pdf"))
      plot( plotMat[,"trueAbund.perBrain"],
            plotMat[,"clusterSize"],
            xlab = "true abundance (cells per brain)",
            ylab = "single cell cluster size",
            main="",
            pch=20,
            type="n")

      text( plotMat[,"trueAbund.perBrain"],
            plotMat[,"clusterSize"],
            rownames(plotMat),
            cex=0.8)

      dev.off()
      }


      obs.relabund <- plotMat[,"clusterSize"] / plotMat["T1","clusterSize"]
      exp.relabund <- plotMat[,"trueAbund.perBrain"] / plotMat["T1","trueAbund.perBrain"]
      obs.v.exp.relabund <- obs.relabund / exp.relabund
      obs.v.exp.relabund <- sort(obs.v.exp.relabund)

      pdf(paste0(dat$specs$outDir,"/",figName,"_relAbundBarPlot.pdf"),
          width=5,height=3)
      par(mar=c(4.5,5,1,1))
#      barplot(rbind(obs.relabund,exp.relabund),
      plot(1:length(obs.v.exp.relabund),
           obs.v.exp.relabund,
           xlab = "",
           ylab = "Observed / Expected\nabundance",
           xaxt="n",
           main="",
           las=1,
           pch=20,
           cex=1.5,
           bty="n",
           log="y")
      abline(h=1,lwd=2,col="darkgray")
      barnames <- names(obs.v.exp.relabund)
      barnames <- gsub("/","/\n",barnames)
      mtext(barnames, 1, at=1:length(obs.v.exp.relabund),
            line=1.5,
            cex=0.8)
      mtext("Labeled single cell clusters", 1,
            line=3, cex=1)
      text(7,0.8,"Normalized to T1 cluster",cex=1,col="darkgray",adj=c(0,1))
      dev.off()
      print(obs.v.exp.relabund)

   }

}


loadExpr.Konstantinides18 <- function(dat) {

# Purpose: Load locally re-processed kallisto results for Konst's FACS-seq
   
   konstantinides18 <- loadExpr(dat$specs,
                                sampleInfo = dat$specs$konstantinides18$sampleInfoFn ) 
   return(konstantinides18)

}


analyzeRepPonConcordance <- function(dat, figName = "msFigS4A") {
# Purpose: analyze replicate concordance in P calls;

# simple def: on in one (>= 0.8), off in other ( <= 0.2)
# denominator: number of genes with <= 0.2 or >= 0.8 calls.

   modSamples <- gsub("tpm.","",colnames(dat$expr$rawpon$samples$p))
   modDrivers <- sort(colnames(dat$expr$rawpon$drivers$raw))

   concords <- c()
   for (driver in modDrivers) {

      curSamples <- intersect(modSamples,
                     dat$specs$sampleInfo$sample_name[
                        dat$specs$sampleInfo$driver %in% c(driver)])
      nSamples <- length(curSamples)

      for (i in 1:(nSamples - 1)) {
         sample.i <- paste0("tpm.",curSamples[i])
         pon.i <- exp(dat$expr$rawpon$samples$p[,
                     paste0("tpm.",curSamples[i])])

         for (j in (i + 1):nSamples) {
            sample.j <- paste0("tpm.",curSamples[j])
            pon.j <- exp(dat$expr$rawpon$samples$p[,
                        paste0("tpm.",curSamples[j])])

            samplePair <- paste0(sample.i,"-",sample.j)

            curConcord.numer <- sum((pon.i <= 0.2 & pon.j <= 0.2) |
                                    (pon.i >= 0.8 & pon.j >= 0.8))
            curConcord.denom <- sum((pon.i <= 0.2 | pon.i >= 0.8) &
                                    (pon.j <= 0.2 | pon.j >= 0.8))


            concords[samplePair] <- sum((pon.i < 0.5 & pon.j < 0.5) |
                                        (pon.i > 0.5 & pon.j > 0.5)) /
                                    length(pon.i)
            concords[samplePair] <- c(curConcord.numer / curConcord.denom)
         }
      }
   }

   pdf(paste0(dat$specs$outDir,"/",figName,"_replicate_concordance.pdf"),
          height=2,width=2.5)
   par(mar=c(3,3.5,0.5,1),ps=10)
   raw.hist <- hist(concords,plot=FALSE)
   plot(raw.hist,
         main="",xlab="",ylab="",xaxt="n",yaxt="n",
         lwd=2,
         las = 1)
   axis(1, cex.axis=0.8, mgp=c(3,0.5,0), las=1)
   axis(2, cex.axis=0.8, las=1)
   mtext("Concordance of inferred state", 1, line=1.5)
   mtext("Replicate pairs", 2, line=1.8)
        
   dev.off()
   return(concords)

}
