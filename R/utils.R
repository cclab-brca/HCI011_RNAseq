################################################
# Utils for HCI011 RNA-seq analysis
################################################

library(svglite)
library(org.Hs.eg.db)
library(edgeR)
library(scales)
library(GSEABase)
library(limma)
library(tibble)

###############################################
# wrap for opening a plot (png, ps, pdf or svg)
###############################################
plot_start = function(fn, w, h, device=NULL, res=72, pointsize=10)
{
  curr_ext = gsub(".+\\.", "", fn)
  if (is.null(device)) {
    device = curr_ext
  }
  
  if (device == "png") {
    png(filename=sub(paste0(curr_ext, "$"), "png", fn), width=w, height=h, res=res, pointsize = pointsize)
  }
  else if (device == "ps") {
    postscript(file=sub(paste0(curr_ext, "$"), "ps", fn), width=w/res, height=h/res)
  }
  else if (device == "pdf") {
    pdf(file=sub(paste0(curr_ext, "$"), "pdf", fn), width=w/res, height=h/res)
  }
  else if (device == "svg") {
    svglite::svglite(file=sub(paste0(curr_ext, "$"), "svg", fn), width=w/res, height=h/res, pointsize = pointsize)
  }
  
  else {  
    stop(sprintf("unknown output device type: %s", device))
  }
}


###############
# DE with edgeR
###############
rnaseq_de_by_edgeR <- function(fit, contrast, gene2entrez, odir, base_name, n_top_go = 20, min_fc=1.2, plot_max_fdr=0.05, plot_min_abs_lfc=2, plot_max_genes=100, gen_csvs=T, gen_plots=T, pv_cols = colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100), device="png")
{
  tr = glmTreat(fit, contrast=contrast, lfc=log2(min_fc))
  
  # GO enrichments
  go = goana(tr, geneid = gene2entrez[rownames(tr$genes)])
  go_up = topGO(go, sort="up", number=n_top_go)
  go_down = topGO(go, sort="down", number=n_top_go)
  
  if (gen_csvs) {
    if (nrow(go_up) > 0) {
      write.csv(go_up, file=sprintf("%s/%s_GO_up.csv", odir, base_name), row.names = F)
    }
    
    if (nrow(go_down) > 0) {
      write.csv(go_down, file=sprintf("%s/%s_GO_down.csv", odir, base_name), row.names = F)
    }
  }
  
  # DE genes
  xx = summary(decideTests(tr))
  message(sprintf("%s:\n\t\t\t\t%s", base_name, paste(rownames(xx), xx[,1], sep=": ", collapse="\t")))
  
  de_df = topTags(tr, n=1e6, p.value = 1) # info on all genes
  if (nrow(de_df$table) > 0) {
    if (gen_csvs) {
      write.csv(de_df$table %>% tibble::rownames_to_column(var="ensembl_gene_id"), file=sprintf("%s/%s.csv", odir, base_name), row.names = F)
    }
    
    de_df_f = de_df$table %>% 
      tibble::rownames_to_column(var="ensembl_gene_id") %>% 
      filter(abs(logFC) >= plot_min_abs_lfc & FDR <= plot_max_fdr) %>% 
      arrange(-logFC)
    
    if (nrow(de_df_f) > 0 & gen_plots) {
      bh = 12
      de_df_f$col = pv_cols[rescale(pmin(-log10(de_df_f$FDR), 5), to = c(1, 100))]
      pv_col_breaks = c(plot_max_fdr, 10**seq(-2, -5))
      de_up = de_df_f[de_df_f$logFC > 0, ] %>% arrange(FDR) %>% head(plot_max_genes) %>% arrange(-logFC)
      de_down = de_df_f[de_df_f$logFC < 0, ] %>% arrange(FDR) %>% head(plot_max_genes) %>% arrange(-logFC)
      
      nb = max(nrow(de_up), nrow(de_down))

      plot_start(sprintf("%s/%s_lfc%.2f_FDR%e.%s", odir, base_name, plot_min_abs_lfc, plot_max_fdr, device), 600, (nb + 0.2) * bh + 100, device=device)
      layout(matrix(1:2, nrow=1))
      
      par(mar=c(4,10,2,0))
      if (nrow(de_up) > 0) {
        barplot(de_up$logFC, horiz=T, col = de_up$col, ylim = c(0, nb * 1.2), xlab="Up")
        axis(2, seq(0.7, by=1.2, length.out = nrow(de_up)), labels = de_up$Symbol, tick = F, las=2, cex.lab=0.5)
      } else {
        plot.new()
        plot.window(1:2, 1:2)
        text(1, 1.5, "No signif. up genes")
      }
      legend(x=plot_min_abs_lfc + 0.5, y=nb*1.2, legend=sprintf("%s%.2e", ifelse(pv_col_breaks == min(pv_col_breaks), "<= ", ""), pv_col_breaks), fill = pv_cols[round(seq(1, 100, len=length(pv_col_breaks)))], title="FDR", cex=0.8)
      
      par(mar=c(4,0,2,10))
      if (nrow(de_down) > 0) {
        barplot(de_down$logFC, horiz=T, col = de_down$col, ylim = c(0, nb * 1.2), xlab="Down")
        axis(4, seq(0.7, by=1.2, length.out = nrow(de_down)), labels = de_down$Symbol, tick = F, las=2, cex.lab=0.5)
      } else {
        plot.new()
        plot.window(1:2, 1:2)
        text(1, 1.5, "No signif. down genes")
      }
      title(main=base_name, cex=3, outer=T, line=-2)
      dev.off()      
    }
  }
  invisible(list(de_df=de_df, de_df_f=de_df_f, go_up=go_up, go_down=go_down))
}

#########################################################################
# GSEA with limma's camera function (adapted from Alistair Martin's code)
#########################################################################
rnaseq_gsea <- function(v, design_matrix, contrast, gset, FDR_thresh=0.05, spec="human", igc=0.01, ...)
{
  if (gset == "GSEA_HM" && spec == "human") {
    gmt = GSEABase::getGmt("../data/h.all.v7.1.entrez.gmt")
    gs = sapply(gmt, function(g) g@geneIds)
    names(gs) = names(gmt)
  } else {
    stop("Currently only supporting human GSEA HM")
  }
  ids <- mapIds.ensembl(rownames(v), spec)
  idx <- ids2indices(gs, id=ids) 
  gsea <- as.data.table(camera(v, idx, design_matrix, contrast = contrast, inter.gene.cor=igc,...), keep.rownames="pathway")[FDR < FDR_thresh]
  gsea
}

mapIds.ensembl <- function(x, spec = "human")
{ 
  if(spec == "human") { 
    db <- org.Hs.eg.db 
  } else if(spec == "mouse") { 
    db <- org.Mm.eg.db 
  } else { 
    stop("Species not found")
  }
  mapIds(db, keys=x, keytype = "ENSEMBL",column = "ENTREZID", multiVals = "first")
}
