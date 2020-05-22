## RNA-seq analysis of HCI011 PDX for the paper: 
## FOXM1 is a biomarker of resistance to PI3Kα inhibition in ER+ breast cancer that is detectable 
## using metabolic imaging

source("utils.R")

library(data.table)
library(biomaRt)

library(plyr)
library(dplyr)

library(RColorBrewer)
library(pheatmap)

# Init parameters
dev = "png"
nm = "HCI011"
subm_cols = c(HCI011="#377EB8", HCI011R="#E41A1C")
treat_cols = c(gdc='purple', vehicle='darkgrey')

gns = brewer.pal(n=9, 'YlGn')
fdr_breaks = c(0, 1e-3, 1e-2, 5e-2, 1)
fdr_cols = c("<= 0.001" = gns[9], "<= 0.01" = gns[6], "<= 0.05" = gns[3], "N.S." = 'white')

figs_dir = "../output"
data_dir = "../data"

dir.create(figs_dir, recursive = T, showWarnings = F)

# Read samples metadata and counts
exp_md = read.csv(paste0(data_dir, "/20200304_Metadata_for_HCI011_seq.csv"), header=T) %>%
  mutate(Submodel = factor(gsub("-", "", Submodel)))

md = read.table(paste0(data_dir, "/HCI011_featureCounts.meta.tsv"), header=T, sep="\t", stringsAsFactors = F)
r_counts = fread(paste0(data_dir, "/HCI011_featureCounts.counts.tsv.gz"), check.names = F, data.table = F)
rownames(r_counts) = r_counts$Geneid

hg_ind = grepl("^ENSG", r_counts$Geneid)

# Build lists of gene counts splited to human and mouse per model. Note that the human genome build is Ensembel 75.
ens_url = biomaRt::listEnsemblArchives() %>%
  dplyr::filter(version == 75) %>% 
  dplyr::select(url)

ens = biomaRt::useMart(host=ens_url$url, biomart='ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')
ens2nm = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene'), filters = "ensembl_gene_id", values = as.character(r_counts[hg_ind, 'Geneid']), mart = ens) %>%
  group_by(ensembl_gene_id) %>%
  group_by(ensembl_gene_id) %>% 
  summarise(hgnc_symbol=paste0(unique(hgnc_symbol), collapse="-"), entrezgene=paste0(unique(entrezgene), collapse="-")) %>%
  rename(Geneid = "ensembl_gene_id")

ens2gene = ifelse(ens2nm$hgnc_symbol == "", ens2nm$Geneid, ens2nm$hgnc_symbol)
names(ens2gene) = ens2nm$Geneid

ens2entrez = ifelse(ens2nm$entrezgene == "", ens2nm$Geneid, ens2nm$entrezgene)
names(ens2entrez) = ens2nm$Geneid

gene2ens = names(ens2gene)
names(gene2ens) = ens2gene

hg_counts = r_counts[hg_ind, ]
hg_counts = merge(hg_counts, ens2nm) 
rownames(hg_counts) = as.character(hg_counts$Geneid) 
  
mm_counts = r_counts[!hg_ind, ]
  
c_md = merge(exp_md, md[, c('Sample.name', 'fname')], sort=F)
counts = list(meta = exp_md,
              human = list(data = hg_counts[, c_md$fname],
                           gene.meta = hg_counts[, c(1:6, ncol(hg_counts))]),
              mouse = list(data = mm_counts[, c_md$fname],
                           gene.meta = mm_counts[, 1:6]))
colnames(counts$human$data) = colnames(counts$mouse$data) = c_md$Sample.name

# Prepare edgeR object
targets = counts$meta
group = factor(paste0(targets$Submodel, ".", targets$Treatment))
  
y = DGEList(counts$human$data, group=group, genes = data.frame(Symbol=ens2gene[counts$human$gene.meta$Geneid]))

# Filter weakly expressed genes
keep = filterByExpr(y)
y = y[keep, , keep.lib.sizes=FALSE]
y = calcNormFactors(y)

# DE analysis
message("Model + treamtnet DE analysis...")
group = y$samples$group

design = model.matrix(~ 0 + group)
colnames(design) = levels(group)

y = estimateDisp(y, design, robust=TRUE)
fit = glmQLFit(y, design, robust=TRUE)

# sample and model + treatment mean
message("Calculating normalised sample (and condition) mean values...")
y_n = cpm(y, log = F)
y_ln = cpm(y, log = T)
y_ln_enr = y_ln - log2(rowMeans(2 ** y_ln))

y_n_geom = t(apply(y_n, 1, function(v) { tapply(X=v, INDEX=y$samples$group, FUN=function(y) { exp(mean(log(1 + y))) - 1 }) }))


reg_c = 2
y_n_geom_enr = log2((y_n_geom + reg_c) / rowMeans(y_n_geom + reg_c))
  
# Checking DE by model: resistent vs sensitive.
message("HCI011 vs HCI011R DE analysis...")
plot_min_abs_lfc = 1
plot_max_fdr = 0.05
plot_max_genes = 50

m_group = factor(targets$Submodel)
m_design = model.matrix(~ 0 + m_group)
colnames(m_design) = levels(m_group)
  
y_m = estimateDisp(y, m_design, robust=TRUE)
fit_m = glmQLFit(y_m, m_design, robust=TRUE)
  
contrast = makeContrasts(contrasts = sprintf("%sR - %s", nm, nm), levels=m_design)
  
models_de = rnaseq_de_by_edgeR(fit_m, contrast, ens2entrez, figs_dir, sprintf("%sR_vs_%s", nm, nm), plot_max_fdr = plot_max_fdr, plot_min_abs_lfc = plot_min_abs_lfc, plot_max_genes= plot_max_genes, device=dev)

## PI3K pathway genes
# 
# Focusing on 43 PI3K resistant genes from "Systematic Functional Characterization
# of Resistance to PI3K Inhibition in Breast Cancer", Le et al. 2016.
message("PI3K pathway genes...")

foc_genes = c(read.csv(sprintf("%s/PI3K_res_genes_from_Le2016.csv", data_dir), header=T, stringsAsFactors = F)$Gene, 'PIK3CA', 'PIK3CB')
missed_foc = setdiff(foc_genes, ens2gene[rownames(y_n_geom_enr)])
if (length(missed_foc) > 0) {
  message(sprintf("%s: Missing %d genes (%s)", nm, length(missed_foc), paste0(missed_foc, collapse=", ")))
}
  
g = levels(y$samples$group)

cond_ann = data.frame(row.names=g, Submodel=gsub("\\..*", "", g), Treatment=gsub(".*\\.", "", g))

ann_cols = list(Treatment=treat_cols, Submodel=subm_cols)

cw = 12
zlim = 0.6

y_disp = t(y_n_geom_enr[gene2ens[setdiff(foc_genes, missed_foc)], ])

colnames(y_disp) = setdiff(foc_genes, missed_foc)
hc = hclust(dist(cor(y_disp)), method='ward.D2')

y_disp = y_disp[paste0(nm, c('R.gdc', 'R.vehicle', '.gdc', '.vehicle')), hc$order]

gene_fdr = data.frame(row.names=colnames(y_disp),  FDR=as.character(cut(pmax(models_de$de_df$table[gene2ens[colnames(y_disp)], 'FDR'], min(fdr_breaks)), fdr_breaks, include.lowest = T, labels=names(fdr_cols))))

plot_start(sprintf("%s/%s_Fig6G.png", figs_dir, nm), ncol(y_disp) * cw + 300, nrow(y_disp) * cw + 300, dev)
pheatmap(pmin(pmax(y_disp, -zlim), zlim), breaks=seq(-zlim, zlim, len=101), cluster_rows=F, cluster_cols=F, cellwidth=cw, cellheight = cw, annotation_row=cond_ann, annotation_col = gene_fdr, annotation_colors = c(ann_cols, list(FDR=fdr_cols)), main=nm)
dev.off()


# FOXM1 Target genes on GDC resistant vs GDC vehicle
message("FOXM1 target genes...")
zlim = 0.6

foxm1_tgts = read.csv(sprintf("%s/FOXM1_targets.csv", data_dir), header=T, stringsAsFactors = F)$Gene
missed_foxm1 = setdiff(foxm1_tgts, ens2gene)

if (length(missed_foxm1) > 0) {
  message(sprintf("FOXM1 target genes: %d not found (%s)", length(missed_foxm1), paste0(missed_foxm1, collapse = ", ")))
  foxm1_tgts = setdiff(foxm1_tgts, missed_foxm1)
}

missed_foxm1 = setdiff(foxm1_tgts, ens2gene[rownames(y_n_geom_enr)])
if (length(missed_foxm1) > 0) {
  message(sprintf("FOXM1 target genes: %d filtered because of low counts (%s)", length(missed_foxm1), paste0(missed_foxm1, collapse = ", ")))
  foxm1_tgts = setdiff(foxm1_tgts, missed_foxm1)
}

contrast = makeContrasts(contrasts = sprintf("%sR.gdc - %s.gdc", nm, nm), levels=design)
 
c_gsea = rnaseq_gsea(y, design, contrast, "GSEA_HM", FDR_thresh=0.05, spec="human")
if (nrow(c_gsea) > 0) {
  write.csv(c_gsea, sprintf("%s/%sR.gdc_vs_%s.gdc_GSEA_HM.csv", figs_dir, nm, nm), row.names = F)
}

c_de = rnaseq_de_by_edgeR(fit, contrast, ens2entrez, figs_dir, sprintf("%sR.gdc_vs_%s.gdc", nm, nm), gen_csvs = T, gen_plots=F)

y_disp = t(y_n_geom_enr[gene2ens[foxm1_tgts], ])
colnames(y_disp) = foxm1_tgts
  
y_disp = y_disp[paste0(nm, c('R.gdc', '.gdc')), ]
y_disp = y_disp[, order(apply(y_disp, 2, diff))]
cond_ann = data.frame(row.names=rownames(y_disp), Submodel=gsub("\\..*", "", rownames(y_disp)), Treatment=gsub(".*\\.", "", rownames(y_disp)))

gene_fdr = c_de$de_df$table[gene2ens[colnames(y_disp)], 'FDR']
gene_ann = data.frame(row.names=colnames(y_disp),  FDR=as.character(cut(pmax(gene_fdr, min(fdr_breaks)), fdr_breaks, include.lowest = T, labels=names(fdr_cols))))
fdr_filt = 0.05
c_y_disp = y_disp[, gene_fdr <= fdr_filt]
plot_start(sprintf("%s/%s_Fig6F_%s.%s", figs_dir, nm, ifelse(fdr_filt == 1, "all", "filt"), dev), ncol(c_y_disp) * cw + 300, nrow(c_y_disp) * cw + 300, dev)
pheatmap(pmin(pmax(c_y_disp, -zlim), zlim), breaks=seq(-zlim, zlim, len=101), cluster_rows=F, cluster_cols=F, cellwidth=cw, cellheight = cw, annotation_row=cond_ann, annotation_col = gene_ann, annotation_colors = c(ann_cols, list(FDR=fdr_cols)), main=nm)
dev.off()

# Sup Fig 5A: Focus on FOXM1, HK-II, LDHA, PTEN 
gs = c('FOXM1', 'HK2', 'LDHA', 'PTEN')

pairs = combn(rev(levels(y$samples$group)), 2)

de_res = NULL

for (i in 1:ncol(pairs)) {
  contrast = makeContrasts(contrasts = sprintf("%s - %s", pairs[1, i], pairs[2, i]), levels=design)
  
  c_de = rnaseq_de_by_edgeR(fit, contrast, ens2entrez, figs_dir, paste0(pairs[, i], collapse="_vs_"), gen_csvs = F, gen_plots = F)
  
  c_res = c_de$de_df$table %>% filter(Symbol %in% gs) %>% mutate(v1=pairs[1, i], v2=pairs[2, i]) 
  de_res = rbind(de_res, c_res)
}

write.csv(de_res, file = sprintf("%s/%s_SFig5A_de.csv", figs_dir, nm), row.names = F, quote=F)
  
glob_ylim = range(y_ln_enr[gene2ens[gs], ])
for (g in gs) {
  plot_start(sprintf("%s/%s_SFig5A_%s.%s", figs_dir, nm, g, dev), 400, 500, dev)
  par(mar=c(8,4,4,2))
  v = y_ln_enr[gene2ens[g], ]
  f = y$samples$group
  mv = tapply(v, f, mean)
  stripchart(v ~ f, las=2, main=paste(nm, g), method='jitter', pch=19, jitter=0.1, cex=0.7, vertical=T, ylab='Expr enr (log2)', ylim=glob_ylim)
  grid(nx=NA, ny=NULL, col="darkgray", lty=3, lwd=1)
  abline(h=0, col="darkgray", lty=3, lwd=2)
  segments(x0=seq(0.75, by=1, length=length(mv)), x1=seq(1.25, by=1, len=length(mv)), y0=mv, y1=mv, col='red', lwd=2)
  
  dev.off()
}

