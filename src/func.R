CalcCorrWithMirna <- function(expr.per.gene, mirnas.expr.list){
  mirna <- expr.per.gene$mirna
  cor.res <- cor.test(as.numeric(mirnas.expr.list[[mirna]][, -c('Symbol')]),
                      as.numeric(expr.per.gene[, -c('mirna', 'Symbol')]),
                      method = 'spearman')
  res <- data.table(mirna = mirna,
                    symbol = expr.per.gene$Symbol,
                    `corr(miRNA, gene)` = cor.res$estimate,
                    p.value = cor.res$p.value)
  return(res)
}

CalcCorrsAllMirnasVsAllGenes <- function(associated.genes.expr, mirnas.expr.list){
  corrs <- associated.genes.expr[, CalcCorrWithMirna(.SD, mirnas.expr.list), by = 1:nrow(associated.genes.expr)]
  corrs <- corrs[complete.cases(corrs)]
  corrs[, ('p.value adj') := p.adjust(.SD$p.value, method = 'BH'), by=mirna]
}

GetSignifCorrsAllMirnasVsAllGenes <- function(enh.filtered, corrs.all, alpha){
  corrs <- corrs.all[`p.value adj`<alpha & `corr(miRNA, gene)`>0]
  enh.filtered <- merge(enh.filtered,
                        corrs,
                        by.x = c('mirna', 'Gene.Name'),
                        by.y = c('mirna', 'symbol'),
                        all.y = T)
  enh.filtered[, nrow:=NULL]
}

PrepareAllCorrsTable <- function(dt, get.max, rna.column = 'mirna', corr.coef.column = "corr(miRNA, gene)"){
  dt[,enhancer:=paste0(enh.chr, ':', enh.start, '-', enh.end)]
  if (get.max == T){
    dt.max <- LeaveMaxCorrs(dt, by.cols=c(rna.column, "enhancer", "Gene.Name"), 
                            cols = c(corr.coef.column))
    # dt.max[,enhancer.mm10:=dt[match(dt.max$enhancer, dt$enhancer)]$enhancer.mm10]
    dt = unique(dt.max)
  }
  return(dt)
}

LeaveMaxCorrs <- function(dt, by.cols, cols){
  dt.max <- 
    dt[,  lapply(.SD, max), 
       by=by.cols,
       .SDcols=cols ]
  cols.pval <- colnames(dt)[grep('p.val', colnames(dt))]
  cols.to.leave <- c(colnames(dt.max), cols.pval)
  dt.max <- merge(dt,
                  dt.max,
                  by = colnames(dt.max))[, cols.to.leave, with=F]
  cols.2 <- cols[grep('corr', cols)]
  dt.max[,(cols.2) := round(.SD,2), .SDcols=cols.2]
  dt.max[,(cols.pval) := round(.SD,3), .SDcols=cols.pval]
  dt.max[,`p.value`:=NULL]
  dt.max
}
