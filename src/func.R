require(qdapTools)
require(stringr)

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

ReadFilesPerMirna <- function(out.dir, cur.dir){
  setwd(out.dir)
  temp <- list.files(pattern="*.bed")
  tables <- lapply(temp, fread)
  names(tables) <-  temp
  rm(temp)
  setwd(cur.dir)
  return(tables)
}

ConvertToDt <- function(list){
  names(list) <- tstrsplit(names(list),"_",fixed=TRUE)[[1]]
  return(as.data.table(list_df2df(list, col1 = 'mirna')))
}

PrepareTriplexatorOutput <- function(dt){
  dt <- SplitRegionColumn(dt = dt[,1:2],
                          region.col = colnames(dt)[1])
  dt[, `Sequence-ID`:=RenameLongMirnaName(column = `Sequence-ID`)]
  setnames(dt, colnames(dt),
           c('mirna', 'enh.chr', 'enh.start', 'enh.end'))
  dt$enh.start <- as.numeric(dt$enh.start)
  dt$enh.end <- as.numeric(dt$enh.end)
  return(dt)
}
SplitRegionColumn <- function(dt, region.col){
  dt[, c("chr", 'coord') := tstrsplit(get(region.col), ":", fixed=TRUE)]
  dt[, c("start", 'end') := tstrsplit(`coord`, "-", fixed=TRUE)]
  dt <- unique(dt[, -c('coord', region.col), with=F])
  return(dt)
}

RenameLongMirnaName <- function(column){
  # converts 'mmu-miR-223-3p' style into 'Mir223'
  res <- str_replace(column, 'mmu-miR-', '')
  res <- paste0('Mir', sub('-.*', '', res))
  return(res)
}

PrepareTriplexatorOutput <- function(dt){
  dt <- SplitRegionColumn(dt = dt[,1:2],
                          region.col = colnames(dt)[1])
  dt[, `Sequence-ID`:=RenameLongMirnaName(column = `Sequence-ID`)]
  setnames(dt, colnames(dt),
           c('mirna', 'enh.chr', 'enh.start', 'enh.end'))
  dt$enh.start <- as.numeric(dt$enh.start)
  dt$enh.end <- as.numeric(dt$enh.end)
  return(dt)
}
SplitRegionColumn <- function(dt, region.col){
  dt[, c("chr", 'coord') := tstrsplit(get(region.col), ":", fixed=TRUE)]
  dt[, c("start", 'end') := tstrsplit(`coord`, "-", fixed=TRUE)]
  dt <- unique(dt[, -c('coord', region.col), with=F])
  return(dt)
}

RenameLongMirnaName <- function(column){
  # converts 'mmu-miR-223-3p' style into 'Mir223'
  res <- str_replace(column, 'mmu-miR-', '')
  res <- paste0('Mir', sub('-.*', '', res))
  return(res)
}

PrepareTriplexatorOutput <- function(dt){
  dt <- SplitRegionColumn(dt = dt[,1:2],
                          region.col = colnames(dt)[1])
  dt[, V2:=RenameLongMirnaName(column = V2)]
  setnames(dt, colnames(dt),
           c('mirna', 'enh.chr', 'enh.start', 'enh.end'))
  dt$enh.start <- as.numeric(dt$enh.start)
  dt$enh.end <- as.numeric(dt$enh.end)
  return(dt)
}
