library(data.table)
library(R.utils)
source('src/func.R')


alpha = 0.05

#parse args from the main script
method.enh <- commandArgs(T)[1]
out.dir <- commandArgs(T)[2]
gene.expr.path <- commandArgs(T)[3]
mirnas.expr.path <- commandArgs(T)[4]
enh.gene.assoc.sign.path <- commandArgs(T)[5]

if (!isAbsolutePath(out.dir)){
  out.dir <- file.path(getwd(), out.dir)
}

gene.expr <- fread(gene.expr.path)
mirnas.expr <- fread(mirnas.expr.path)
enh.gene.assoc.sign <- fread(enh.gene.assoc.sign.path)

if(method.enh == 'seed_match_needle'){
  enh.filtered <- ReadFilesPerMirna(dir = out.dir)
  
  enh.filtered <- ConvertToDt(enh.filtered)
  setnames(enh.filtered,
           c('V1', 'V2', 'V3'),
           c('enh.chr', 'enh.start', 'enh.end'))
}else if(method.enh == 'miranda'){
  enh.miranda.path <- file.path(out.dir, 'miranda.enh.tsv')
  enh.filtered <- fread(enh.miranda.path)
  enh.filtered <- unique(enh.filtered[,-c('max.energy', 'max.score')])
}else if(method.enh == 'triplexator'){
  enh.triplexator.path <- file.path(out.dir, 'triplex_search.summary')
  enh.filtered <- fread(enh.triplexator.path)[,-c('V11')]
  enh.filtered <- PrepareTriplexatorOutput(dt = enh.filtered)
}

enh.found <- unique(enh.gene.assoc.sign[match(enh.filtered$enh.start, enh.gene.assoc.sign$enh.start)][, c('Enhancer')])
enh.found.genes <- unique(enh.gene.assoc.sign[Enhancer %in% enh.found$Enhancer])
enh.filtered <- merge(enh.found.genes,
                      enh.filtered,
                      by = c('enh.chr', 'enh.start', 'enh.end'),
                      allow.cartesian=T
)[,c('mirna', 'enh.chr', 'enh.start', 'enh.end', 'Gene.Name', 'Enhancer')]


# STATE IN README that gene expression table should be already aggregated if you want to sum expression of transcripts per gene and 
# must include same nomenclature of gene symbols in Symbol column as in Gene.Name column in file of enh-gene interaction
associated.genes.expr <- enh.filtered[, gene.expr[Symbol %in% .SD$Gene.Name], by=mirna]

mirnas.expr.list <- sapply(mirnas.expr$Symbol, function(x){mirnas.expr[Symbol == x]}, USE.NAMES = T,simplify = F)

corrs.all <- CalcCorrsAllMirnasVsAllGenes(associated.genes.expr, mirnas.expr.list)
trio.dt <- GetSignifCorrsAllMirnasVsAllGenes(enh.filtered,
                                             corrs.all,
                                            alpha
                                            )
trio.dt.prepared <- PrepareAllCorrsTable(trio.dt, get.max=T)
write.table(trio.dt.prepared,
            file.path(out.dir, 'mir_enh_gene_trios.tsv'),
            sep = '\t', quote = F, row.names = F)


if(method.enh == 'seed_match_needle'){
    out.dir.enh.active <- file.path(out.dir, 'enh_active')
    dir.create(out.dir.enh.active, showWarnings = FALSE)

    enh.active.bed <- sapply (unique(trio.dt$mirna), function(de.mirna){
      print(de.mirna)
      as.data.table(unique(trio.dt[mirna == de.mirna, c('enh.chr', 'enh.start', 'enh.end')]))
    }, simplify = F, USE.NAMES = T)

    lapply(names(enh.active.bed), function(de.mirna){
      write.table(enh.active.bed[[de.mirna]],
                  file.path(out.dir.enh.active, paste0(de.mirna, '_enh_active.bed')),
                  quote = F,
                  col.names = F,
                  row.names = F,
                  sep = '\t')
    })

}