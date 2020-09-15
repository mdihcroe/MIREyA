library(data.table)
source('src/func.R')


alpha = 0.05

#parse args from the main script
method.enh <- commandArgs(T)[1]
out.dir <- commandArgs(T)[2]
gene.expr.path <- commandArgs(T)[3]
mirnas.expr.path <- commandArgs(T)[4]
enh.gene.assoc.sign.path <- commandArgs(T)[5]

if(!file.exists(gene.expr.path)){
  print('gene.expr.path does not exist. Please provide valid full path to the file with -ge or --gene_expression argument')
  stop("Not enough observations in 'x': n < 500")
}

# method.enh = 'miranda'
# out.dir <- '/home/anna/Projects/MIREyA/out/miranda_out'
# gene.expr.path <- '../MIREyA/data/DE_gene_expression.tsv'
# mirnas.expr.path <- '../MIREyA/data/DE_mirnas_expression.tsv'
enh.gene.assoc.sign <- fread(enh.gene.assoc.sign.path)

if(method.enh == 'seeds'){
  enh.filtered <- ReadFilesPerMirna(dir = file.path(out.dir, 'enh_with_seeds'))
  
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

# gene.expr.all = fread('/home/anna/Projects/mirna.targets/out/gene_expr_FANTOM.txt')
# load(paste0('data/', type, 'degs.unique.RData'))
# gene.expr.all <- GetDEGenesExpr(gene.expr.all,
#                                 degs)
# write.table(gene.expr.all[,-c(1:4, 6)],
#             '../MIREyA/data/DE_gene_expression.tsv',
#             sep = '\t', quote = F, row.names = F)
gene.expr <- fread(gene.expr.path)


# STATE IN README that gene expression table should be already aggregated if you want to sum expression of transcripts per gene and 
# must include same nomenclature of gene symbols in Symbol column as in Gene.Name column in file of enh-gene interaction
associated.genes.expr <- enh.filtered[, gene.expr[Symbol %in% .SD$Gene.Name], by=mirna]

# extract expression of miRNAs of interest
# since names differ everywhere let the user supply it as argument
# la  =as.data.table(list_df2df(de.mirna.expr, col1='mirna'))
# la= la[,-c('mirna', 'chr', 'strand')]
# write.table(la,
#             '../MIREyA/data/DE_mirnas_expression.tsv',
#             sep = '\t', quote = F, row.names = F)
mirnas.expr <- fread(mirnas.expr.path)
mirnas.expr.list <- sapply(mirnas.expr$Symbol, function(x){mirnas.expr[Symbol == x]}, USE.NAMES = T,simplify = F)

corrs.all <- CalcCorrsAllMirnasVsAllGenes(associated.genes.expr, mirnas.expr.list)
trio.dt <- GetSignifCorrsAllMirnasVsAllGenes(enh.filtered,
                                             corrs.all,
                                            alpha
                                            )
trio.dt <- PrepareAllCorrsTable(trio.dt, get.max=T)
write.table(trio.dt,
            file.path(out.dir, 'mir_enh_gene_trios.tsv'),
            sep = '\t', quote = F, row.names = F)