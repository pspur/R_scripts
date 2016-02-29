### This script performs Mutation Set Enrichment Analysis

setwd('/home/paulspur/R_scripts/MSEA')
#setwd('/Users/Paul/Documents/MSEA')
library(doParallel)
library(plyr)
library(ggplot2)
library(grid)
library(dplyr)

source('./bin/MSEA.R')

cfg = list(project.stem = 'MSEA_SIM_WGS',
           infile.refgene = './input/ref.gene.final.txt',
           infile.promoter = './input/prom.plot_ready.txt',
           infile.variant = './input/deleterious_nonsilent_SNVs_Indels.plot_ready.txt',
           infile.domain = './input/protein.gbk.regions.symbol.txt',
           infile.tfbs = './input/tfbs.plot_ready.txt',
           output.directory = './output'
)
DATE <- Sys.Date()

cl <- makeCluster(2, outfile='')
registerDoParallel(cl)

match.tfbs.to.promoter <- function(x) {
    chr <- as.character(x[1])
    mstart <- as.numeric(x[2])
    mend <- as.numeric(x[3])
    mstrand <- as.character(x[5])
    chr.match <- prom.to.test[which(chr==prom.to.test$chrom), ]
    chr.strand.match <- chr.match[which(mstrand==chr.match$strand), ]
    chr.strand.match <- chr.strand.match[match(unique(chr.strand.match$Symbol), chr.strand.match$Symbol), ]
    chr.strand.match$pstart <- chr.strand.match$start-chr.strand.match$mut_pos 
    cs.match <- chr.strand.match[which(chr.strand.match$pstart<=mstart & mend<=chr.strand.match$pstart+1000), ]
    if (length(cs.match[,1])==0) {
        symbol = 'n'
        domain.start <- -1
        domain.end <- -2
    } else {
        symbol = cs.match$Symbol
        domain.start <- mstart-cs.match$pstart
        domain.end <- domain.start+(mend-mstart)
    } 
    return(list('symbol'=symbol,'domain.start'=domain.start,'domain.end'=domain.end))
}

var.freq <- function(x) {
    var.present <- apply(x, c(1,2), function(u) {
        gt <- strsplit(u,split=':')[[1]][1];
        ifelse(gt=='0/0'|gt=='./.', u <- 0, u <- 1) 
    })
    freq <- apply(var.present,1,sum)
    return(freq)
}

get.mutloc.and.strand <- function(x) {
    chr <- as.character(x[1])
    mloc <- as.numeric(x[2])
    gene <- as.character(x[7])
    gene.match <- ref.gene[which(gene==ref.gene$V13), ]
    gene.chr.match <- gene.match[which(chr==gene.match$V3), ]
    gene.chr.match$pstart <- ifelse(gene.chr.match$V4=='+', 
                                    gene.chr.match$V5-1000, 
                                    gene.chr.match$V6)
    gene.chr.match$pend <- ifelse(gene.chr.match$V4=='+',
                                  gene.chr.match$V5,
                                  gene.chr.match$V6+1000)
    match <- filter(gene.chr.match, pstart <=mloc, mloc <= pend)
    match <- filter(match, row_number() == 1)
    if (length(row(match))==0) {
        mloc <- -1
        strand <- '0'
    } else {
        mloc <- abs(mloc - match$pstart)
        strand <- match$V4
    }
    output <- list('mloc'=mloc,'strand'=strand)
    return(output)
}

get.ref.gene.length <- function(x) {
    cds.start <- as.numeric(x[7])
    cds.end <- as.numeric(x[8])
    exon.starts <- as.numeric(unlist(strsplit(x[10], ',')))
    exon.ends <- as.numeric(unlist(strsplit(x[11], ',')))
    
    # get the index of first coding exon
    for (i in 1:length(exon.starts)) {
        if (cds.start >= exon.starts[i] & cds.start <= exon.ends[i]) {
            exon.starts[i] <- cds.start
            start.exon <- i
            break
        }
    }
    
    # get the index of the last coding exon
    for (i in 1:length(exon.ends)){
        if (cds.end >= exon.starts[i] & cds.end <= exon.ends[i]) {
            exon.ends[i] <- cds.end
            end.exon <- i
            break
        }
    }
    
    exon.starts <- exon.starts[start.exon:end.exon]
    exon.ends <- exon.ends[start.exon:end.exon]
    
    len <- sum(exon.ends - exon.starts)
    return(len/3) # convert to a.a. length
}

####################################################################################################
# get a.a lengths for RefSeq gene
####################################################################################################
ref.gene <- read.table(cfg$infile.refgene, sep = '\t', header=T, stringsAsFactors=F)

refseq.length <- apply(ref.gene, 1, get.ref.gene.length)
names(refseq.length) <- ref.gene$V2

####################################################################################################
# read promoter variant data
####################################################################################################
prom.pass <- read.table(cfg$infile.promoter, sep ='\t', header=T, stringsAsFactors=F)

nMut.per.prom <- array(rowsum(prom.pass$freq, prom.pass$Symbol, reorder=F))
names(nMut.per.prom) <- unique(prom.pass$Symbol)

# decide on a proper number here
prom_to_test <- names(which(nMut.per.prom >= 4))
prom.to.test <- filter(prom.pass, Symbol %in% prom_to_test)

# print some numbers
print(paste('# input promoters: ', length(nMut.per.prom), 
            '; eligible promoters to test: ', length(prom_to_test), sep=''))

prom.num <- length(unique(prom.to.test$Symbol))
refseq.len.dummy <- rep(1000, prom.num)
names(refseq.len.dummy) <- unique(prom.to.test$Symbol)

####################################################################################################
# read gene variant data
####################################################################################################
mutations <- read.table(cfg$infile.variant, header=T, stringsAsFactors=F, sep='\t')

nMut.per.gene <- tapply(mutations$amino_acid_change, mutations$RefSeq.ID, length) 

# only include genes with 4+ mutations
genes_to_test <- names(which(nMut.per.gene >= 4))
genes_to_test <- genes_to_test[genes_to_test %in% names(refseq.length)]

idx <- match(mutations$RefSeq.ID, genes_to_test)
mutations.to.test <- mutations[!is.na(idx), ]

# print some numbers
print(paste('# input genes: ', length(nMut.per.gene), 
            '; eligible genes to test: ', length(genes_to_test), sep=''))

####################################################################################################
# read domain and transcription factor binding site data (used as domains for promoters)
####################################################################################################
domain <- read.delim(cfg$infile.domain, header=T, stringsAsFactors=F, sep='\t')
tfbs <- read.table(cfg$infile.tfbs, header=T, stringsAsFactors=F, sep='\t')

####################################################################################################
# switch to output directory
####################################################################################################
setwd(cfg$output.directory)
output.file <- paste(cfg$project.stem, DATE, "MSEA_clust.txt", sep = '_')
#save.image(file = paste(cfg$project.stem, DATE, "Data4EsCum.RData", sep = ''))

# run MSEA-clust on genes
result.gene <- MSEA.clust(mutations.to.test, refseq.length, domain)

write.table(result.gene, file=output.file, row.names=F, col.names=T, quote=F, sep='\t') 

# run MSEA-clust on promoters
result.promoter <- MSEA.clust(prom.ten[1:5,], refseq.len.dummy, tfbs)

####################################################################################################
# run MSEA-domain
####################################################################################################
M1.output <- paste(cfg$project.stem, DATE, 'M1.output.txt', sep = '_')
M2.output <- paste(cfg$project.stem, DATE, 'M2.output.txt', sep = '_')
M3.output <- paste(cfg$project.stem, DATE, 'M3.output.txt', sep = '_')

result.domain <- MSEA.domain(mutations.to.test, domain, refseq.length, M1.output=M1.output, M2.output=M2.output, M3.output=M3.output)

### WRITE SESSION INFO TO FILE #####################################################################
sink(file = paste(DATE, cfg$project.stem, 'Session_Info.txt', sep='_'));

## write memory usage to file
cat('### Memory #########################################################################################\n');
print(gc());

## write process time to file
cat('\n### Time ###########################################################################################\n');
print(proc.time());

## write list of objects
cat('\n### ls() ###########################################################################################\n');
print(ls(pos = 1));

## write sessionInfo to file
cat('\n### Session Info ###################################################################################\n');
print(sessionInfo());

## close the file
sink();

stopCluster(cl)
