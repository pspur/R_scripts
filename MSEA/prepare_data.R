### The following contains the procedures to reproduce MSEA input data.
### All files used should reside in MSEA/data

setwd('/home/paulspur/R_scripts/MSEA/')
#setwd('/Users/Paul/Documents/MSEA/')
library(doParallel)
library(plyr)
library(ggplot2)
library(grid)
library(dplyr)

source('./bin/MSEA.R')

###################################################################################################
### Procedure to create 'protein.gbk.regions.symbol.txt' 

### 1. Download 'ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/protein/protein.gbk.gz'
### This is a flat file of all human refseq proteins, including predicted (XP_) ones

### 2. Run 'makeProteinID_list.py'. This creates 'pids.txt', a file containing all 
### established (not predicted) protein ids, based on the flat file.

### 3. Go to http://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi and upload
### 'pids.txt' then submit the search. This step has a run time of ~2 hours.

### 4. Save the resulting file locally. Use default name: hitdata.txt.

### 5. Run 'makeDomainTable.py'. This creates protein.gbk.regions.symbol.txt 
### from the flat file and hitdata.txt.

### 6. Place protein.gbk.regions.symbol.txt in MSEA/input

###################################################################################################
### Procedure to create 'ref.gene.final.txt'.

### 1. Download 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz'
### This contains hg19 mRNA info, including exon intervals.

ref.gene <- read.table('./data/refGene.txt', sep='\t', header=F, stringsAsFactors=F)

# remove entries with non standard chroms
ref.gene <- ref.gene[grep('_', ref.gene$V3, invert=T), ]

write.table(ref.gene, file='./input/ref.gene.final.txt', quote=F, sep='\t')

###################################################################################################
### Procedure to create 'tfbs.plot_ready.txt', the promoter equivalent of domains

### 1. Download 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/tfbsConsFactors.txt.gz'
### and 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/tfbsConsSites.txt.gz'

tf_sites.raw <- read.table('./data/tfbsConsSites.txt', sep='\t', header=F, stringsAsFactors=F)
tf_factors.raw <- read.table('./data/tfbsConsFactors.txt', sep='\t', header=F, stringsAsFactors=F)
tf_factors <- filter(tf_factors.raw, V3=='human')
tf_sites <- tf_sites.raw[ ,c(-1,-6)]
tf_sites <- filter(tf_sites, V5 %in% tf_factors$V1)    

# assign tfbs to promoter regions, calculate motif starts and ends relative to promoter start
# run time 6-9 hours
system.time(tfbs_gene_info <- apply(tf_sites, 1, match.tfbs.to.promoter))[3]

# assign returned values to respective columns, messy b/c tfbs_gene_info is list of lists of (sometimes) list
tf_sites$symbol <- lapply(tfbs_gene_info, function(u) { u[[1]]})
tf_sites$domain.start <- lapply(tfbs_gene_info, function(u) { u[[2]]})
tf_sites$domain.end <- lapply(tfbs_gene_info, function(u) { u[[3]]})

# remove all sites that didn't match up with a promoter region
tf_sites<- filter(tf_sites, symbol!='n')

# convert new columns from lists of lists to vectors of characters
tf_sites$symbol <- sapply(tf_sites$symbol, paste0, collapse=",")
tf_sites$domain.start <- sapply(tf_sites$domain.start, paste0, collapse=",")
tf_sites$domain.end <- sapply(tf_sites$domain.end, paste0, collapse=",")

# store new columnns as lists of characters to facilitate expanding any with multiple entries
a <- strsplit(tf_sites$symbol, ',')
b <- strsplit(tf_sites$domain.start, ',')
c <- strsplit(tf_sites$domain.end, ',')

# create final tfbs df with all rows that contain multiple gene entries expanded
tfbs <- data.frame(symbol = unlist(a),
                   refseq.ID = unlist(a),
                   protein.ID = NA,
                   length = 1000, 
                   domain.start = as.integer(unlist(b)),
                   domain.end = as.integer(unlist(c)),
                   domain.source = NA,
                   domain.name = rep(tf_sites$V5, sapply(a, length)),
                   domain.anno = NA,
                   domain.type = 'placehold',
                   zscore = rep(tf_sites$V8, sapply(a, length)),
                   stringsAsFactors=F)

write.table(tfbs, file='./input/tfbs.plot_ready.txt', quote=F, sep='\t')

###################################################################################################
### Procedure to create 'prom.plot_ready.txt', a file containing variants in gene promoter regions

### 1. Run 'remove_alt_alleles.py [input annotated vcf] [output]'

### 2. Repeat step 1 four more times on each subsequent output
### This creates an annotated vcf file with all alternate alleles removed.  

### 3. Execute the following command to collect all variants in promoter regions:
### grep " upstream" /path/to/annovaredvcf.txt > 20160105_sim_wgs.hg19_multianno_1kb_promoters.txt
### (That is a tab followed by 'upstream' in the grep expression.)
### By default annovar marks any variants within 1000 bases upstream of genes as "upstream",
### so redefining promoter regions using exactly this procedure requires rerunning 
### annovar with an altered promoter length parameter.

prom.raw <- read.table('./data/20160105_sim_wgs.hg19_multianno_1kb_promoters.txt', 
                       sep='\t', header=F, stringsAsFactors=F, quote=NULL)
prom.pass <- filter(prom.raw, V162=='PASS')

# drop downstream genes from 'upstream;downstream' rows
prom.pass$V7 <- as.character(lapply(prom.pass$V7, function(u) { strsplit(u, split=';')[[1]][1]}))

prom.pass$V1 <- paste0('chr', prom.pass$V1)

# slice range is all sample genotype info columns
gt.fields <- prom.pass[ ,165:214]

# calculate each variant's frequency to expand later
prom.pass$freq <- var.freq(gt.fields)

# keep chr, start, end, ref, alt, gene, freq
prom.pass <- prom.pass[ ,c(1:5,7,215)]

# expand rows with multiple genes
a <- strsplit(prom.pass$V7, ',')
prom.pass <- data.frame(V1 = rep(prom.pass$V1, sapply(a, length)), 
                        V2 = rep(prom.pass$V2, sapply(a, length)),
                        V3 = rep(prom.pass$V3, sapply(a, length)),
                        V4 = rep(prom.pass$V4, sapply(a, length)),
                        V5 = rep(prom.pass$V5, sapply(a, length)),
                        freq = rep(prom.pass$freq, sapply(a, length)),
                        V7 = unlist(a),
                        stringsAsFactors=F)

# remove genes from promoters that don't exist in refGene
prom.pass <- filter(prom.pass, V7 %in% ref.gene$V13)

# calculate mut locs referenced to gene start locations and get strand info
system.time(prominfo <- apply(prom.pass, 1, get.mutloc.and.strand))[3]
prominfo <- matrix(unlist(prominfo), nrow=length(prominfo), byrow=T)
prominfo <- data.frame(prominfo, stringsAsFactors=F)
prom.pass$mut_pos <- as.integer(prominfo$X1)
prom.pass$strand <- prominfo$X2

colnames(prom.pass) <- c('chrom', 'start', 'end', 'ref', 'alt', 'freq', 'Symbol', 'mut_pos', 'strand')

# remove genes from promoters that had mismatched refGene info (mut_loc = -1)
prom.pass <- filter(prom.pass, mut_pos!=-1)

nt <- c('A','T','C','G')
ref.snp <- prom.pass$ref %in% nt
alt.snp <- prom.pass$alt %in% nt
prom.pass$snp <- ref.snp & alt.snp
prom.pass$snp <- ifelse(prom.pass$snp=='TRUE','SNP','INDEL')

write.table(prom.pass, file='./input/prom.plot_ready.txt', quote=F, sep='\t')

###################################################################################################
### Procedure to create the various variant input files:

### 1. Execute 'expand.row.pl' on the de-alternate-alleled annotated variant text file.
### This produces 6 files:
### 'all_nonsilent_SNVs.txt', 'deleterious_nonsilent_SNVs.txt',
### 'all_nonsilent_SNVs_Indels.txt', 'deleterious_nonsilent_SNVs_Indels.txt',
### 'all_silent_SNVs.txt', 'all_silent_benign_nonsilent_SNVs.txt'

### 2. Run the following script, files created are placed in ./input 

nsSNV <- read.table('./data/all_nonsilent_SNVs.txt', sep='\t', header=F, stringsAsFactors=F)
nsSNVIndel <- read.table('./data/all_nonsilent_SNVs_Indels.txt', sep='\t', header=F, stringsAsFactors=F)
sbnsSNV <- read.table('./data/all_silent_benign_nonsilent_SNVs.txt', sep='\t', header=F, stringsAsFactors=F)
sSNV <- read.table('./data/all_silent_SNVs.txt', sep='\t', header=F, stringsAsFactors=F)
dnsSNV <- read.table('./data/deleterious_nonsilent_SNVs.txt', sep='\t', header=F, stringsAsFactors=F)
dnsSNVIndel <- read.table('./data/deleterious_nonsilent_SNVs_Indels.txt', sep='\t', header=F, stringsAsFactors=F)

process.vars(nsSNV,'all_nonsilent_SNVs')
process.vars(nsSNVIndel,'all_nonsilent_SNVs_Indels')
process.vars(sbnsSNV,'all_silent_benign_nonsilent_SNVs')
process.vars(sSNV,'all_silent_SNVs')
process.vars(dnsSNV,'deleterious_nonsilent_SNVs')
process.vars(dnsSNVIndel,'deleterious_nonsilent_SNVs_Indels')

process.vars <- function(variants,fname) {
    
    #################### check this logic, it ignores rows with multiple genes in gene field
    var_info <- apply(variants, 1, function(u) {v <- strsplit(u[10], split = ',')[[1]]})
    var_info <- lapply(var_info, function(w) {rv <- strsplit(w, split = ':')})
    mutations <- as.data.frame(matrix((unlist(var_info)), ncol=5, byrow=T), stringsAsFactors=F)
    mutations$freq <- variants$V11
    
    ### why comment these out?
    #mutations = mutations[which(mutations$DNA_change!=""), ]  ### remove occasional errors
    #mutations = mutations[!is.na(match(mutations$RefSeq.ID, names(refseq.length))), ]  ### remove refseq genes with no gene length
    
    nt <- c('A','T','C','G')
    ref.snp <- variants$V4 %in% nt
    alt.snp <- variants$V5 %in% nt
    mutations$snp <- ref.snp & alt.snp
    mutations$snp <- ifelse(mutations$snp=='TRUE','SNP','INDEL')
    colnames(mutations) <- c('Symbol', 'RefSeq.ID', 'Exon.index', 'DNA_change', 
                             'amino_acid_change', 'freq', 'snp')
    
    ### extract mutation positions
    mut_pos <- as.numeric(unlist(lapply(mutations$amino_acid_change, function(u){ 
        if(grepl("_", u)){
            u1 <- substr(u, 3, regexec("_", u)[[1]][1]-1  )
            m <- regexec("[0-9]+", u1)
            a <- regmatches(u1,m)
        } else {
            m <- regexec("[0-9]+", u)
            m <- regmatches(u,m)
        } } ) ))
    mutations$mut_pos <- mut_pos
    
    write.table(mutations, file=paste('./input/', fname, '.plot_ready.txt', sep=''), quote=F, sep='\t')
}

# free up memory
gc()
