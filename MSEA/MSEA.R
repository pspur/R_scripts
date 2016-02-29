MSEA.clust = function(mutations, refseq.length, domains, expand=T){
	
    # hacky way to make code referencing RefSeq.ID work with promoter info
    if (length(mutations)==10) { 
        mutations$RefSeq.ID <- mutations$Symbol
	}
    
    genes_to_test <- unique(mutations$RefSeq.ID)
	
    #expand each row to match number of mutations (freq)
	if (expand) {
        mutations <- mutations[rep(seq_len(nrow(mutations)), mutations$freq), ]
	    row.names(mutations) <- NULL
	}
    
    result <- foreach(k=1:length(genes_to_test)) %dopar% {
        print(paste(genes_to_test[k], ' ', k, '/', length(genes_to_test), sep = ''))
        refseq.ID <- genes_to_test[[k]]
		exonic <- mutations[which(mutations$RefSeq.ID==refseq.ID), ]
		gene.length <- refseq.length[refseq.ID]
		mut_pos <- exonic$mut_pos
        snp <- exonic$snp
        freq <- exonic$freq
	    gene.domains <- domains[(which(domains$refseq.ID==refseq.ID)), ]	
        
		### es.random
		es.random <- vector(length = gene.length*10, mode = 'numeric')
		for(ii in 1:(gene.length*10)){
			Nm <- nrow(exonic)
        
            # randomly select the same number of variants as in input data
			mut_pos.pai <- sample(1:gene.length, Nm, replace=T)
			inc <- 1/length(mut_pos.pai) # same as Nm
			nMut.per.location <- table(mut_pos.pai)
			dec <- 1/(gene.length-length(nMut.per.location))
			inc.1 <- rep(0, gene.length)
			inc.1[as.numeric(names(nMut.per.location))] <- inc*nMut.per.location
			dec.1 <- rep(-dec, gene.length)
			dec.1[as.numeric(names(nMut.per.location))] <- 0
			ss <- inc.1 + dec.1
			es.cum <- cumsum(ss)
			es.pai <- (max(es.cum)-min(es.cum))
			es.random[ii] <- es.pai
		}
		
		### es true
		inc <- 1/length(mut_pos)
		nMut.per.location <- table(mut_pos)
		dec <- 1/(gene.length-length(nMut.per.location))
		inc.1 <- rep(0, gene.length)
		inc.1[as.numeric(names(nMut.per.location))] <- inc*nMut.per.location
		dec.1 <- rep(-dec, gene.length)
		dec.1[as.numeric(names(nMut.per.location))] <- 0
		ss <- inc.1 + dec.1
		es.cum <- cumsum(ss)
		true.es.cum <- es.cum
		es.true <- (max(es.cum)-min(es.cum))
        pvalue <- sum(es.random>=es.true)/(length(es.random)+1)
		nes <- (es.true-mean(es.random))/sd(es.random)
		p1 <- sum(es.random>=es.true)/(length(es.random)+1)
        rv <- c(as.character(exonic$Symbol[1]), refseq.ID, es.true, nes, p1, nrow(exonic), 
                which.min(es.cum), min(es.cum), which.max(es.cum), max(es.cum) )
        list(RefGene=refseq.ID, Symbol=exonic$Symbol[1], es.cum=es.cum, mut_pos=mut_pos, freq=freq, 
             snp=snp, pvalue=pvalue, gene_len=gene.length, domains=gene.domains, rv=rv)
    }

    clust.res <- NULL
    for (i in 1:length(result)) {
        clust.res <- rbind(clust.res, result[[i]]$rv)
    }
    
    clust.res <- as.data.frame(clust.res, stringsAsFactors=FALSE)
    #clust.res <- as.data.frame(matrix(unlist(result), ncol=10), stringsAsFactors = FALSE)
	colnames(clust.res) = c('symbol', 'refseq', 'es', 'nes', 'pvalue', 'nMut', 
	                        'es.min.location', 'es.min', 'es.max.location', 'es.max')
    # adjust for mutiple test
    clust.res$Q <- p.adjust(clust.res$pvalue, method='fdr')
    
    # CONDITIONAL FOR CLUST.RES[Q] IS SET HIGH SO PLOT FUNCTION IS ALWAYS CALLED
    for (i in 1:length(result)) {
        #genes with domains
        if ((clust.res$Q[i]<01.001) & (clust.res$symbol[i]!=clust.res$refseq[i]) & (isTRUE(!nrow(result[[i]]$domains)==0))) {
            output.file <- paste(cfg$project.stem, '_', DATE, '_', result[[i]]$Symbol, '_', result[[i]]$RefGene, sep='')
            generate.mas.plot(expand,
                              result[[i]]$Symbol, 
                              result[[i]]$RefGene, 
                              result[[i]]$gene_len, 
                              result[[i]]$es.cum, 
                              max(clust.res$es.min[i]),
                              max(clust.res$es.max[i]),
                              result[[i]]$mut_pos,
                              result[[i]]$freq,
                              result[[i]]$snp,
                              clust.res$Q[i], 
                              output.file, 
                              result[[i]]$domains)
        }
        #promoters with motifs
        if ((clust.res$Q[i]<01.001) & (clust.res$symbol[i]==clust.res$refseq[i]) & (isTRUE(!nrow(result[[i]]$domains)==0))) {
            output.file <- paste(cfg$project.stem, '_', DATE, '_promoter_', result[[i]]$Symbol, sep='')
            generate.mas.plot(expand,
                              result[[i]]$Symbol, 
                              result[[i]]$RefGene, 
                              result[[i]]$gene_len, 
                              result[[i]]$es.cum, 
                              max(clust.res$es.min[i]),
                              max(clust.res$es.max[i]),
                              result[[i]]$mut_pos,
                              result[[i]]$freq,
                              result[[i]]$snp,
                              clust.res$Q[i], 
                              output.file,
                              result[[i]]$domains)
        }
    }
    #output list of domains containing mutations
    sink(file=paste(cfg$project.stem, '_', DATE, '_Domain_Hits.txt', sep=''))
    for (i in 1:length(result)) {
        doms <- result[[i]]$domains
        # only process results with domains present
        if (isTRUE(!nrow(doms)==0)) {
            symbol <- result[[i]]$Symbol
            rgene <- result[[i]]$RefGene
            ml <- unique(result[[i]]$mut_pos)
            d.hits <- lapply(ml, function(q) {
                for (ii in 1:nrow(doms)) {
                    if (q >= doms$domain.start[ii] & q <= doms$domain.end[ii]) {
                        return(doms$domain.name[ii])
                    }}})
            if (!is.null(unlist(d.hits))) {
                if ('placehold' %in% doms$domain.type) {
                    cat('Domains in',paste0(symbol,' promoter with mutations:'),paste(unlist(unique(d.hits)),collapse=', '),'\n')
                } else {
                    cat('Domains in',paste0(symbol,' (', rgene, ') with mutations:'),paste(unlist(unique(d.hits)),collapse=', '),'\n')
                }    
            } else {
                if ('placehold' %in% doms$domain.type) {
                    cat('No domains in',paste0(symbol,' promoter contain mutations\n'))
                } else {
                    cat('No domains in',paste0(symbol,' (', rgene, ') contain mutations\n'))
                }
            }
        }
    }
    sink()
    
    ## sort the dataframe by Q values
    clust.res <- clust.res[order(clust.res$Q, decreasing=FALSE),]
    return(clust.res)
}

MSEA.domain <- function(mutations, domain, refseq.length, M1.output=NULL, M2.output=NULL, M3.output=NULL){
	genes_to_test <- unique(mutations$RefSeq.ID)
	
	M3.res <- c()
	all.anno <- c()
	M1.res <- c()
	M2.res <- c()
    

	for(k in 1:length(genes_to_test)){
        print(paste(k, '/', length(genes_to_test), sep = ''))
		refseq.ID <- genes_to_test[[k]]
		exonic <- mutations[which(mutations$RefSeq.ID==refseq.ID),]
		gene.length <- refseq.length[refseq.ID]
		mut_pos <- exonic$mut_pos
        
		### if no domain info or amino acid length does not match, skip
		domain.info <- domain[domain$refseq.ID==refseq.ID,]
		domain.info <- domain.info[which(gene.length == as.numeric(domain.info[1,4]) + 1), ]
        
    # skip if there is no domain info available for this gene
		if(nrow(domain.info) < 1) {
            next
        }
		
		### H0
		nMut.per.location <- table(mut_pos)
		n_mut <- rep(0, gene.length)
		n_mut[as.numeric(names(nMut.per.location))] <- nMut.per.location
        
        possibleError <- tryCatch(
		    h0 <- glm.nb(n_mut ~ 1),
            error <- function(e) {
            message(paste("refseq.ID", e))
            } )
        
	    if (inherits(possibleError, "error")) next
        
		########################
		### M1.p
		########################
		domains.X <- c()
		this.M1.res <- c()
        
        # testing for multiple domain in a single gene
		for(k1 in 1:nrow(domain.info)){
			domain.start  <- domain.info$domain.start[k1];
			domain.end    <- domain.info$domain.end[k1];
			domain.length <- domain.end - domain.start
			if( gene.length - domain.length < 10) next  # same as MSEA.clust, so don't do it again
			
			### keep record of domain regions
			X <- rep(0, gene.length )
			X[domain.start:domain.end] <- 1
			domains.X <- rbind(domains.X, X)
			
			if (sum(mut_pos %in% domain.start:domain.end)!=0) {
				h1 <- glm.nb(n_mut ~ 1+as.factor(X) )
				fit <- anova(h0, h1, test='Chisq')
				M1.res      <- rbind(M1.res,      c(domain.info[k1, ], M1.p=fit[2,'Pr(Chi)']) )
				this.M1.res <- rbind(this.M1.res, c(domain.info[k1, ], M1.p=fit[2,'Pr(Chi)']) )
			}
		}
        
        # get the domain which showed more significant p value
		if(!is.null(this.M1.res)) {
			idx <- which.min(as.numeric(this.M1.res[, 'M1.p']))
			M1.p <- min(as.numeric(this.M1.res[, 'M1.p']))
			M1.site <- this.M1.res[idx, 'domain.start']
		} 
		else {
			M1.p <- 1 # default to 1
			M1.site <- 'NULL'
            next
	    }
		
		########################
		### M3.p
		########################
		domains.X.M3 <- apply(domains.X, 2, function(u)sum(u)!=0)
		if (length(unique(domains.X.M3))==1) next
		
		h3 <- glm.nb(n_mut ~ 1 + as.factor(domains.X.M3)) 
		M3.anova = anova(h0, h3, test='Chisq')
		M3.p <- M3.anova[2,'Pr(Chi)']
		
		########################
		### M2.p
		########################	
		cc <- 1
		new <- c(cc)
		for (nn in 2:length(domains.X.M3)) {
			if (domains.X.M3[nn] == domains.X.M3[nn-1]) {
				new <- c(new, cc)
			} 
		    else {
				cc <- cc+1
				new <- c(new, cc)
			}
		}
		new[which(domains.X.M3==0)] <- 0
		new.tab <- table(new)
		
		this.M2.res <- c()
		M2.anova <- list()
		for (k1 in 2:length(new.tab)) {
			tag <- names(new.tab)[k1]
			domains.X.M2 <- ifelse(new==tag, 1, 0)
			
			h2 <- glm.nb(n_mut ~ 1 + as.factor(domains.X.M2))
			
			fit <- anova(h0, h2, test='Chisq')
			M2.res <- rbind(M2.res, c(exonic[1,1], refseq.ID, 
			                          paste(min(which(domains.X.M2==1)), ' - ', 
			                                max(which(domains.X.M2==1)), sep=''), 
			                          fit[2,'Pr(Chi)']))
			this.M2.res <- rbind(this.M2.res, M2.res) 
		}
		if (!is.null(this.M2.res)) {
			idx <- which.min(as.numeric(this.M2.res[,4]))
			M2.p <- min(as.numeric(this.M2.res[,4]))
			M2.site <- 'M2'
		} 
		else {
			M2.p <- 1
			M2.site <- 'NULL'
		}
		########################
		
		results <- c(M1.p, M2.p, M3.p)
		min.site <- c(M1.site, 'M2', 'M3')[which.min(results)]
        min.p <- min(results)
		M3.res <- rbind(M3.res, c(exonic[1,1], refseq.ID, results, min.site, min.p ) )
	}
    colnames(M3.res) <- c('Gene', 'refseq.ID', 'P1', 'P2', 'P3', 'min.site', 'min.p')

	write.table(M1.res, file=M1.output, row.names=F, col.names=F, quote=F, sep='\t')
	write.table(M2.res, file=M2.output, row.names=F, col.names=F, quote=F, sep='\t')
	write.table(M3.res, file=M3.output, row.names=F, col.names=T, quote=F, sep='\t');
}

generate.mas.plot <- function(expand, symbol, refGene, gene.len, mas, y.floor, y.ceiling, 
                              mut_loc, freq, snp, p.adjusted, output.file, gene.info) {
    
    y.floor <- round_any(min(as.numeric(y.floor)), .01, f=floor)
    y.ceiling <- round_any(max(as.numeric(y.ceiling)), .1, f=ceiling)
    x.ceiling <- round_any(gene.len, 100, f=ceiling)
    
    mas <- as.data.frame(mas)
    names(mas)[1] <- 'y'
    
    #ensure mut_loc is in increasing order and snp matches up with it
    if (mut_loc[1] > rev(mut_loc)[1]) {
        mut_loc <- rev(mut_loc)
        snp <- rev(snp)
        freq <- rev(freq)
    }
    
    mut_loc <- as.data.frame(mut_loc)
    names(mut_loc)[1] <- 'x'
    mut_loc$snp <- factor(snp, level=c('INDEL','SNP'))
    mut_loc$freq <- freq
    
    if (expand) {
        output.file <- paste0(output.file, '_expand.png')
        z <- data.frame(table(mut_loc$x))
        count <- apply(z[-1],1,function(u){seq(from=0, to=u-1)})
        if (is.matrix(count)) {
            count <- as.vector(count)
        }
        mut_loc$count <- unlist(count)
        
        mut_loc$interval <- abs(mas[(mut_loc$x)-1,1] - mas[mut_loc$x,1])
        
        mut_loc$sub_interval <- ifelse(mut_loc$freq>1, 
                                       mut_loc$interval/((mut_loc$freq-1)), 
                                       mut_loc$interval/2)
        
        mut_loc$y <- ifelse(mut_loc$freq>1,
                            mas[(mut_loc$x)-1,1] + (mut_loc$count * mut_loc$sub_interval), 
                            mas[(mut_loc$x)-1,1] + mut_loc$sub_interval)
        
    } else {
        # find midpoint of each vertical line created by mutation so can put geom_point there
        #mut_loc$y <- mas[(mut_loc$x)-1,1] + ((mas[(mut_loc$x),1] - mas[(mut_loc$x)-1,1]) / 2)
        
        # put each variant triangle directly on sequence line
        mut_loc$y <- y.floor
        output.file <- paste0(output.file,'_noexpand.png')
    }
      
    palette <- c('#CC79A7','#D55E00','#0072B2','#F0E442','#009E73','#56B4E9','#E69F00','#999999')
    gene.info$y <- y.floor-0.075
    arrow = arrow(angle=15,length=unit(.15,'inches'),end='both',type='closed')
    bottom.dashed <- geom_segment(data=mas,linetype='dashed',size=.25,
                                  aes(x=which.min(y)+3,xend=length(y)+3,y=min(y),yend=min(y)))
    top.dashed <- geom_segment(data=mas,linetype='dashed',size=.25,
                               aes(x=which.max(y)+3,xend=length(y)+3,y=max(y),yend=max(y)))
    #top.dashed <- annotate('segment', linetype='dashed', 
    #                       x=which(mas$y==max(mas[which(mas$y==min(mas$y)):length(mas$y),])),
    #                       xend=length(mas$y),
    #                       y=max(mas$y[which(mas$y==min(mas$y)):length(mas$y)]),
    #                       yend=max(mas$y[which(mas$y==min(mas$y)):length(mas$y)]))
    mes.line <- geom_segment(data=mas, aes(x=length(y)+5,xend=length(y)+5,y=min(y),yend=max(y)), arrow=arrow)
    mes.text <- annotate('text',label='MES',angle=90,size=4,vjust=2,x=gene.len,y=((min(mas)+max(mas))/2))
    seq.line <- annotate('segment', x=1, xend=gene.len, y=y.floor-0.065, yend=y.floor-0.065, color='gray', size=3)
    
    q <- ggplot()
    q <- q + seq.line
    q <- q + mes.line
    q <- q + mes.text
    q <- q + bottom.dashed
    q <- q + top.dashed
    q <- q + theme_classic()
    q <- q + theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
    if ('placehold' %in% gene.info$domain.type) {
        plot.title <- paste0(symbol, ' (promoter), ', format(p.adjusted, digits = 6)) 
        q <- q + geom_rect(data=gene.info, fill='blue', aes(xmin=domain.start,xmax=domain.end,ymin=y,ymax=y+0.020,alpha=zscore))
        q <- q + scale_alpha_continuous(name='Motif\nZ score', range=c(0.5,1.0), limits=c(1,5))
        q <- q + labs(title = plot.title, x='Nucleotide Sequence', y='Mutation Accumulation Score (MAS)')
    } else {
        gene.info$y <- ifelse(gene.info$domain.type=='multi-dom',y.floor-0.085,y.floor-0.065)
        if (length(unique(gene.info$domain.type))==1) {
            gene.info$y <- y.floor-0.075
        }
        plot.title <- paste0(symbol, ' (', refGene, '), ', format(p.adjusted, digits = 6))
        q <- q + geom_rect(data=gene.info, aes(xmin=domain.start,xmax=domain.end,ymin=y,ymax=y+0.020,fill=domain.name))
        q <- q + scale_fill_discrete(name='Domain\nName')    
        q <- q + labs(title = plot.title, x='Amino Acid Sequence', y='Mutation Accumulation Score (MAS)') 
    }
    q <- q + geom_line(data=mas, aes(x=seq_along(y), y=y), color='blue', stat='identity')
    q <- q + geom_point(data=mut_loc, aes(x=x,y=y,color=snp), shape=17, size=1.5) 
    #q <- q + scale_y_continuous(limits=c(y.floor-0.1,y.ceiling),breaks=seq(y.floor,y.ceiling,.1))
    q <- q + scale_y_continuous(limits=c(y.floor-0.1,y.ceiling),
                                breaks=seq(round_any(y.floor,.1,f=floor),y.ceiling,.1))
    if (x.ceiling <= 2000) {
        q <- q + scale_x_continuous(limits=c(1,x.ceiling+50),breaks=seq(0,x.ceiling,100))
    } else if (x.ceiling <= 5000) {
        q <- q + scale_x_continuous(limits=c(1,round_any(x.ceiling+50, 200, f=ceiling)),
                                    breaks=seq(0,round_any(x.ceiling, 200, f=ceiling),200))
    } else if (x.ceiling <= 10000) {
        q <- q + scale_x_continuous(limits=c(1,round_any(x.ceiling+50, 500, f=ceiling)),
                                    breaks=seq(0,round_any(x.ceiling, 500, f=ceiling),500))
    } else if (x.ceiling <= 20000) {
        q <- q + scale_x_continuous(limits=c(1,round_any(x.ceiling+50, 1000, f=ceiling)),
                                    breaks=seq(0,round_any(x.ceiling, 1000, f=ceiling),1000))
    } else {
        q <- q + scale_x_continuous(limits=c(1,round_any(x.ceiling+50, 2000, f=ceiling)),
                                    breaks=seq(0,round_any(x.ceiling, 2000, f=ceiling),2000))
    }
    q <- q + scale_color_manual(drop=T,limits=levels(mut_loc$snp),name='Variant\nType',
                                values=c('springgreen4','red'))#, labels=c('SNP','INDEL'), breaks=c('TRUE','FALSE'))
    q <- q + guides(color = guide_legend(order=1), fill = guide_legend(order=2))
    
    png(output.file,type='cairo',width=1200,height=1200,res=144)
    print(q)
    dev.off()
}
