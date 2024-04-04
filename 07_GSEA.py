library(ggplot2)
library(dplyr)
library(fgsea)

makeGeneList <- function(filename){
    gl <- readr::read_tsv(filename)
    #y <- grepl('^RPS|^RPL|^MRPL|^MRPS|^MT-', gl$X1)
    y <- grepl('XXXXXXXX', gl$X1)
    gl <- gl[!y, ]
    gl <- gl %>% dplyr::select(X1, logfoldchanges , pvals)
    gl$neglog10pval <- -log10(gl$pvals)
    rank <- unlist(gl$neglog10pval*sign(gl$logfoldchanges))
    rank[rank == Inf] = 300
    rank[rank == -Inf] = -300
    names(rank) <- gl$X1
    rank <- rev(sort(rank))
    gl <- rank
}

celltypes = c('ASDC',
 'B intermediate',
 'B memory',
 'B naive',
 'CD14 Mono',
 'CD16 Mono',
 'CD4 CTL',
 'CD4 Naive',
 'CD4 Proliferating',
 'CD4 TCM',
 'CD4 TEM',
 'CD8 Naive',
 'CD8 Proliferating',
 'CD8 TCM',
 'CD8 TEM',
 'cDC2',
 'dnT',
 'Eryth',
 'gdT',
 'HSPC',
 'ILC',
 'MAIT',
 'NK Proliferating',
 'NK_CD56bright',
 'NK',
 'pDC',
 'Plasmablast',
 'Platelet',
 'Treg'
             )
groups = c('DHF', 'DF')

comparisons = list()
for (c in celltypes){
    for (g in groups){
        comparisons[[c]][[g]] = makeGeneList(paste0('/home/jovyan/AsymDENV/08_DGE/wilcoxon/',g,'_',c, '_wilcoxon.txt'))
    }
}

parse_gmt <- function(file, header = TRUE, sep = "\t", ...) {

  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)

  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }

  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
}

h <- as.list(parse_gmt('h.all.v7.4.symbols.gmt'))

c7 <- as.list(parse_gmt('c7.immunesigdb.v7.4.symbols.gmt'))

c5bp <- as.list(parse_gmt('c5.go.bp.v7.4.symbols.gmt'))

c2 <- as.list(parse_gmt('c2.cp.kegg.v7.4.symbols.gmt'))

res = list()
for (c in celltypes){
    for (g in groups){
    
        res[[c]][[g]] <- fgsea(pathways = h, stats = comparisons[[c]][[g]], nperm = 10000, minSize = 0, maxSize =1000)
    }
}

for (c in celltypes){
    for (g in groups){
        res[[c]][[g]] <- as.data.frame(res[[c]][[g]])
    }
}

result = res

# names(result) <- gsub('.csv', "",files)
for(i in 1:length(comparisons)){
    result[[i]] <- lapply(result[[i]], function(x){
    x$ranking <- -log10(x$pval)*sign(x$NES) 
    x <- x[order(x$pathway), ]
    return(x)
})
}

result <- lapply(result ,function(x){
    x[['DHF']]$group = "DHF"
    x[['DF']]$group = "DF"

    return(x)
})

result2 <- lapply(result, function(x) {
    y <- do.call(rbind, x)
    return(y)
})

plotGSEA_Hallmark <- function(gsea, group_ref = NULL, cols = NULL, newlabels = NULL, keep_significant_only = TRUE) {
    require(ggplot2)
    gsea$NES[which(is.na(gsea$NES))] <- 0
    gsea$pval[which(is.na(gsea$pval))] <- 1
    gsea$padj[which(is.na(gsea$padj))] <- 1
    gsea$ranking[which(is.na(gsea$ranking))] <- 0
    gsea <- gsea[order(gsea$ranking),]      
    gsea_spl <- split(gsea, gsea$group)
    if(!is.null(group_ref)){
        gsea_spl[[group_ref]] <- gsea_spl[[group_ref]][order(gsea_spl[[group_ref]]$ranking),]
        gsea_spl[[group_ref]]$ranking <- gsea_spl[[group_ref]]$ranking*999
    } else {
        gsea_spl[[2]] <- gsea_spl[[2]]$ranking*999
    }
    names(gsea_spl) <- NULL

    gsea <- do.call(rbind, gsea_spl)
    
    if (keep_significant_only){
        gseax <- split(gsea, gsea$pathway)
        for (i in 1:length(gseax)){
            if (all(gseax[[i]]$pval >= 0.05)|all(gseax[[i]]$padj >=0.05)){
                gseax[[i]] <- NA        
            }
        }
        gseax <- gseax[!is.na(gseax)]
        gsea <- do.call(rbind, gseax)
    }
    gsea <- gsea[order(gsea$ranking), ]
    gsea$pathway <- gsub("KEGG_|", "", gsea$pathway)
    gsea$group[which(gsea$pval >= 0.05 & gsea$padj >= 0.05)] <- 'NotSig'
    gsea$group[which(gsea$pval < 0.05 & gsea$padj >= 0.05)] <- 'NotSig'
    gsea$group[which(gsea$pval >= 0.05)] <- 'NotSig'
    gsea$group <- factor(gsea$group, levels = c('NotSig', 'DF', 'DHF'))

    x_lim_min <- abs(ceiling(min(-log10(gsea$padj))))
    x_lim_max <- abs(ceiling(max(-log10(gsea$padj))))
        
    if(x_lim_min > x_lim_max){
        xval1 <- x_lim_min * -1
        xval2 <- x_lim_min
    } else {
        xval1 <- x_lim_max * -1
        xval2 <- x_lim_max
    }

    if(!is.null(cols)){
        gg_color_hue <- function(n) {
            hues = seq(15, 375, length = n + 1)
            hcl(h = hues, l = 65, c = 100)[1:n]
        }
        cols. = gg_color_hue(dplyr::n_distinct(gsea$group, na.rm = TRUE))
    } else {
        cols. = cols
    }    
    
    g <- ggplot(gsea, aes(x = -log10(padj)*sign(NES), y = reorder(pathway, ranking), col = group, size = abs(NES))) + 
        geom_point() + 
        labs(x = expression(paste("Signed", " -log" ["10"], "adjusted pval")), y = "KEGG") +
        theme_bw() +
        geom_vline(xintercept = 0) +
        geom_vline(xintercept = -log10(0.05)) +
        geom_vline(xintercept = -log10(0.05)*-1) +
        xlim(xval1, xval2) +
        scale_size_area(oob = scales::squish, max_size = 3, limits = c(0,2)) +
        theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(), 
            axis.line = element_blank(), 
            axis.ticks = element_blank())
    
    g$data <- g$data[order(g$data$group, na.last = TRUE), ]
    return(g)
}

for(i in 1:29){
tryCatch({p <- plotGSEA_Hallmark(result2[[i]], group_ref = "DF", keep_significant_only = TRUE) + scale_color_manual(values = c('#e7e7e7' , '#fd8d3c','#e31a1c'), drop=FALSE) + theme(axis.text.y=element_text(size=5))
ggsave(paste0('08_DGE/gsea',names(result2)[i],".png"), plot = p, w = 12)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#save results
GSEA.list <- list()
for(i in 1:27){
tryCatch({
df <- result2[[i]] 
df$CellType <- celltypes[[i]]
row.names(df) <- with(df, paste(pathway, group, CellType, sep = "_"))
GSEA.list[[i]] <- df
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}

GSEA_df <- do.call(rbind, GSEA.list)

getwd()

#library(data.table)
#fwrite(GSEA_df, 'GSEA_KEGG.csv')
GSEA_df <- read.csv('/home/jovyan/AsymDENV/08_DGE/GSEA_results_tables/GSEA_Hallmark.csv', header=T)
#head(GSEA_df)



library(stringr)

GSEA_df 

#extract specific pathway
pathway <- GSEA_df %>% filter(pathway == 'HALLMARK_INFLAMMATORY_RESPONSE')
pathway[1,]

plotGSEA <- function(gsea, group_ref = NULL, cols = NULL, newlabels = NULL, keep_significant_only = TRUE) {
    require(ggplot2)
    gsea$NES[which(is.na(gsea$NES))] <- 0
    gsea$pval[which(is.na(gsea$pval))] <- 1
    gsea$padj[which(is.na(gsea$padj))] <- 1
    gsea$ranking[which(is.na(gsea$ranking))] <- 0
    gsea <- gsea[order(gsea$ranking),]      
    gsea_spl <- split(gsea, gsea$group)
    if(!is.null(group_ref)){
        gsea_spl[[group_ref]] <- gsea_spl[[group_ref]][order(gsea_spl[[group_ref]]$ranking),]
        gsea_spl[[group_ref]]$ranking <- gsea_spl[[group_ref]]$ranking*999
    } else {
        gsea_spl[[2]] <- gsea_spl[[2]]$ranking*999
    }
    names(gsea_spl) <- NULL

    gsea <- do.call(rbind, gsea_spl)
    
    if (keep_significant_only){
        gseax <- split(gsea, gsea$pathway)
        for (i in 1:length(gseax)){
            if (all(gseax[[i]]$pval >= 0.05)|all(gseax[[i]]$padj >=0.5)){
                gseax[[i]] <- NA        
            }
        }
        gseax <- gseax[!is.na(gseax)]
        gsea <- do.call(rbind, gseax)
    }
    gsea <- gsea[order(gsea$ranking), ]
    gsea$pathway <- gsub("Hallmark_|", "", gsea$pathway)
    gsea$group[which(gsea$pval >= 0.05 & gsea$padj >= 0.05)] <- 'NotSig'
    gsea$group[which(gsea$pval < 0.05 & gsea$padj >= 0.05)] <- 'NotSig'
    gsea$group[which(gsea$pval >= 0.05)] <- 'NotSig'
    gsea$group <- factor(gsea$group, levels = c('NotSig', 'DF', 'DHF'))

    x_lim_min <- abs(ceiling(min(-log10(gsea$padj))))
    x_lim_max <- abs(ceiling(max(-log10(gsea$padj))))
        
    if(x_lim_min > x_lim_max){
        xval1 <- x_lim_min * -1
        xval2 <- x_lim_min
    } else {
        xval1 <- x_lim_max * -1
        xval2 <- x_lim_max
    }

    if(!is.null(cols)){
        gg_color_hue <- function(n) {
            hues = seq(15, 375, length = n + 1)
            hcl(h = hues, l = 65, c = 100)[1:n]
        }
        cols. = gg_color_hue(dplyr::n_distinct(gsea$group, na.rm = TRUE))
    } else {
        cols. = cols
    }    
    
    g <- ggplot(gsea, aes(x = -log10(padj)*sign(NES), y = reorder(CellType, desc(CellType)), col = group, size = abs(NES))) + 
        geom_point() + 
        labs(x = expression(paste("Signed", " -log" ["10"], "adjusted pval")), y = "CellType") +
        theme_bw() +
        geom_vline(xintercept = 0) +
        geom_vline(xintercept = -log10(0.05)) +
        geom_vline(xintercept = -log10(0.05)*-1) +
        xlim(xval1, xval2) +
        #scale_size_area(oob = scales::squish, max_size = 4, limits = c(0,3)) +
        theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(), 
            axis.line = element_blank(), 
            axis.ticks = element_blank())
        
    
    g$data <- g$data[order(g$data$group, na.last = TRUE), ]
    return(g)
}

#extract specific pathway
pathway <- GSEA_df %>% filter(pathway == 'HALLMARK_TNFA_SIGNALING_VIA_NFKB')
pathway[1,]

plotGSEA(pathway, group_ref = "DF", keep_significant_only = TRUE) + scale_color_manual(values = c('#e7e7e7' , '#fd8d3c','#e31a1c'), drop=FALSE) + theme(axis.text.y=element_text(size=8))#Only select significant celltypes


pathway.filter <- pathway %>% filter(CellType %in% c('B naive',
                                                    'CD14 Mono',
                                                    'B memory',
                                                    'B intermediate',
                                                    'NK',
                                                    'MAIT',
                                                    'CD4 TCM',
                                                    'CD4 TEM',
                                                    'CD16 Mono'))

plotGSEA(pathway.filter, group_ref = "DF", keep_significant_only = TRUE) + scale_color_manual(values = c('#e7e7e7' , '#fd8d3c','#e31a1c'), drop=FALSE) + theme(axis.text.y=element_text(size=8))

ggsave('GOBP_FC_EPSILON_RECEPTOR_SIGNALING_PATHWAY.pdf', width = 5, height = 4)

