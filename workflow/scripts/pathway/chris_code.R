library(data.table)
library(dplyr)
library(magrittr)

##' read in lead snps
reader=function(nm) {
    x=fread(nm)
    setnames(x, make.names(names(x)))
    setnames(x, "Gene.s.","genes")
    x[,genes:=gsub(" ","",genes)]
    x
}
iga=reader("iga-250730.tsv")
igg=reader("igg-250730.tsv")
igm=reader("igm-250730.tsv")
sep="," # use this later - it is the character that separates gene names

process=function(iga,igg, igm) {
    iga[,what:="IgA"]
    igg[,what:="IgG"]
    igm[,what:="IgM"]
    x=rbind(iga,igg,igm)
    x[,vid:=paste0(Chromosome,"_",gsub(",","",Position),"_",gsub(".*>","",rsID),"_",gsub(".*:|>.*","",rsID))]
    x[,avid:=paste0(Chromosome,"_",gsub(",","",Position),"_",gsub(".*:|>.*","",rsID),"_",gsub(".*>","",rsID))]

    ## merging colocalised Ig snps & their annotations
    setnames(igcolocs,make.names(names(igcolocs)))
    igcolocs=merge(igcolocs,unique(x[,.(rsID=sub(":.*","",rsID),vid1=vid)]),by.x="First.isotype.s.lead.SNP",by.y="rsID")
    igcolocs=merge(igcolocs,unique(x[,.(rsID=sub(":.*","",rsID),vid2=vid)]),by.x="Second.isotype.s.lead.SNP",by.y="rsID")
    head(igcolocs,2)

    ## decision threshold
    thr=seq(0,1,by=0.01)
    fdr=sapply(thr, function(a) with(igcolocs[PP.H4.abf>a & Filtered==FALSE], mean(1-PP.H4.abf)))
    ## plot(thr,fdr)
    which(fdr < 0.05) # 81 up
    thr[which(fdr < 0.05)] # 0.80 up
    ## hist(igcolocs[Filtered==FALSE]$PP.H4.abf)
    ## use 0.8 to be safe
    setnames(igcolocs,"First.isotype.s.lead.SNP", "rsid1")
    setnames(igcolocs,"Second.isotype.s.lead.SNP", "rsid2")
    igcolocs=igcolocs[PP.H4.abf >= 0.8 & Filtered==FALSE, .(vid1,vid2,rsid1,rsid2,First.isotype,Second.isotype,Pearson.correlation,PP.H4.abf)]
    ## igcolocs # 5, all H4 > 0.93, all IgA/IgM

    ## merge rows of x
    xm=copy(x)
    for(i in 1:nrow(igcolocs)) {
        message(i)
        snps=c(igcolocs$vid1[i],igcolocs$vid2[i])
        w=which(x$vid %in% snps)
        if(length(w)!=2) # unanticipated
            stop()
        cat("merging ", snps[1], " ", snps[2], "\n")
        xm$vid[w[1]]=paste0(unique(snps), collapse=",")
        xm$what[w[1]]=paste0(unique(xm$what[w]), collapse=",")
        xm$Beta[w[1]]=paste0(unique(xm$Beta[w]), collapse=",")
        xm$genes[w[1]]=unique(xm$genes[w])  %>% strsplit(., sep)  %>% unlist()  %>% unique()  %>% paste(., collapse=sep)
        xm$Novel[w[1]]=any(xm$Novel[w])
        xm=xm[-w[2]]
    }
    xm
}


## run and plot pathway analysis
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(cowplot)
f=function(kk) {
    tmp=as.data.table(kk)[qvalue < 0.05]
    tmp=tmp[order(-qvalue)]
    tmp[,Description:=factor(Description,levels=Description)]
                                        # create a dotplot of tmp showing x=FoldEnrichment, y=Description, and the fill of the point p.adjust
    cols=enrichplot:::get_enrichplot_color(2)
    ggplot(tmp, aes(x=FoldEnrichment,y=Description,fill=-log10(qvalue),size=Count)) +
        geom_point(shape=21) +
        scale_fill_gradient(low=cols[2], high=cols[1]) +
        scale_x_continuous(limits=c(1,max(tmp$FoldEnrichment)+.5)) +
        theme_cowplot() +
        labs(y="") +
        background_grid(major="y")#, breaks=seq(1,13,2)) 
}

## xm$genes[w[1]]=unique(xm$genes[w])  %>% strsplit(., sep)  %>% unlist()  %>% unique()  %>% paste(., collapse=sep)

run_react=function(entrez) {
    rct <- enrichPathway(gene=entrez[,"ENTREZID"],
                         pvalueCutoff = 0.05, readable=TRUE)
    print(as.data.table(rct)[qvalue < 0.05,.(ID,Description,Count,FoldEnrichment,qvalue)])
    print(f(rct))
    invisible(rct)
}
run_kk=function(entrez) {
   kk <- enrichKEGG(gene      = entrez[,"ENTREZID"],
                 organism     = 'hsa',
                 qvalueCutoff = 0.1)
   print(as.data.table(kk)[qvalue < 0.05,.(ID,Description,Count,FoldEnrichment,qvalue)])
   print(f(kk))
   invisible(kk)
}

## all genes
get_entrez=function(xm) {
    symbol=strsplit(xm$genes,sep)  %>% unlist()  %>% unique() 
    entrez=bitr(symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
}


x=process(iga,igg,igm)
kk=run_kk(get_entrez(x))
rct=run_react(get_entrez(x))


## make data structure and group variants by pathway analysis
library(pheatmap)
getM=function(kk, trans=TRUE) {
    dt=as.data.table(kk)[qvalue < 0.05]
    if(trans) {
        entrez=get_entrez(x)
        getsym=function(entrezid) 
            entrez[ entrez[,"ENTREZID"] %in% entrezid, "SYMBOL"]
        memgenes=strsplit(dt$geneID,"/") %>%
            lapply(., getsym)
    } else {
        memgenes=strsplit(dt$geneID,"/")
    }
    allgenes=unlist(memgenes) %>% unique()
    ## maps variants 2 genes
    G=matrix(0,nrow(x),length(allgenes),dimnames=list(x$rsID,allgenes))
    Gi=strsplit(gsub(" ","",x$genes),",") 
    for(i in 1:nrow(x)) 
        G[i,]=as.numeric(allgenes %in% Gi[[i]])
    ## maps genes 2 pathways
    P=matrix(0,length(allgenes),length(memgenes),dimnames=list(allgenes,dt$Description))
    for(i in seq_along(memgenes))
        P[,i]=as.numeric(allgenes %in% memgenes[[i]])
    ## combine
    M= G %*% P # will count through how many genes a variant is linked to a pathway
    M=pmin(M,1) # binarises this
    M=M[rowSums(M)>0,colSums(M)>4]
}

    
g=function(kk,trans=TRUE, ...) {
    M=getM(kk, trans=trans)
                                        # M=t(M)
    on=colnames(M)
    cann=x[rsID %in% rownames(M),.(rsID,IgA=grepl("IgA",what), IgG=grepl("IgG",what), IgM=grepl("IgM",what))] %>%
        as.data.frame()
    cann$IgA %<>% as.numeric()
    cann$IgM %<>% as.numeric()
    cann$IgG %<>% as.numeric()
    rownames(cann)=cann$rsID
    cann=cann[,c("IgA","IgG","IgM")]

    om=rownames(M)  
    summary(nchar(om))

    m=sprintf("%-11s\t%-22s",sub(":.*","",om), x$genes[ match(om, x$rsID) ])
    rownames(M)=m
    rownames(cann)=m
    print(dim(M))
    pheatmap(M,show_rownames=TRUE,show_colnames=TRUE,
             color=colorRampPalette(c("grey90", "grey20"))(50),
             clustering_method="ward.D",
             ## labels_col=substring(colnames(M),1,nchar(colnames(M))),
             labels_col=stringr::str_wrap(colnames(M),50, exdent=6),
             annotation_row=cann,
                                        # annotation_row=rann,
             cutree_row=7, cutree_col=5,
             legend=FALSE, annotation_legend=FALSE, ...)

}

g(rct, FALSE, main="Reactome pathways")
g(kk,main="Kegg pathways")
## g(mf,FALSE)
## g(bp,FALSE)

## combine kk, rct
getsym=function(entrezid) 
    entrez[ entrez[,"ENTREZID"] %in% entrezid, "SYMBOL"]
entrez=get_entrez(x)
memgenes=strsplit(kk$geneID,"/") %>%
    lapply(., getsym)  %>%
    sapply(., paste, collapse="/")
kk2=as.data.table(kk)
kk2$geneID=memgenes
m=rbind(kk2, as.data.table(rct), fill=TRUE)
g(m, FALSE, main="Reactome + Kegg", filename="combined_pathways.tiff", width=12, height=10, fontsize=8)
