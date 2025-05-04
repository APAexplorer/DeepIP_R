#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual element_text facet_wrap theme element_text
#' @importFrom IRanges width start end subsetByOverlaps
#' @importFrom movAPA readPACds faFromPACds
#' @importFrom utils read.table write.table read.csv write.csv
NULL

options(stringsAsFactors=F)


## --------- utils ----------
#' formatPathFileName converts Windows file name to Linux.
#' @param aPath a string
#' @param addBar whether to add / in the end of the string, default is FALSE
#' @param noExt default is FALSE
#' @param onlyFileName default is FALSE, T to get only filename
#' @return a new string in Linux format
#' @examples
#' formatPathFileName('data\\a.bam') # "data/a.bam"
#' formatPathFileName('data/a.bam') # "data/a.bam"
#' formatPathFileName('data/a.bam', addBar=T) # "data/a.bam/"
#' formatPathFileName('data/a.bam/', addBar=T) # "data/a.bam/"
#' formatPathFileName('data/a.bam/', addBar=T, onlyFileName = TRUE) # a.bam
#' formatPathFileName('data/a.bam/', noExt=T, onlyFileName = TRUE) # a
#' @export
formatPathFileName <- function(aPath, addBar=FALSE, noExt=FALSE, onlyFileName=FALSE) {
  new_path <- gsub(pattern = "\\", replacement = "/", x = aPath, fixed = TRUE)
  if (addBar) {
    if (substr(new_path, nchar(new_path), nchar(new_path)) != '/')
      new_path=paste0(new_path, '/')
  }

  if (onlyFileName) new_path=basename(new_path) #"data/a.bam/" --> "a.bam"

  if (noExt) return(tools::file_path_sans_ext(new_path))

  return(new_path)
}

#' trimPAseqs trims each sequence in fa files with start~stop
#' @param faFiles fa file names.
#' @param start start position for all sequences.
#' @param stop stop position for all sequences, final sequence will be from start(include)..end (include).
#' @param sufix output file name be like <faFiles><_trim>.fa
#' @param ofiles must the same length as faFiles. If provided, then will ignore sufix.
#' @param outputDir if NULL, then output each file in the dir of corresponding `faFiles`.
#' @param basePlot if TURE, to plot first 6 output files' base compositions.
#' @return output file names.
#' @examples
#' \dontrun{
#' ## start=21, stop=220 to trim 200nt from a 241nt seq with PA in 101
#' trimPAseqs('3utr_nonPA_seq_train_per70_deepPASTA_5000T5000F_train.fa',
#'            sufix='_trim200PA101',
#'            start=21, stop=220, basePlot=TRUE)
#' }
#' @export
trimPAseqs<-function(faFiles, start=NULL, stop=NULL, sufix='_trim',
                     ofiles=NULL, outputDir=NULL, basePlot=FALSE) {

  if (is.null(start)) start=1

  if (is.null(ofiles)) {
    ofiles=gsub('.fa', '', faFiles)
    ofiles=paste0(ofiles, sufix, '.fa')
  } else {
    if (length(ofiles)!=length(faFiles))
      stop("ofiles should be the same length as faFiles!")
  }

  plotFile=paste0(sufix,'.base_profile')

  if (!is.null(outputDir)) {
    if (!dir.exists(outputDir)) dir.create(outputDir)
    ofiles=basename(ofiles)
    ofiles=paste0(formatPathFileName(outputDir, addBar=TRUE), ofiles)
    plotFile=paste0(formatPathFileName(outputDir, addBar=TRUE), paste0(sufix,'.base_profile'))
  }

  for (f in faFiles) {
    cat('trim', f, '\n')
    fa=seqinr::read.fasta(f, forceDNAtolower = F, as.string = TRUE)
    fa2=fa
    for (i in 1:length(fa)) {
      if (is.null(stop)) stop=length(fa[i])
      fa2[i]=substr(fa[i], start=start, stop=stop)
    }
    seqinr::write.fasta(sequences = fa2,
                        name=names(fa),
                        nbchar = 1000,
                        file.out = ofiles[faFiles==f],
                        open="w")

  }
  if (basePlot)
    movAPA::plotATCGforFAfile(ofiles[1:min(6, length(ofiles))], filepre = plotFile, mergePlots=TRUE)
  return(ofiles)
}


#' resetNonPATitle resets sequence title in the fa file to make sure the title reflects the position of the sequence
#' @param faFiles fa file names.
#' Such fa files are sequences like "nonPA-pos", or "IP=".
#' @param sufix output file name be like <faFiles><_ts>.fa
#' @return output file names.
#' @export
resetNonPATitle<-function(faFiles,
                            sufix="_ts") {


  .reset<-function(title) {
    t1=unlist(strsplit(title, split=';'))
    st=as.numeric(t1[3])
    end=as.numeric(t1[4])
    pos=as.numeric(t1[5])
    if (anyNA(c(st,end,pos))) stop("error fa title, should be like: >PA51504:chr11;-;84093923;84098722;nonPA-pos=3901... OR IP=3901...")

    if (t1[istrand]=='+') {
      newpos=st+pos-1
    } else {
      newpos=end-pos+1
    }
    newtitle=paste(t1[1],t1[2], newpos, title, sep=";")
    return(newtitle)
  }

  ofiles=gsub('.fa', '', faFiles)
  ofiles=paste0(ofiles, sufix, '.fa')

  for (f in faFiles) {
    fa=seqinr::read.fasta(f, forceDNAtolower = F, as.string = TRUE)

    # >PA51504:chr11;-;84093923;84098722;nonPA-pos=3901;TATAAA:0;TATAAA:0

    # >PA46953:chr10;-;92531949;92534377;IP=1232;ACTAAA:0;ACTAAA:0
    titles=names(fa)

    t1=unlist(strsplit(titles[1], split=';'))
    istrand=which(t1 %in% c('+','-'))[1]

    if (grepl('nonPA-pos=', titles[1], fixed=TRUE)) {
      str='nonPA-pos='
    } else if (grepl('IP=', titles[1], fixed=TRUE)) {
      str='IP='
    } else {
      stop("nonPA-pos= or IP= not in the title of fa file", f,"\n")
    }


    titles=gsub(str,'',titles, fixed=TRUE)

    titles=unlist(lapply(titles, .reset))

    seqinr::write.fasta(sequences = fa,
                        names = titles,
                        nbchar = 1000,
                        file.out =ofiles[faFiles==f],
                        open="w")

    cat('reset title -->', ofiles[faFiles==f], '\n' )

  }

}


#' sampleNseqs randomly select N seqs from given faFiles and output each file.
#' @param faFiles fa file names.
#' @param N number of sequences to be sampled.
#' @param sufix output file name be like <faFiles><_Ns>.fa, if NULL, then sufix is like '_5000s'.
#' @param ofiles must the same length as faFiles. If provided, then will ignore sufix.
#' @param outputDir if NULL, then output each file in the dir of corresponding `faFiles`.
#' @param addLabel if NULL, then add <addLabel> to seq title, like <title>:1.
#' @param seed random seed, if NULL then is the system time.
#' @return output file names.
#' @examples
#' \dontrun{
#' sampleNseqs(files, N=5000, outputDir='select5000sdir')
#' }
#' @export
sampleNseqs<-function(faFiles, N, sufix=NULL, ofiles=NULL,
                     outputDir=NULL, addLabel=NULL, seed=NULL) {

  if (is.null(sufix)) sufix=paste0('_',N,'s')

  if (is.null(ofiles)) {
    ofiles=gsub('.fa', '', faFiles)
    ofiles=paste0(ofiles, sufix, '.fa')
  } else {
    if (length(ofiles)!=length(faFiles))
      stop("ofiles should be the same length as faFiles!")
  }

  if (!is.null(outputDir)) {
    if (!dir.exists(outputDir)) dir.create(outputDir)
    ofiles=basename(ofiles)
    ofiles=paste0(formatPathFileName(outputDir, addBar=TRUE), ofiles)
  }

  for (f in faFiles) {
    cat('randomly select seqs', f, '\n')
    fa=seqinr::read.fasta(f, forceDNAtolower = F, as.string = TRUE)

    if (is.null(seed)) seed=Sys.time()
    set.seed(seed)
    idx1=sample(1:length(fa), size=N, replace = FALSE)

    fa=fa[idx1]

    if (!is.null(addLabel)) names(fa)=paste0(names(fa), addLabel)
    seqinr::write.fasta(sequences = fa,
                        names = names(fa),
                        nbchar = 1000,
                        file.out =ofiles[faFiles==f],
                        open="w")
  }
  return(ofiles)
}


## --------- DL-model data functions ----------

.formatChrs<-function(gr, chrs=NULL) {
  if (is.data.frame(gr)) {
    cname='seqnames'
    if (!(cname %in% colnames(gr))) cname='chr'
    if (!is.null(chrs)) {
      achrs=unique(gr[, cname])
      if (!any(achrs %in% chrs)) {
        if (grepl('Chr', chrs[1])) {
          cat(sprintf("Add [Chr] to chr names\n"))
          gr[, cname]=paste0('Chr', gr[, cname]) #chr1
        }else if (grepl('chr', chrs[1])) {
          cat(sprintf("Add [chr] to chr names\n"))
          gr[, cname]=paste0('chr', gr[, cname]) #Chr1
        }
        achrs=unique(gr[, cname])
        if (!any(achrs %in% chrs)) {
          cat( sprintf("Not any chr names in data.frame[%s] (%s) in chrs (%s)!\n", cname, toString(achrs), toString(chrs)))
        }
      }
      gr=gr[gr[, cname] %in% chrs, ]
    }
    return(gr)
  }

  if (!is.null(chrs)) {
    achrs=unique(gr@seqnames)
    if (!any(achrs %in% chrs)) {
      if (grepl('Chr', chrs[1])) {
        cat(sprintf("Add [Chr] to chr names\n"))
        gr@seqnames=paste0('Chr', gr@seqnames) #chr1
      }else if (grepl('chr', chrs[1])) {
        cat(sprintf("Add [chr] to chr names\n"))
        gr@seqnames=paste0('chr', gr@seqnames) #Chr1
      }
      achrs=unique(gr@seqnames)
      if (!any(achrs %in% chrs)) {
        cat( sprintf("Not any chr names in GRanges (%s) in chrs (%s)!\n", toString(achrs), toString(chrs)))
      }
    }
    gr=gr[gr@seqnames %in% chrs]
  }

  return(gr)
}



#' readBedFiles2GR reads BED files of pA sites and removes duplicates.
#' @details
#' This function will check the chr name consistency, and plot base composition for each BED file if plotFA=TRUE.
#' @param bedfiles one or more PA/BED files without header, with >=4 columns or =3 columns (chr/strand/start).
#' @param chrs a string vector storing chrs to subset, default is NULL.
#' @param plotFA TRUE to plot base composition for each BED file.
#' `N` PAs will be sampled plot each bed file, and output into BED_files_base_profile.pdf.
#' @param N number of PAs sampled to plot each BED file.
#' @param bsgenome if plotFA, bsgenome should be provided.
#' @param outputDir output dir for .freq, .fa and .pdf files (if plotFA=TRUE). NULL to output to the path of bedfiles.
#' @return GRanges storing pAs.
#' @export
readBedFiles2GR<-function(bedfiles, chrs=NULL,
                          plotFA=TRUE, N=5000,
                          bsgenome=NULL,
                          outputDir=NULL) {

  if (!is.null(outputDir)) outputDir=formatPathFileName(outputDir, addBar = TRUE)
  polyAsites=NULL
  fafiles=c()
  for (bedfile in bedfiles) {

    if (!file.exists(bedfile)) stop(sprintf("BED file %s not exists!\n", bedfile))
    polyAsite=read.table(bedfile, header = FALSE)

    nc=ncol(polyAsite)
    if (nc>=6) { #.bed polyAsite.bed
      cat("Ncol>=6, use the 1/2/3/6 columns in", bedfile, '\n')
      polyAsite=polyAsite[, c(1,2,3,6)]
    #} else if (nc>6) {
    #  cat("Ncol>6, use the firt 4 columns in", bedfile, '\n')
    #  polyAsite=polyAsite[, 1:4]
    } else if (nc==5) {
      cat("Ncol=5, use the 1/2/3/5 columns in", bedfile, '\n')
      polyAsite=polyAsite[, c(1,2,3,5)]
    } else if (nc==3) {
      colnames(polyAsite)=c('chr','strand','start')
      polyAsite$end=polyAsite$start
      polyAsite=polyAsite[, c('chr','start', 'end', 'strand')]
    } else if (nc==4) {
      colnames(polyAsite)=c('chr','start','end','strand')
      polyAsite=polyAsite[, c('chr','start', 'end', 'strand')]
    } else {
      stop("Number of columns in ", bedfile, "is ", nc,", should be >=4 or =3 (chr/strand/start)!")
    }
    colnames(polyAsite)=c('chr','start','end','strand')

    # remove dups
    polyAsite=unique(polyAsite)

    if (!is.null(chrs))
      polyAsite=.formatChrs(polyAsite, chrs=chrs)

    cat(sprintf("Read BED file %s: %d PAs\n", bedfile, nrow(polyAsite)))

    # plotFA
    if (plotFA) {
      ps=polyAsite[sample(1:nrow(polyAsite), size=min(N, nrow(polyAsite))), c('chr','strand','start')]
      colnames(ps)=c('chr','strand','coord')
      ps=movAPA::readPACds(ps)
      fapre=bedfile
      if (!is.null(outputDir)) {
        fapre=formatPathFileName(bedfile, noExt = TRUE, onlyFileName = TRUE)
        fapre=paste0(outputDir, fapre)
      }
      fafiles=c(fafiles, movAPA::faFromPACds(ps, bsgenome, what='updn', fapre=fapre))
    }

    if (is.null(polyAsites)) {
      polyAsites=polyAsite
    }  else  {
      polyAsites=rbind(polyAsites, polyAsite)
    }
  }

  polyAsites=with(polyAsites, GRanges(seqnames=chr, ranges=IRanges(start, end), strand=strand))
  if (length(bedfiles)>1)
    cat(sprintf("Read BED file total: %d PAs\n", length(polyAsites)))

  if (plotFA) {
    names(fafiles)=formatPathFileName(bedfiles, noExt=T, onlyFileName = TRUE)
    #filepre='BED_files_base_profile'
    movAPA::plotATCGforFAfile(faFiles=fafiles, ofreq=TRUE, opdf=T, refPos=301, mergePlots=FALSE)
  }
  return(polyAsites)
}


# expandGRPA expands up and down region of pA ranges.
# gr: GRanges
# extPA: expand upstream and downstream PA.
expandGRPA<-function(gr, extPA=200) {
  grExt=GenomicRanges::resize(gr, width=width(gr)+extPA, fix="end")
  grExt=GenomicRanges::resize(grExt, width=width(grExt)+extPA, fix="start")
  return(grExt)
}


#' getTxdbUTRs gets 3UTR ranges in TXDB.
#' @details
#' This function will output a list of utrsPA and utrsExt, which are used as reference regions for geting IP and real pAs.
#' @param txdb a TXDB object, or a data.frame with at least columns: seqnames, start, end ,strand
#' @param extLen length to extend 3UTR downstream for searching A-rich positions
#' @param plotFA TRUE to plot base composition for the 3UTR end position.
#' `N` PAs will be sampled plot each bed file, and output into BED_files_base_profile.pdf.
#' @param N number of PAs sampled to plot each BED file.
#' @param bsgenome if plotFA, bsgenome should be provided.
#' @param chrs a string vector storing chrs to subset, default is NULL.
#' If provided, will try to format the seqnames in txdb to make consistent with chrs.
#' @param outputDir output dir for .freq, .fa and .pdf files (if plotFA=TRUE). NULL to output to the path of bedfiles.
#' @return a list of utrsPA and utrsExt.
#' @export
getTxdbUTRs<-function(txdb, extLen=5000,
                      plotFA=TRUE,
                      N=5000, bsgenome=NULL, chrs=NULL,
                      outputDir=NULL) {
  if (!is.null(outputDir)) outputDir=formatPathFileName(outputDir, addBar = TRUE)

  if (class(txdb)=='data.frame') {
    if (!all(c('seqnames','start','end','strand') %in% colnames(txdb)))
      stop("getTxdbUTRs error: txdb is a data.frame, but seqnames,start,end,strand not all in txdb!\n")
    utrs=txdb[, c('seqnames','start','end','strand')]
  } else {
    utrs=GenomicFeatures::threeUTRsByTranscript(txdb)
    utrs=GenomicRanges::reduce(utrs, ignore.strand=FALSE)
    utrs=as.data.frame(utrs)
    #head(utrs)
    utrs=utrs[, c('seqnames','start','end','strand')]
  }

  utrs=unique(utrs)
  # format chrs
  utrs=.formatChrs(utrs, chrs=chrs)

  ## validate the 3UTR end, to see if they are like pAs
  if(plotFA & !is.null(bsgenome)) {
    ps=utrs[utrs$strand=='+', c('seqnames','strand','end')]
    colnames(ps)=c('chr','strand','coord')
    ps=ps[sample(1:nrow(ps), size=min(N, nrow(ps))), ]
    ps=movAPA::readPACds(ps)

    fapre='TXDB_3UTR'
    if (!is.null(outputDir)) {
      fapre=paste0(outputDir, 'TXDB_3UTR')
    }
    fafiles=movAPA::faFromPACds(ps, bsgenome, what='updn', fapre=fapre)

    filepre='TXDB_3UTR_base_profile'
    movAPA::plotATCGforFAfile (faFiles=fafiles, ofreq=TRUE, opdf=T, refPos=301, filepre=filepre)
  }

  ## extend txdb's 3UTR for 5000bp
  utrs=GenomicRanges::GRanges(utrs)
  utrsExt=GenomicRanges::resize(utrs, width=width(utrs)+extLen, fix="start")
  utrsExt=GenomicRanges::reduce(utrsExt)

  utrsPA=GenomicRanges::resize(utrs, width=1, fix="end")

  cat(sprintf("utrsPA: %d\n", length(utrsPA)))
  cat(sprintf("utrsExt (reduced): %d\n", length(utrsExt)))

  return(list(utrsPA=utrsPA, utrsExt=utrsExt))
}


#' getTxdbIntrons gets intron ranges in TXDB.
#' @param txdb a TXDB object , or a data.frame with at least columns: seqnames, start, end ,strand, to get intron ranges.
#' @param chrs a string vector storing chrs to subset, default is NULL.
#' If provided, will try to format the seqnames in txdb to make consistent with chrs.
#' @return a GRanges of introns
#' @export
getTxdbIntrons<-function(txdb, chrs=NULL) {

  if (class(txdb)=='data.frame') {
    if (!all(c('seqnames','start','end','strand') %in% colnames(txdb)))
      stop("getTxdbIntrons error: txdb is a data.frame, but seqnames,start,end,strand not all in txdb!\n")
    ranges=txdb[, c('seqnames','start','end','strand')]
  } else {
    ranges=GenomicFeatures::intronsByTranscript(txdb)
    ranges=GenomicRanges::reduce(ranges, ignore.strand=FALSE)
    ranges=as.data.frame(ranges)
    ranges=ranges[, c('seqnames','start','end','strand')]
  }

  ranges=unique(ranges)
  ranges=GenomicRanges::GRanges(ranges)

  cat(sprintf("introns#: %d\n", length(ranges)))

  if (!is.null(chrs)) {
    ranges=.formatChrs(ranges, chrs=chrs)
    cat(sprintf("introns in chrs#: %d\n", length(ranges)))
  }

  return(ranges)
}

#' getTxdb5UTRs gets 5UTR ranges in TXDB.
#' @param txdb a TXDB object , or a data.frame with at least columns: seqnames, start, end ,strand, to get 5UTR ranges.
#' @param chrs a string vector storing chrs to subset, default is NULL.
#' If provided, will try to format the seqnames in txdb to make consistent with chrs.
#' @return a GRanges of 5UTRs.
#' @export
getTxdb5UTRs<-function(txdb, chrs=NULL) {

  if (class(txdb)=='data.frame') {
    if (!all(c('seqnames','start','end','strand') %in% colnames(txdb)))
      stop("getTxdb5UTRs error: txdb is a data.frame, but seqnames,start,end,strand not all in txdb!\n")
    ranges=txdb[, c('seqnames','start','end','strand')]
  } else {
    ranges=GenomicFeatures::fiveUTRsByTranscript(txdb)
    ranges=GenomicRanges::reduce(ranges, ignore.strand=FALSE)
    ranges=as.data.frame(ranges)
    ranges=ranges[, c('seqnames','start','end','strand')]
  }

  ranges=unique(ranges)

  ranges=GenomicRanges::GRanges(ranges)

  cat(sprintf("5UTRs#: %d\n", length(ranges)))

  if (!is.null(chrs)) {
    ranges=.formatChrs(ranges, chrs=chrs)
    cat(sprintf("5UTRs in chrs#: %d\n", length(ranges)))
  }

  return(ranges)
}


#' getTxdbCDSs gets CDS ranges in TXDB.
#' @param txdb a TXDB object , or a data.frame with at least columns: seqnames, start, end ,strand, to get CDS ranges.
#' @param chrs a string vector storing chrs to subset, default is NULL.
#' If provided, will try to format the seqnames in txdb to make consistent with chrs.
#' @return a GRanges of CDS.
#' @examples
#' \dontrun{
#' getTxdbCDSs(txdb)
#' }
#' @export
getTxdbCDSs<-function(txdb, chrs=NULL) {

  if (class(txdb)=='data.frame') {
    if (!all(c('seqnames','start','end','strand') %in% colnames(txdb)))
      stop("getTxdbCDSs error: txdb is a data.frame, but seqnames,start,end,strand not all in txdb!\n")
    ranges=txdb[, c('seqnames','start','end','strand')]
  } else {
    ranges=GenomicFeatures::cdsBy(txdb,'gene')
    ranges=GenomicRanges::reduce(ranges, ignore.strand=FALSE)
    ranges=as.data.frame(ranges)
    ranges=ranges[, c('seqnames','start','end','strand')]
  }

  ranges=unique(ranges)

  ranges=GenomicRanges::GRanges(ranges)

  cat(sprintf("CDS#: %d\n", length(ranges)))

  if (!is.null(chrs)) {
    ranges=.formatChrs(ranges, chrs=chrs)
    cat(sprintf("CDS in chrs#: %d\n", length(ranges)))
  }

  return(ranges)
}


#' getTxdbIntergenics gets intergenic ranges in TXDB.
#' @param txdb a TXDB object , or a data.frame with at least columns: seqnames, start, end ,strand, to get intergenic ranges.
#' @param chrs a string vector storing chrs to subset, default is NULL.
#' If provided, will try to format the seqnames in txdb to make consistent with chrs.
#' @param extLen length to extend transcript downstream to narrow intergenic ranges.
#' @return a GRanges of intergenics, whose both strands not within any transcripts+extLen.
#' @examples
#' \dontrun{
#' itgs=getTxdbIntergenics(txdb)
#' }
#' @export
getTxdbIntergenics<-function(txdb, chrs=NULL, extLen=5000) {

  if (class(txdb)=='data.frame') {
    if (!all(c('seqnames','start','end','strand') %in% colnames(txdb)))
      stop("getTxdbIntergenics error: txdb is a data.frame, but seqnames,start,end,strand not all in txdb!\n")
    ranges=txdb[, c('seqnames','start','end','strand')]
  } else {

    ranges=GenomicFeatures::transcripts(txdb)
    ranges=GenomicRanges::reduce(ranges, ignore.strand=TRUE)

    ranges=ranges[width(ranges)>2*extLen]

    ranges=IRanges::resize(ranges, width=IRanges::width(ranges)-extLen, fix='start')
    ranges=IRanges::resize(ranges, width=IRanges::width(ranges)-extLen, fix='end')

    ranges=as.data.frame(ranges)
    ranges=ranges[, c('seqnames','start','end','strand')]
    ranges=unique(ranges)
    ranges$strand='+'
    ranges2=ranges
    ranges2$strand='-'
    ranges=rbind(ranges, ranges2)
  }

  ranges=GenomicRanges::GRanges(ranges)

  cat(sprintf("intergenics#: %d\n", length(ranges)))

  if (!is.null(chrs)) {
    ranges=.formatChrs(ranges, chrs=chrs)
    cat(sprintf("intergenics in chrs#: %d\n", length(ranges)))
  }

  return(ranges)
}


#' getFreePARangesIn3UTR gets genomic regions including 3'UTRs and downstream regions but free of pAs.
#' @details
#' The output GRanges contains 3UTRs and their downstream 5000bp regions, but withou any pAs.
#' These ranges do not include potential PA or near the end of 3UTR [-200,+200], which are completely PA-free regions.
#' Each region in the GRanges is at least 500 nt.
#' @param paBedFiles one or more BED files.
#' @param txdb a TXDB object or a data.frame with columns (seqnames, strand, start, end) to get 3UTR ranges.
#' @param extUTRLen length to extend 3UTR downstream for searching A-rich positions.
#' Larger extUTRlen to get more utrsNoPA for subsequent IP scanning.
#' @param extPA extend PA regions.
#' @param outputRds the output .RDS file name, default is `utrsNoPA.RDS`.
#' @param chrs a string vector storing chrs to subset, default is NULL.
#' If provided, will try to format the seqnames in txdb to make consistent with chrs.
#' @param plotFA TRUE to plot base composition for the 3UTR end position.
#' `N` PAs will be sampled plot each bed file, and output into BED_files_base_profile.pdf.
#' @param N number of PAs sampled to plot each BED file.
#' @param bsgenome if plotFA, bsgenome should be provided.
#' @param outputDir the output dir for .RDS file and other .freq/.pdf files (if plotFA=TRUE)
#' @return a GRanges object, and output the .RDS file to `outputRds`.
#' @export
getFreePARangesIn3UTR<-function(paBedFiles, txdb, extUTRLen=5000, extPA=200,
                          outputRds='utrsNoPA.RDS', chrs=NULL,
                          plotFA=TRUE, N=5000, bsgenome=NULL,
                          outputDir='./') {

  outputDir=formatPathFileName(outputDir, addBar = TRUE)
  ## read polyA data
  cat('###### readBedFiles2GR...\n')
  polyAsite=readBedFiles2GR(paBedFiles, chrs=chrs,
                            plotFA=plotFA, N=N, bsgenome=bsgenome,
                            outputDir=outputDir)

  cat('###### expandGRPA...\n')
  polyAsite=expandGRPA(polyAsite, extPA=extPA)

  cat('###### getTxdbUTRs...\n')
  # Extract the 3UTR region of TXDB as PA and utrExt
  utrs=getTxdbUTRs(txdb, extLen=extUTRLen,
                   plotFA=plotFA, N=N, bsgenome=bsgenome,
                   chrs=chrs,
                   outputDir=outputDir)
  utrsPA=utrs$utrsPA
  utrsExt=utrs$utrsExt
  utrs=NULL
  utrsPA=expandGRPA(utrsPA, extPA=extPA)

  ## Obtain the 3UTR region without PA
  ## Merge polyAdb and polyAsite, and use the extended area of UTRs as the PA regions
  ## these are potential PA containing regions that need to be excluded
  cat('###### reduce PAs+UTRPAs...\n')
  polyA=c(polyAsite, utrsPA)
  polyA=GenomicRanges::reduce(polyA)

  ## Remove the polyA region from UTRSExt to obtain the remaining region that definitely does not contain PA:
  ## first disjoin UTRs, and then remove the region that overlaps with polyA
  cat('###### disjoin PAs+UTRPAs+UTRExts...\n')
  all=c(utrsExt, polyA)
  all=GenomicRanges::disjoin(all) # cut utrs by pA ranges
  utrsNoPA=subsetByOverlaps(all, polyA, minoverlap=1, type='any', invert=TRUE)

  n=GenomicRanges::countOverlaps(utrsNoPA, polyA) #Verification, there is indeed no OVP
  if (sum(n)>0) cat("Warning: countOverlaps(utrsNoPA, polyA)>0!\n")

  cat('###### output final PA-free UTRs >500nt...\n')
  ##Take at least 500nt of the regions as a candidate for scanning final nonPAs
  utrsNoPA=utrsNoPA[width(utrsNoPA)>=500]
  cat("utrsNoPA:", length(utrsNoPA),'\n')

  outputRds=formatPathFileName(outputRds,onlyFileName = TRUE)
  saveRDS(utrsNoPA, file=paste0(outputDir, outputRds))
  return(utrsNoPA)
}


## .getFreePARangesInIntron gets genomic regions of introns but free of pAs.
## @details
## The output GRanges contains only introns withou any pAs.
## Each region in the GRanges is at least 500 nt.
## @param paBedFiles one or more BED files.
## @param txdb a TXDB object or a data.frame with columns (seqnames, strand, start, end) to get intron ranges.
## @param outputRds the output .RDS file name, default is `intronsNoPA.RDS`.
## @param chrs a string vector storing chrs to subset, default is NULL.
## If provided, will try to format the seqnames in txdb to make consistent with chrs.
## @return a GRanges object, and output the .RDS file to `outputRds`.
#' @export
.getFreePARangesInIntron<-function(paBedFiles, txdb,
                                outputRds='intronsNoPA.RDS', chrs=NULL) {

  cat('###### readBedFiles2GR...\n')
  polyAsite=readBedFiles2GR(paBedFiles, chrs=chrs, plotFA=FALSE)

  cat('###### getTxdbIntrons...\n')
  ranges=getTxdbIntrons(txdb)
  ranges=.formatChrs(ranges, chrs=chrs)

  cat('###### get PA-free introns...\n')
  rangesNoPA=subsetByOverlaps(ranges, polyAsite, minoverlap=1, type='any', invert=TRUE)

  n=GenomicRanges::countOverlaps(rangesNoPA, polyAsite) # validate, indeed no ovp
  if (sum(n)>0) cat("Warning: countOverlaps(introns, polyAsite)>0!\n")

  cat('###### output final PA-free introns >500nt...\n')
  rangesNoPA=rangesNoPA[width(rangesNoPA)>=500]
  cat("intronsNoPA:", length(rangesNoPA),'\n')

  saveRDS(rangesNoPA, file=outputRds)
  return(rangesNoPA)
}


#' getFreePARangesInNon3UTRs gets PA-free regions of 5UTR/CDS/intron/intergenic.
#' @details
#' The output GRanges contains only ranges without any pAs.
#' Each region in the GRanges is at least 500 nt.
#' For long ranges, will tile into `tileWidth` bp ranges to split long regions.
#' @param paBedFiles one or more BED files.
#' @param txdb  a TXDB object.
#' @param region genomic region to output, including all, cds, intron, 5utr, intergenic.
#' @param outputRds the output .RDS file name, default is `non3UTRsNoPA.RDS`.
#' If region=all, then the output RDS will be like <non3UTRsNoPA_5utr/intron/intergneic/cds.RDS>.
#' @param chrs a string vector storing chrs to subset, default is NULL.
#' If provided, will try to format the seqnames in txdb to make consistent with chrs.
#' @param tileWidth tile larger region to smaller pieces.
#' @return a GRanges object, and output the .RDS file(s) to `outputRds`.
#' @export
getFreePARangesInNon3UTRs<-function(paBedFiles, txdb, region='all',
                                    chrs=NULL, tileWidth=1000,
                                    outputRds='non3UTRsNoPA.RDS') {

  cat('###### readBedFiles2GR...\n')
  polyAsite=readBedFiles2GR(paBedFiles, chrs=chrs, plotFA=FALSE)

  region=tolower(region)
  if (region=='all') {
    region=c('intron','cds','5utr','intergenic')
    outputRds=gsub('.RDS','',outputRds, fixed=TRUE)
    outputRds=paste0(outputRds, '_', region, '.RDS')
  }

  if (!all(region %in% c('intron','cds','5utr','intergenic')))
    stop("region should be in all/intron/cds/5utr/intergenic!")

  for (rg in region) {
    cat('###### getTxdb', rg, '...\n')
    if (rg=='intron') {
      ranges=getTxdbIntrons(txdb)
    }
    if (rg=='5utr') {
      ranges=getTxdb5UTRs(txdb)
    }
    if (rg=='intergenic') {
      ranges=getTxdbIntergenics(txdb)
    }
    if (rg=='cds') {
      ranges=getTxdbCDSs(txdb)
    }

    # for large regions, split to multiple 10000bp ranges, to avoid long ranges to be removed only because containing one pA
    # but only have effects on intron/intergenics because cds and 5utr are short.
    ranges=unlist(GenomicRanges::tile(ranges, width=tileWidth))

    ranges=.formatChrs(ranges, chrs=chrs)

    cat('###### get PA-free', rg, '\n')
    rangesNoPA=subsetByOverlaps(ranges, polyAsite, minoverlap=1, type='any', invert=TRUE)

    n=GenomicRanges::countOverlaps(rangesNoPA, polyAsite)
    if (sum(n)>0) cat("Warning: countOverlaps(region, polyAsite)>0!\n")

    cat('###### output final PA-free', rg , '(>500nt)...\n')
    rangesNoPA=rangesNoPA[width(rangesNoPA)>=500]
    cat("PA-free:", rg, ':', length(rangesNoPA),'\n')

    saveRDS(rangesNoPA, file=outputRds[region==rg])
  }
  return(outputRds)
}



## ----- A-rich PA seqs ------

#' getArichPAseqs gets A-rich pA sequences within `granges`.
#' @details
#' The function first gets reference pAs from `paBedFiles`, subsetting by `chrs`. And then extract A-rich pAs.
#' @param paBedFiles one or more BED files storing PAs.
#' @param granges a GRanges object to subset pAs.
#' @param bsgenome if plotFA, bsgenome should be provided.
#' @param chrs a string vector storing chrs to subset, default is NULL.
#' @param grams a string vector of grams (pA signal) to search within -10 to -50 of pAs.
#' pA sequences containing each gram are saved to one .fa file.
#' If grams=NULL, then polyA signal is not scanned, just output all A-rich seqs.
#' @param gramsPriority same length as grams, to set the priority of searching, e.g., c(1, 2, 3, 3, 3)
#' @param outputSeqPre the output .fa file is like <outputSeqPre>.AATAAA.fa.
#' @param shift  0 not to shift PA, >0 shift to 3' end, <0 shift to 5'end.
#' @return output fa file names, and output sequences of each gram to a .fa file, like <realPA_Arich_seq.AATAAA/NOPAS.fa>.
#' @examples
#' \dontrun{
#' getArichPAseqs(paBedFiles='Human_hg38.PolyA_DB.bed',
#'                granges, bsgenome=bsgenome, chrs=chrs,
#'                grams=gramsMM, gramsPriority=gramsPriority,
#'                outputSeqPre="realPA_Arich_seq")
#' }
#' @export
getArichPAseqs<-function(paBedFiles, granges,
                         bsgenome, chrs=NULL,
                         grams=NULL, gramsPriority=NULL,
                         outputSeqPre="realPA_Arich_seq",
                         shift=0) {

  outputSeqPre=formatPathFileName(outputSeqPre, addBar = FALSE)
  dname=dirname(outputSeqPre)
  if (!dir.exists(dname)) dir.create(dname)

  ## read polyA data
  cat('###### readBedFiles2GR...\n')
  polyAdb=readBedFiles2GR(paBedFiles, chrs=chrs, plotFA=FALSE)

  cat('###### subset PAs in granges...\n')
  polyAdbinUTR=subsetByOverlaps(polyAdb, granges)

  pacds=as.data.frame(polyAdbinUTR)
  colnames(pacds)=c('chr','UPA_start','UPA_end','width','strand')

  if (shift!=0) {
    cat(sprintf('###### shift PA position (%d) nt [>0 right, <0 left]...\n', shift))
    pacds$coord=pacds$UPA_start+shift
    pacds$coord[pacds$strand=='-']=pacds$UPA_start[pacds$strand=='-']-shift
  } else {
    pacds$coord=pacds$UPA_start
  }

  pacds=movAPA::readPACds(pacds)

  cat('###### get Arich PAs by removePACdsIP...\n')
  pacdsIP=movAPA::removePACdsIP(pacds, bsgenome, returnBoth=TRUE, up=-10, dn=10, conA=6, sepA=7)
  pacdsIP=pacdsIP$ip #Real PA similar to IP

  # Find PAS upstream of PA
  cat('###### search upstream PA signals...\n')
  if (!is.null(grams)) {
    pacdsIP=movAPA::annotateByPAS(pacdsIP, bsgenome, grams=grams, from=-50, to=-10, priority=gramsPriority, label='PAS')

    # only 1 gram in grams
    if (!('PAS_gram' %in% colnames(pacdsIP@anno))) {
      pacdsIP@anno$PAS_gram=grams
      pacdsIP@anno$PAS_gram[is.na(pacdsIP@anno$PAS_dist)]=NA
    }

    pacdsIP@anno$PAS_gram[is.na(pacdsIP@anno$PAS_gram)]='NOPAS'
    print(sort(table(pacdsIP@anno$PAS_gram), decreasing = TRUE))

  } else {
    pacdsIP@anno$PAS_gram='NOPAS'
  }


  cat('###### Output 200nt seq (PA is 101nt) for each gram (or all to NOPAS if gram=NULL) ...\n')
  # Output a PA sequence containing the given gram
  fafiles=movAPA::faFromPACds(pacdsIP, bsgenome, what='updn', fapre=outputSeqPre, up=-100, dn=99, byGrp='PAS_gram')

  if (!is.null(grams)) {
    ng=min(5, length(grams))
    movAPA::plotATCGforFAfile(faFiles=paste0(outputSeqPre, '.', c(grams[1:ng], 'NOPAS'), '.fa'), ofreq=FALSE, opdf=TRUE,
                      filepre=paste0(outputSeqPre,'_',ng+1,'grams_base_profile'), mergePlots=TRUE)
  } else {
    movAPA::plotATCGforFAfile(faFiles=paste0(outputSeqPre, '.NOPAS.fa'), ofreq=FALSE, opdf=TRUE,
                              filepre=paste0(outputSeqPre, '.NOPAS.base_profile'), mergePlots=TRUE)
  }

  return(fafiles)
}



#' getArichPAseqsIn3UTR gets A-rich pA sequences in 3'UTRs.
#' @details
#' The function first gets reference pAs from `paBedFiles`, subsetting by `chrs`. And then extract A-rich pAs in 3UTRs.
#' @param paBedFiles one or more BED files.
#' @param txdb  a TXDB object or a data.frame with columns (seqnames, strand, start, end) to get 3UTR ranges.
#' @param extUTRLen length to extend 3UTR downstream for searching A-rich positions
#' @param bsgenome if plotFA, bsgenome should be provided.
#' @param chrs a string vector storing chrs to subset, default is NULL.
#' If provided, will try to format the seqnames in txdb to make consistent with chrs.
#' @param grams a string vector of grams (pA signal) to search. pA sequences containing each gram are saved to one .fa file.
#' @param gramsPriority same length as grams, to set the priority of searching, e.g., c(1, 2, 3, 3, 3)
#' @param outputSeqPre the output .fa file is like <outputSeqPre>.AATAAA.fa. Default is realPA_Arich_seq.
#' @param shift  0 not to shift PA, >0 shift to 3' end, <0 shift to 5'end.
#' @return output fa file names, and output sequences of each gram to a .fa file.
#' @export
getArichPAseqsIn3UTR<-function(paBedFiles, txdb, extUTRLen=5000,
                               bsgenome, chrs=NULL,
                               grams=NULL, gramsPriority=NULL,
                               outputSeqPre="realPA_Arich_seq",
                               shift=0) {

  cat('###### getTxdbUTRs...\n')
  # Extract the 3UTR region of TXDB as PA and UTrext
  utrs=getTxdbUTRs(txdb, extLen=extUTRLen, plotFA=FALSE, chrs=chrs)
  utrsExt=utrs$utrsExt
  utrs=NULL

  fafiles=getArichPAseqs(paBedFiles, granges=utrsExt,
                         bsgenome=bsgenome, chrs=chrs,
                         grams=grams, gramsPriority=gramsPriority,
                         outputSeqPre=outputSeqPre,
                         shift=shift)

  return(fafiles)
}


#' getArichPAseqsInIntron gets A-rich pA sequences in introns.
#' @details
#' The function first gets reference pAs from `paBedFiles`, subsetting by `chrs`. And then extract A-rich pAs in introns.
#' @param paBedFiles one or more BED files.
#' @param txdb  a TXDB object or a data.frame with columns (seqnames, strand, start, end) to get intron ranges.
#' @param bsgenome if plotFA, bsgenome should be provided.
#' @param chrs a string vector storing chrs to subset, default is NULL.
#' If provided, will try to format the seqnames in txdb to make consistent with chrs.
#' @param grams a string vector of grams (pA signal) to search. pA sequences containing each gram are saved to one .fa file.
#' @param gramsPriority same length as grams, to set the priority of searching, e.g., c(1, 2, 3, 3, 3)
#' @param outputSeqPre the output .fa file is like <outputSeqPre>.AATAAA.fa. Default is intron_realPA_Arich_seq.
#' @param shift  0 not to shift PA, >0 shift to 3' end, <0 shift to 5'end.
#' @return output fa file names, and output sequences of each gram to a .fa file.
#' @export
getArichPAseqsInIntron<-function(paBedFiles, txdb, bsgenome, chrs=NULL,
                               grams=NULL, gramsPriority=NULL,
                               outputSeqPre="intron_realPA_Arich_seq",
                               shift=0) {

  cat('###### getTxdbIntrons...\n')
  ranges=getTxdbIntrons(txdb, chrs=chrs)

  fafiles=getArichPAseqs(paBedFiles, granges=ranges,
                         bsgenome=bsgenome, chrs=chrs,
                         grams=grams, gramsPriority=gramsPriority,
                         outputSeqPre=outputSeqPre,
                         shift=shift)

  return(fafiles)
}



## ----- A-rich IP seqs -----

#' getArichIPseqs gets A-rich internal priming (IP) sequences based on given PA-free regions.
#' @details
#' Given regionsNoPARDS (from getFreePARangesIn.. functions) with no PA at all, this function outputs the IP sequence of A-rich.
#' It takes about 30m to output about 9w sequences (500 sequences/10s) on a normal notebook.
#' @param regionsNoPARDS a file from from getFreePARangesIn..function, which stores regions without any pAs.
#' @param bsgenome if plotFA, bsgenome should be provided.
#' @param grams a string vector of grams (pA signal) to search. pA sequences containing each gram are saved to one .fa file.
#' @param N the number of sequences needed. If NULL, output all.
#' @param outputSeqPre the output .fa file is like <outputSeqPre>.AATAAA.fa. Default is IP_Arich_seq.
#' @return output fa file names, and output sequences of each gram (or NOPAS) to a .fa file.
#' @examples
#' \dontrun{
#' getArichIPseqs(utrsNoPARDS='utrsNoPA.RDS', bsgenome, grams=gramsMM,
#'                  N=NULL, outputSeqPre="IP_Arich_seq")
#' getArichIPseqs(utrsNoPARDS='utrsNoPA.RDS', bsgenome, grams=NULL,
#'                  N=2000, outputSeqPre="IP_Arich_seq")
#' }
#' @export
getArichIPseqs<-function(regionsNoPARDS, bsgenome,
                         grams=NULL,
                         N=NULL,
                         outputSeqPre="IP_Arich_seq") {

  outputSeqPre=formatPathFileName(outputSeqPre, addBar = FALSE)
  dname=dirname(outputSeqPre)
  if (!dir.exists(dname)) dir.create(dname)

  utrsNoPA=readRDS(regionsNoPARDS)
  cat('###### load regionsNoPARDS:', length(utrsNoPA),'\n')

  ## First, build PACDS and extract the sequence
  pacds=as.data.frame(utrsNoPA)
  colnames(pacds)=c('chr','UPA_start','UPA_end','width','strand')
  pacds$coord=pacds$UPA_start
  pacds=movAPA::readPACds(pacds)

  cat('###### get regionsNoPA seqs...\n')
  seq=movAPA::faFromPACds(pacds, bsgenome, what='pac', fapre=NULL, up=NULL, dn=NULL)

  if (!is.null(grams)) {
    cat('###### scan PAS (grams)...\n')
    ## scan grams and An
    vms=list()

    for(g in grams) {
      vm=Biostrings::vmatchPattern(pattern=g, subject=seq, fixed=TRUE)
      vms[[g]]=vm
    }
  } else {
    vms=NULL
  }

  cat('###### scan An...\n')
  invisible(gc())
  vmIP=Biostrings::vmatchPattern(pattern='AAAAAA', subject=seq, fixed=TRUE)

  # vms[[1]] -- noPA region's AATAAA positions
  # MIndex object of length 92545
  # $`PA1:chr1;+;70009;71384`
  # IRanges object with 2 ranges and 0 metadata columns:
  #   start       end     width
  # <integer> <integer> <integer>
  #   [1]       927       932         6
  # [2]      1236      1241         6

  # > vmIP  -- noPA region's AAAAAA positions
  # MIndex object of length 92545
  # $`PA1:chr1;+;70009;71384`
  # IRanges object with 15 ranges and 0 metadata columns:
  #   start       end     width
  # <integer> <integer> <integer>
  #   [1]       344       349         6
  # [2]       345       350         6
  # [3]       346       351         6


  #lapply(vms, length) #all nseq

  NCHAR=100 # The length of the upstream/downstream sequence to be extracted
  gramSeqs=list() # Store the output seq for each type of gram
  #The final truncated sequence consists of the IP position and 100nt on each side, totaling 201nt.
  #The IP position starts with An, meaning An is downstream of the IP

  # There are over 9w candidate sequences in total, and the output of all of them is relatively slow
  # 5000 output IP sequences can be randomly selected here.
  # If not enough, all sequences can be used: sid=1:length(vmIP)
  if (is.null(N)) {
    sid=1:length(vmIP)
  } else  {
    sid=sample(1:length(vmIP), min(N, length(vmIP)), replace = FALSE)
  }
  cat(sprintf("Output N=%d from total %d seqs\n", N, length(vmIP)))

  #For each seq (noPA region), find the IP location and determine the PAS pairing for each IP
  counter=0
  for (i in sid) {
    counter=counter+1
    vip=vmIP[[i]]
    if (length(vip)==0) next
    seqlen=width(seq[i])

    if (length(vip)>1) vip=IRanges::reduce(vip) #combine nearing As

    # Ensure that there are>=100nt sequences on both sides of the IP
    vip=vip[start(vip)>NCHAR & end(vip)<seqlen-NCHAR, ]
    if (length(vip)==0) next

    # For these IPs, extract the upstream 10-50 range of the IP,
    # sequentially determine whether each gram is within this range,
    # and output the extracted seq
    pas=vip
    start(pas)=start(pas)-10
    end(pas)=start(pas)
    pas=IRanges::resize(pas, width=40, fix='end')

    if (counter %% 500==0) cat(counter, 'of',length(sid), 'done ...\n')

    noPAS=1


    .getIPseq<-function(id, ig) {
      ct=start(vip)[id][1]+floor(width(vip)[id][1]/2) #An's center
      IPpos=sample((ct-10):(ct+10), 1) #Randomly select a location within the range of 10nt upstream and downstream of the An center as the PA

      if (IPpos-NCHAR<1 | IPpos+NCHAR>Biostrings::nchar(seq[i])) return(c())
      oseq=substring(seq[i], IPpos-NCHAR, IPpos+NCHAR-1)
      otitle=paste0('>', names(seq)[i], ';IP=', IPpos, ';', ig, ':0') #add label
      return(c(otitle, oseq))
    }

    if (!is.null(grams)) {
      for (g in grams) { # For each type of gram, determine whether the position of the gram in the seq is within the IP ranges of that seq
        vg=vms[[g]][[i]]
        if (length(vg)==0) next
        ovp=IRanges::findOverlaps(vg, pas, type='within', select='first')

        if (sum(!is.na(ovp))==0) next

        ivg=which(!is.na(ovp))
        ipas=ovp[ivg] #The IP index corresponding to PAS is stored

        #Extract the upstream and downstream sequences corresponding to each iPad in VIP
        #Each original seq here has multiple IPs paired with PAS, and only the pairing result of the first IP is output,
        #so that the output sequences do not overlap

        #20230813 Set the fake PA position to a random position within the upstream and downstream range of An
        # (the original code was that the IP position started with An)

        ss=.getIPseq(id=ipas, ig=g)
        if (length(ss)>0) {
          gramSeqs[[g]]=c(gramSeqs[[g]], ss)
        } else {
          next
        }

        # ct=start(vip)[ipas][1]+floor(width(vip)[ipas][1]/2) #center of An
        # IPpos=sample((ct-10):(ct+10), 1) #Randomly select a location within the range of 10nt upstream and downstream of the An center as the PA
        #
        # if (IPpos-NCHAR<1 | IPpos+NCHAR>nchar(seq[i])) next
        # oseq=substring(seq[i], IPpos-NCHAR, IPpos+NCHAR-1)
        # otitle=paste0('>', names(seq)[i], ';IP=', IPpos, ';', g, ':0') #add label
        # gramSeqs[[g]]=c(gramSeqs[[g]], otitle, oseq)

        noPAS=0
      }
    }

    # No gram, output to NOPAS, select the first one from the VIP of this seq and output
    if (noPAS) {
      ss=.getIPseq(id=1, ig='NOPAS')
      if (length(ss)>0) {
        gramSeqs[['NOPAS']]=c(gramSeqs[['NOPAS']], ss)
      }
    }

  }

  ## Output the IP sequence for each type of gram
  for (g in names(gramSeqs)) {
    cat(g, 'nseq=', length(gramSeqs[[g]])/2, '\n')
    write.table(gramSeqs[[g]], file=paste0(outputSeqPre,'.', g, '.fa'), sep="\n", col.names = F, row.names = F, quote=F)
  }

  ## Create the ATCG distribution for the first 6 grams
  if (!is.null(grams)) {
    ng=min(5, length(grams))
    movAPA::plotATCGforFAfile(faFiles=paste0(outputSeqPre, '.', c(grams[1:ng], 'NOPAS'), '.fa'), ofreq=FALSE, opdf=TRUE,
                              filepre=paste0(outputSeqPre,'_',ng+1,'grams_base_profile'), mergePlots=TRUE)
  } else {
    movAPA::plotATCGforFAfile(faFiles=paste0(outputSeqPre, '.NOPAS.fa'), ofreq=FALSE, opdf=TRUE,
                              filepre=paste0(outputSeqPre, '.NOPAS.base_profile'), mergePlots=TRUE)
  }

}



## read filename like 'realPA_Arich_seq.ATTATA.fa', 'NOPAS.fa'
.readFas<-function(path, filePre, readNOPAS=TRUE, verbose=FALSE) {
  files=list.files(path, pattern=paste0(filePre,'.*fa$'))
  if (length(files)==0) {
    stop("No ", filePre," file in ", path, '\n')
  }
  fas=list()
  for (f in files) {
    gram=f
    if (filePre!='') gram=gsub(filePre, '', f, fixed=TRUE)
    gram=gsub('.fa', '', gram, fixed=TRUE)
    #if (nchar(gram)!=6) next
    if (!readNOPAS & gram=='NOPAS') next
    fa=seqinr::read.fasta(paste0(path,'/', f), forceDNAtolower = F)
    fas[[gram]]=fa
    if (verbose) cat(sprintf("%s %s %d seqs\n", filePre, gram, length(fa)))
  }
  return(fas)
}

##
#' statCntFas stats seq number in filename like 'realPA_Arich_seq.ATTATA.fa', 'NOPAS.fa'
#' @param path path of the fa files
#' @param filePre prefix of the fa file names, e.g., for file "realPA_Arich_seq.NOPAS.fa", the filepre is "realPA_Arich_seq.".
#' The gram should be present in the file name, that is the file name is <filePre>.gram.fa.
#' @param readNOPAS TRUE to read NOPAS file
#' @param justTotal FALSE to not count total number
#' @return NULL, output to screen
#' @export
statCntFas<-function(path, filePre='', readNOPAS=TRUE, justTotal=FALSE) {
  fas=.readFas(path, filePre, readNOPAS=readNOPAS, verbose=FALSE)
  fas=unlist(lapply(fas, length))
  fas=sort(fas, decreasing = TRUE)
  if (justTotal) {
    cat(path, sum(fas), '\n')
    return(invisible(NULL))
  } else {
    print(filePre)
    print(fas)
    cat('Total', sum(fas), '\n')
    return(fas)
  }
}

#' split2TrainTest randomly splits fa files (each gram one file) to two folders (e.g., train and test)
#' @param path path of the fa files
#' @param filePre prefix of the fa file names, e.g., for file "realPA_Arich_seq.NOPAS.fa", the filepre is "realPA_Arich_seq.".
#' The gram should be present in the file name, that is the file name is <filePre>.gram.fa.
#' @param dir1 like 'dir1' or '../xx/dir2'. split `per1` percentage of fa files to dir1. If dir is not exists, will be created.
#' @param dir2 split remaining (`1-per1`) percentage of fa files to dir2
#' @param per1 percentage of sequences to be extracted into dir1, default is 0.5
#' @param label add label to the sequence title, default is ":1".
#' @return NULL, output files to dir1 and dir2.
#' The gram (file name) and the `label` will be added to the sequence title, like <title>;AATAAA:1.
#' @examples
#' \dontrun{
#' split2TrainTest(path='polyAseqTrap/modelData/',
#'                 filePre='realPA_Arich_seq.',
#'                 dir1='realPA_Arich_seq_train', dir2='realPA_Arich_seq_test',
#'                 per1=0.7, label=":1")
#' }
#' @export
split2TrainTest<-function(path, filePre,
                          dir1='train', dir2='test',
                          per1=0.5, label=":1") {

  if (per1>=1) stop("per1 should be <1 (split to dir1 and dir2)!")

  if (dir.exists(dir1)) stop("dir1 exists!\n")
  if (dir.exists(dir2)) stop("dir2 exists!\n")
  dir.create(dir1); dir.create(dir2)

  cat(paste0("split2TrainTest: the new sequence title will like <old title>;AATAAA", label, "\n"))

  fas=.readFas(path, filePre)

  for (g in names(fas)) {
    fa=fas[[g]]
    n1=floor(per1*length(fa))
    idx1=sample(1:length(fa), size=n1, replace = FALSE)
    idx2=(1:length(fa))[-idx1]

    cat(sprintf("%s dir1=%d (%f); dir2=%d; total=%d\n", g, length(idx1), per1, length(idx2), length(fa)))
    fa1=fa[idx1] #AATAAA train real
    fa2=fa[idx2] #AATAAA test real

    # Add \ 1 to the end of true and add a gram tag
    # (previously output sequence did not have a gram tag, IP has a gram tag)
    seqinr::write.fasta(sequences = fa1,
                        names = paste0(names(fa1),';', g, label),
                        nbchar = 1000,
                        file.out = paste0(dir1, '/', g, '.fa') )

    seqinr::write.fasta(sequences = fa2,
                        names = paste0(names(fa2),';', g, label),
                        nbchar = 1000,
                        file.out = paste0(dir2, '/', g, '.fa') )
  }
}


#' split2TrainTestByFiles splits each of fafiles to dir1 (per1) and dir2 (1-per1), adding label to the title.
#' @param fafiles fa file names.
#' @param dir1 like 'dir1' or '../xx/dir2'. split `per1` percentage of fa files to dir1. If dir is not exists, will be created.
#' @param dir2 split remaining (`1-per1`) percentage of fa files to dir2
#' @param per1 percentage of sequences to be extracted into dir1, default is 0.5
#' @param label add label to the sequence title, default is ":1".
#' @return NULL, output files to dir1 and dir2.
#' The `label` will be added to the sequence title, like <title>:1.
#' @examples
#' \dontrun{
#' split2TrainTestByFiles(fafiles,
#'                       dir1='3utr_shiftNonDRSPA_seq_train_per70',
#'                       dir2='3utr_shiftNonDRSPA_seq_test_per30',
#'                       per1=0.7, label=":0")
#' }
#' @export
split2TrainTestByFiles<-function(fafiles,
                          dir1='train', dir2='test',
                          per1=0.5, label=":1") {

  if (per1>=1) stop("per1 should be <1 (split to dir1 and dir2)!")

  if (dir.exists(dir1)) stop("dir1 exists!\n")
  if (dir.exists(dir2)) stop("dir2 exists!\n")
  dir.create(dir1); dir.create(dir2)

  cat(paste0("split2TrainTest: the new sequence title will like <old title>", label, "\n"))

  for (f in fafiles) {

    fa=seqinr::read.fasta(f, forceDNAtolower = F, as.string = TRUE)

    n1=floor(per1*length(fa))
    idx1=sample(1:length(fa), size=n1, replace = FALSE)
    idx2=(1:length(fa))[-idx1]

    cat(sprintf("%s dir1=%d (%f); dir2=%d; total=%d\n", f, length(idx1), per1, length(idx2), length(fa)))
    fa1=fa[idx1] # train real
    fa2=fa[idx2] # test real

    seqinr::write.fasta(sequences = fa1,
                        names = paste0(names(fa1), label),
                        nbchar = 1000,
                        file.out = paste0(dir1, '/', basename(f)))

    seqinr::write.fasta(sequences = fa2,
                        names = paste0(names(fa2), label),
                        nbchar = 1000,
                        file.out = paste0(dir2, '/', basename(f)))
  }
  return(invisible(NULL))

}




## ----- train/test seq files -----


# this function is to sample sequences from a dir N times with replacements
# compareFaFiles(f1, f2)
#' @export
getNseqsWithReplace<-function(seqDir,
                   N,
                   nsplits=10,
                   outputPre,
                   perc=NULL) {

  for (i in 1:nsplits) {
    seed=Sys.time()
    # <outputPre>.split.1.fa
    tmpFile=paste0(outputPre,'_XXX_')
    ofile=getNseqs(seqDir=seqDir,
             N=N,
             nsplits=1,
             outputPre=tmpFile,
             perc=perc,
             seed=seed)
    file.rename(ofile, paste0(outputPre, '.split.',i,'.fa'))
  }
  return(paste0(outputPre, '.split.',1:nsplits,'.fa'))
}


#' getNseqs gets a combined fa file for each split following a given grams distribution.
#' @details
#' getNseqs tries to split given fa files without replacement.
#' However, when the sequence number of a gram is less than the needed number for all splits (defined by `nsplits`), the sequences of that gram in different splits may be duplicated
#' To meet the given `N` in each split, the fa file of the gram with largest sequence number will be used to supplement.
#' This function requires that the number of sequences in each fa file is enough for at least one split, otherwise, please reduce `N` or `nsplits`.
#' @param seqDir dir of fa files with each gram a file, e.g., AATAAA.fa.
#' @param N the number of sequences for each split (fa file).
#' @param nsplits number of splits (fa files), default is 10.
#' @param outputPre output dir and file name prefix, e.g., "./deepDataTrain" to output files like deepDataTrain.split.1/2/...fa.
#' @param perc a vector of gram proportion or gram counts, which stores the fraction of each gram in `seqDir`.
#' The names of the `perc` vector should contain all grams in `seqDir`.
#' If perc=NULL, then will calculate the fraction from `seqDir`.
#' If perc is count, then will calculate the percentages.
#' @param seed a random seed for shuffle the order of sequences of each fa file in `seqDir`.
#' @return a file list based on `outputPre`, e.g., <outputPre>.split.1.fa
#' @examples
#' \dontrun{
#' ## split to 10 splits each with 1000s
#' files=getNseqs(seqDir='realPA_Arich_seq_test/',
#'                N=1000, nsplits=10,
#'                outputPre="./realPA_Arich_seq_test_1000s_10splits",
#'                perc=NULL)
#' plotATCGforFAfile (faFiles=files, ofreq=FALSE, opdf=T, refPos=101,
#'              filepre='realPA_Arich_seq_test_1000s_10splits', mergePlots=TRUE)
#'
#' ## split IP files, using perc from real PA files
#' perc=statCntFas(path='polyAseqTrap/modelData/', filePre='realPA_Arich_seq.')
#' ## split to 2 splits each with 4000s
#' files=getNseqs(seqDir='IP_Arich_seq_train/',
#'                N=4000, nsplits=2,
#'                outputPre="./IP_Arich_seq_train_4000s_2splits",
#'                perc=perc)
#'
#' ## randomly generate 1 split with random seed,
#' ## which can be run multiple times with different seeds
#' ## to get different splits even there are some overlapping.
#' files=getNseqs(seqDir='realPA_Arich_seq_test/',
#'                N=10000, nsplits=1,
#'                outputPre="realPA_Arich_seq_test_10000s",
#'                perc=NULL, seed=9999)
#' }
#' @export
getNseqs<-function(seqDir,
                   N,
                   nsplits=10,
                   outputPre,
                   perc=NULL,
                   seed=123) {

  outputPre=formatPathFileName(outputPre, addBar = FALSE)
  dname=dirname(outputPre)
  if (!dir.exists(dname)) dir.create(dname)

  cat('###### Read PA fa\n')
  # realfas[[AATCCC]]=fa
  realfas=.readFas(seqDir, filePre='')

  nEach=unlist(lapply(realfas, length))
  totalN=sum(nEach)

  if (is.null(perc)) perc=nEach/totalN

  if (any(perc>1)) perc=perc/sum(perc)

  if (!all(names(nEach) %in% names(perc))) {
    cat("seq's grams:", names(nEach),"\n")
    cat("perc's grams:", names(perc),"\n" )
    stop("seq's grams not all in perc's!")
  }

  perc=perc[names(nEach)]

  nNeedEachSplit=floor(N*perc)
  nNeed=nNeedEachSplit*nsplits

  if (all(nNeed>nEach)) {
    g=names(nEach)[nNeed>nEach]
    print("seqs in seqDir:")
    print(nEach[g])
    print("seqs needed:")
    print(nNeed[g])
    stop("seq number for all grams are not enough for sampling!")
  }

  if (any(nNeed>nEach)) {
    g=names(nEach)[nNeed>nEach]
    print("seq number of some grams are not enough for sampling, will sampled with replace for these grams")
    print("seqs in seqDir:")
    print(nEach[g])
    print("seqs needed:")
    print(nNeed[g])
  }

  # Put the most at the end and make up for the missing parts
  nEach=sort(nEach, decreasing = FALSE)
  nNeed=nNeed[names(nEach)]
  nNeedEachSplit=nNeedEachSplit[names(nEach)]
  realfas=realfas[names(nEach)]

  # output file: .../deepDataTrain.splitN.fa
  ofiles=paste0(outputPre, '.split.',1:nsplits,'.fa')
  for (i in 1:nsplits) {
    if (file.exists(ofiles[i])) file.remove(ofiles[i])
  }

  # count each split
  ns=rep(0, nsplits)

  for (g in names(realfas)) { # names(realfas)
    fa=realfas[[g]]

    # random order
    set.seed(seed)
    fa=fa[sample(x=1:length(fa), size=length(fa), replace=FALSE)]

    # Divided into N parts, but each part is not yet nNeedEachSplit
    if (nsplits==1)
      breaks=rep(1, length(fa))
    else
      breaks=cut(seq(1, length(fa)), breaks=nsplits, labels=FALSE)

    for (i in 1:nsplits) {

      ids=which(breaks==i)
      idsAll=ids

      # not enough breaks
      if (length(ids)==0) {
        if (nNeedEachSplit[g]>length(fa)) {
          cat(sprintf("Not enough seqs (%d) for even %s 1 split (need %d), will use all seqs for each split\n", length(fa), g, nNeedEachSplit[g]))
          ids=1:length(fa)
        } else {
          cat(sprintf("Not enough seqs for %s split %d (need %d), will sample from the whole for this split\n", g, i, nNeedEachSplit[g]))
          ids=sample(1:length(fa), size=nNeedEachSplit[g])
        }
      } else {
        # After splitting each original data, take nNeedEachSplit from each split ID
        ids=ids[1:min(length(ids), nNeedEachSplit[g])]
      }

      ns[i]=ns[i]+length(ids)

      # Finally fill in the blanks - take the last few sequences of the split
      if (g==names(realfas)[length(realfas)]) {
        if (ns[i]!=N) {
          nlast=N-ns[i]
          cat(sprintf("Split %d: Get %d seqs from %s to make N=%d\n", i, nlast, g, N))

          if (length(idsAll)<nNeedEachSplit[g]+nlast) {
            cat(sprintf("No extra seqs from %s for split %d to make N=%d, please set [perc]\n", g, i, N))
            stop()
          }
          ids=c(ids, idsAll[(length(idsAll)-nlast+1):length(idsAll)])
        }
      }

      fa1=fa[ids]
      seqinr::write.fasta(sequences = fa1,
                          names = names(fa1),
                          nbchar = 1000,
                          file.out = ofiles[i],
                          open="a")

    }

  }

  cat('>>>', dname, ':\n')
  print(basename(ofiles))
  return(ofiles)
}





##Randomly select train and test sequences from a set of realPA or IP sequences
##If ipPre=NULL, only the files specified by realPre will be extracted, meaning that the output sequence will all end with: 1
##If ipPre provides, draw trainN real and trainN IP addresses; TestN real and IP (ending with real=: 1 and IP=: 0 respectively)
##If ipPre=NULL, draw trainN real and testN real
##Output:< Output Pre>. rain. fa and<output Pre>. test. fa
## 1. get real and IP (trainN 1 and trainN 0; testN 1 and testN 0)
#files=getTrainTestSeqs(trainN=4000, testN=4000, outputPre="deepData",
#                       realPre='realPA_Arich_seq.', ipPre='IP_Arich_seq.',
#                       path='./', seed=123)
#plotATCGforFAfile (faFiles=files, ofreq=TRUE, opdf=T, refPos=101, filepre='deepData_base_profile', mergePlots=TRUE)

## 2. get only real (We will also draw both trainN and testN at the same time, but both are 1)
#files=getTrainTestSeqs(trainN=4000, testN=4000, outputPre="deepData",
#                       realPre='realPA_Arich_seq.', ipPre=NULL,
#                       path='./', seed=123)
#plotATCGforFAfile (faFiles=files, ofreq=TRUE, opdf=T, refPos=101, filepre='deepData_base_profile', mergePlots=TRUE)
.getTrainTestSeqs<-function(trainN=4000,
                           testN=4000,
                           outputPre="deepData",
                           realPre,
                           ipPre=NULL,
                           path='./',
                           seed=123) {

  if (file.exists(paste0(outputPre, '.train.fa'))) file.remove(paste0(outputPre, '.train.fa'))
  if (file.exists(paste0(outputPre, '.test.fa'))) file.remove(paste0(outputPre, '.test.fa'))

  cat('###### Read real PA fa\n')
  # realfas[[AATCCC]]=fa
  realfas=.readFas(path, filePre=realPre)

  if (!is.null(ipPre)) {
    cat('###### Read IP fa\n')
    ipfas=.readFas(path, filePre=ipPre)

    # get common grams between real and ip fa files
    cmgrams=intersect(names(realfas), names(ipfas))
    realfas=realfas[cmgrams]
    ipfas=ipfas[cmgrams]

    cat(sprintf("%d common grams between real~IP\n", length(cmgrams)))

    ipLen=unlist(lapply(ipfas, length))
    ipN=sum(ipLen)
  }


  realLen=unlist(lapply(realfas, length))
  realN=sum(realLen)

  selN=trainN+testN

  if (!is.null(ipPre))
    cat(sprintf("Randomly select %d*2 train seqs and %d*2 test seqs from %d realPAs and %d IPs", trainN, testN, realN, ipN))
  else
    cat(sprintf("Randomly select %d train seqs and %d test seqs from %d realPAs", trainN, testN, realN))

  set.seed(seed)
  # Regarding the real seq: Train N and test N, draw equal proportions of seq for each gram, and add \ 1 at the end
  per1=trainN/realN; per2=testN/realN
  sum1=0; sum2=0
  sampleNs=data.frame(gram=names(realfas), train=0, test=0)
  rownames(sampleNs)=sampleNs$gram
  sampleNs$gram=NULL
  cat('###### Randomly select real PAs (add :1)...\n')
  cat(sprintf("train percent=%f; test percent=%f", per1, per2))

  file1=paste0(outputPre, '.train.fa')
  file2=paste0(outputPre, '.test.fa')

  for (g in names(realfas)) {
    fa=realfas[[g]]
    n1=floor(per1*length(fa))
    idx1=sample(1:length(fa), size=n1, replace = FALSE)

    n2=floor(per2*length(fa))
    idx2=sample((1:length(fa))[-idx1], size=floor(per2*length(fa)), replace = FALSE)

    sum1=sum1+length(idx1)
    sum2=sum2+length(idx2)

    cat(sprintf("%s train=%d test=%d; total train=%d; total test=%d\n", g, length(idx1), length(idx2), sum1, sum2))
    fa1=fa[idx1] #AATAAA train real
    fa2=fa[idx2] #AATAAA test real

    # Add \ 1 to the end of true and add a gram tag (previously output sequence did not have a gram tag, IP has a gram tag)
    seqinr::write.fasta(sequences = fa1,
                names = paste0(names(fa1),';', g,':1'),
                nbchar = 1000,
                file.out = file1,
                open="a")

    seqinr::write.fasta(sequences = fa2,
                names = paste0(names(fa2),';', g,':1'),
                nbchar = 1000,
                file.out = file2,
                open="a")

    # The last union may not exactly train N, so select a few more seqs from the last gram for training and testing
    idx3=c(); idx4=c()
    if (names(realfas)[length(realfas)]==g) {
      if (sum1!=trainN) {
        if (sum1<trainN) {
          nadd=trainN-sum1
          idx3=sample((1:length(fa))[-c(idx1,idx2)], size=nadd, replace = FALSE)
          sum1=sum1+length(idx3)
          cat(sprintf("Add %s train=%d; total train=%d\n", g, length(idx3), sum1))

          fa1=fa[idx3] #AATAAA train real
          seqinr::write.fasta(sequences = fa1,
                      names = paste0(names(fa1),';', g,':1'),
                      nbchar = 1000,
                      file.out = file1,
                      open="a")
        } else {
          stop("Error, sum train > trainN!\n")
        }
      }

      if (sum2!=testN) {
        if (sum2<testN) {
          nadd=testN-sum2
          # To remove training, testing, and supplementary sets that have already been drawn
          idx4=sample((1:length(fa))[-c(idx1,idx2, idx3)], size=nadd, replace = FALSE)
          sum2=sum2+length(idx4)
          cat(sprintf("Add %s test=%d; total test=%d\n", g, length(idx4), sum2))

          fa1=fa[idx4] #AATAAA train real
          seqinr::write.fasta(sequences = fa1,
                      names = paste0(names(fa1),';', g,':1'),
                      nbchar = 1000,
                      file.out = file2,
                      open="a")
        } else {
          stop("Error, sum test > testN!\n")
        }
      }
    } #last add some seqs

    # SampleNs records the number of samples drawn for each pair
    sampleNs[g, ]=c(length(idx1)+length(idx3), length(idx2)+length(idx4))
  }

  if (!is.null(ipPre)) {
    cat('###### Randomly select IPs (add :0)...\n')
    sum1=0; sum2=0
    # Regarding IP seq: TrainN and testN, each gram draws the same sequence as the real seq
    for (g in rownames(sampleNs)) {
      fa=ipfas[[g]]
      n1=sampleNs[g, 'train']
      n2=sampleNs[g, 'test']

      idx1=sample(1:length(fa), size=n1, replace = FALSE)
      idx2=sample((1:length(fa))[-idx1], size=n2, replace = FALSE)
      sum1=sum1+length(idx1)
      sum2=sum2+length(idx2)

      cat(sprintf("%s train=%d test=%d; total train=%d; total test=%d\n", g, length(idx1), length(idx2), sum1, sum2))
      fa1=fa[idx1] #AATAAA train real
      fa2=fa[idx2] #AATAAA test real

      seqinr::write.fasta(sequences = fa1,
                          names = paste0(names(fa1),':0'),
                          nbchar = 1000,
                          file.out = file1,
                          open="a")

      seqinr::write.fasta(sequences = fa2,
                          names = paste0(names(fa2),':0'),
                          nbchar = 1000,
                          file.out = file2,
                          open="a")
    }
    cat(sprintf("Output %d train seqs (1) and %d test seqs (1) to <%s>, <%s>", trainN, testN, paste0(outputPre, '.train.fa'), paste0(outputPre, '.test.fa')))
  } else {
    cat(sprintf("Output %d*2 train seqs (1/0) and %d*2 test seqs (1/0) to <%s>, <%s>", trainN, testN, paste0(outputPre, '.train.fa'), paste0(outputPre, '.test.fa')))
  }
  return(c(file1, file2))
}


#' combineFaFiles combines multiple fa files into one file, allowing adding grams and labels to the sequence title.
#' @param fafiles fa file names, or file names with only sequences (each line one sequence).
#' @param ofile output fa file name.
#' @param addGram TRUE to add the gram in the sequence file name to the sequence title. Default is FALSE.
#' @param label If provided (e.g., ":1"), then add label to the end of the sequence title. Default is NULL.
#' @param verbose if TRUE, then show msg about seq numbers in fa files. Default is FALSE.
#' @return number of sequences in `ofile`.
#' @examples
#' \dontrun{
#' ## combine fa files
#' combineFaFiles(fafiles=list.files('intron_realPA_Arich_seq', full.names=TRUE),
#'                ofile='../modelDataSplits/intron_realPA_Arich_seq.all.fa',
#'                addGram = TRUE, label=":1",
#'                verbose=TRUE)
#' ## combine sequence files (each line one seq)
#' }
#' @export
combineFaFiles<-function(fafiles, ofile,
                         addGram=FALSE, label=NULL,
                         verbose=FALSE) {

  if (length(fafiles)==0) {
    cat("fafiles is empty!\n")
    return(0)
  }

  if (addGram) {
    cms=movAPA:::.autoDetectCommonString(fafiles, sbj=NULL, beAll = TRUE)
    grams=movAPA:::.removeCommonStr(fafiles, cms, fixed=TRUE)
    if (max(nchar(grams))>10) stop("addGram=TRUE, but fafiles seem not contain gram: ", grams[1], '\n')
  }

  if (is.null(label)) label=''
  N=0

  unlink(ofile)

  for (f in fafiles) {

    con <- file(f, "r")
    str=readLines(con, n=1)
    close(con)

    if (substr(str,1,1)=='>') {
      fa=seqinr::read.fasta(f, forceDNAtolower = F, as.string = TRUE)
      names=names(fa)

      if (addGram) {
        names = paste0(names,';', grams[fafiles==f], label)
      } else {
        names = paste0(names, label)
      }

      seqinr::write.fasta(sequences = fa,
                          names = names,
                          nbchar = 1000,
                          file.out = ofile,
                          open='a')

    } else {
      fa=read.table(f, header=F)$V1
      fa=as.list(fa)
      names=paste0('seq', 1:length(fa),label)
      seqinr::write.fasta(sequences = fa,
                          names = names,
                          nbchar = 1000,
                          file.out = ofile,
                          open='a')
    }

    N=N+length(fa)
    if (verbose) cat(f, '\t', length(fa),'\n')
  }

  if(verbose) {
    cat('Total\t', N, '\n')
  }
  return(N)
}



#' combineFaFiles_fast combines multiple fa files into one file using readLines.
#' @param fafiles fa file names.
#' @param ofile output fa file name.
#' @param verbose if TRUE, then show msg about seq numbers in fa files. Default is FALSE.
#' @return number of sequences in `ofile`.
#' @export
combineFaFiles_fast<-function(fafiles, ofile, verbose=FALSE) {
  N=rep(0, length(fafiles))
  names(N)=fafiles
  con <- file(description = ofile, open = 'w')
  for (f in fafiles) {
    txt=readLines(f, n = -1)
    writeLines(txt, con=con)
    N[f]=length(txt)/2
  }
  close(con)
  if(verbose) {
    for (i in 1:length(N)) {
      cat(names(N)[i], '\t', N[i],'\n')
    }
    cat('Total\t', sum(N), '\n')
  }
  return(sum(N))
}


#' compareFaFiles counts overlapping of sequence titles between two fa files.
#' @param f1 the first fa file.
#' @param f2 the second fa file.
#' @return number of sequences and overlapped sequences.
#' @examples
#' \dontrun{
#' compareFaFiles('1.fa','2.fa')
#' }
#' @export
compareFaFiles<-function(f1, f2) {
  fa1=seqinr::read.fasta(f1, forceDNAtolower = F)
  fa2=seqinr::read.fasta(f2, forceDNAtolower = F)
  ovp=length(intersect(names(fa1), names(fa2)))
  cat("n1: ",length(fa1),'\n')
  cat("n2: ",length(fa1),'\n')
  cat("ovp:", ovp,'\n')
  return(c(n1=length(fa1), n2=length(fa2), overlapping=ovp))
}


#' getDLdata gets a combined fa file from true and false seq dir.
#' @details
#' If exactMatch=TRUE, then splits number and index in tSplitDir and fSplitDir should be matched.
#' @param tSplitDir dir of fa files. This function does not check the label (T/F), just combine two dirs.
#' @param fSplitDir dir of fa files.
#' @param outputPre output dir and file name prefix, e.g., "./train.T.F" to output files like train.T.F.1/2/...fa.
#' @param exactMatch default is TRUE.
#' If exactMatch=TRUE, then splits number and index in tSplitDir and fSplitDir should be matched.
#' If it is FALSE, then will rep the split index of tSplitDir/fSplitDir with less splits to make the number of splits equal.
#' @return NULL. Will create the output dir and output the fa files.
#' @examples
#' \dontrun{
#' getDLdata(tSplitDir='realPA_Arich_seq_train_1000s_10splits',
#'           fSplitDir='IP_Arich_seq_train_1000s_10splits',
#'           outputPre="train.1000T.1000F/train.1000T.1000F",
#'           exactMatch=TRUE)
#' }
#' @export
getDLdata<-function(tSplitDir, fSplitDir, outputPre="./train.T.F", exactMatch=TRUE) {

  .getSplitIdx<-function(files) {
    i=as.numeric( gsub("\\D", "", substr(files, nchar(files)-5, nchar(files)) ) )
    names(i)=files
    return(i)
  }

  outputPre=formatPathFileName(outputPre, addBar = FALSE)
  dname=dirname(outputPre)
  if (!dir.exists(dname)) dir.create(dname)

  tfiles=list.files(tSplitDir, full.names = TRUE) # **split.1.fa
  ti=.getSplitIdx(tfiles)

  ffiles=list.files(fSplitDir, full.names = TRUE) # **split.1.fa
  fi=.getSplitIdx(ffiles)

  ti=sort(ti)
  fi=sort(fi)

  if (exactMatch) {
    if (length(ti)!=length(fi)) stop( sprintf("split number not equal: tSplitDir (%d) != fSplitDir (%d)!", length(ti), length(fi)) )
    if (sum(ti-fi)!=0) stop( sprintf("split index of tSplitDir and fSplitDir not matched!"))
  }

  # Not long enough, make up for it
  n=max(length(ti), length(fi))
  ti=rep(ti, length.out=n)
  fi=rep(fi, length.out=n)

  for (i in 1:n) {
	  ofile=paste0(outputPre,'.', i, '.fa')
	  nseq=combineFaFiles(c(names(ti)[i], names(fi)[i]), ofile=ofile)
    cat(sprintf(">>> %s (%d seqs from tSplit %d + fSplit %d)\n", ofile, nseq, ti[i], fi[i]))
  }

}


## ----- stat trained results -----

#' plotFaDeepRes plots single nucleotide profile and output fasta files for TP/FP/TN/FN.
#' @param deepResCsvs the csv file name from DeepIP_test.py.
#' @param fafile the fa file corresponding to deepResCsvs
#' @param outputFa NULL for not outputing fa but just plotting. Otherwise, the prefix for the output fa files.
#' @return fa file names.
#' Output fasta file as follows:<outputFa or deepResCsv>.FN10.fa means that 1 is judged as 0, indicating FN.
#' @examples
#' \dontrun{
#' plotFaDeepRes('human_on_chlamy.csv',
#'               'chlamy.fa',
#'               outputFa='fa/human_on_chlamy')
#' }
#' @export
plotFaDeepRes<-function(deepResCsv, fafile, outputFa=NULL) {
  d=read.csv(deepResCsv)
  d$title=gsub('^>','',d$title)

  deepResCsv=formatPathFileName(deepResCsv, addBar = FALSE, noExt = TRUE)

  fa=seqinr::read.fasta(fafile, forceDNAtolower = F)

  if (!all(d$title %in% names(fa))) {
    n=sum(d$title %in% names(fa))
    stop(sprintf("Only (%d from %d) title from deepResCsv are in fafile!\n", n, nrow(d)))
  }

  fa=fa[d$title]

  str=c(FN='10', TP='11', FP='01', TN='00')
  ofiles=c()
  for (res in unique(d$res)) {
    fa2=fa[d$res==res]
    cat(sprintf("%s: %d\n", res, length(fa2)))

    if (!is.null(outputFa)) {
      ofile=paste0(outputFa,'.', res, str[res], '.fa')
    } else {
      ofile=paste0(deepResCsv,'.', res, str[res], '.fa')
    }
    ofiles=c(ofiles, ofile)
    seqinr::write.fasta(sequences = fa2,
                        name=names(fa2),
                        nbchar = 1000,
                        file.out = ofile,
                        open="w")
  }

  movAPA::plotATCGforFAfile (faFiles=ofiles, ofreq=FALSE, opdf=T, refPos=101, filepre=paste0(deepResCsv,'_confusion.base_profile'), mergePlots=TRUE)
  if (is.null(outputFa)) unlink(ofiles)
  return(ofiles)
}

#' statDeepRes gets metrics for .csv prediction results from DeepIP_test.py.
#' @param deepResCsvs a string vector storing one or multiple csv files from DeepIP_test.py.
#' If the name of the vector is provided, then will use the name as the output row title denoting each csv file.
#' @param ofile output .csv or other file to store the metrics table. If it is NULL, then not output.
#' @param verbose TRUE to print the metrics table.
#' @param shortFileName TRUE to remove common strings in file names of deepResCsvs.
#' @return A dataframe of the metrics, with each row one file of `deepResCsvs` and each column a metric.
#' Metrics are: TP_11	FP_01	TN_00	FN_10	Sensitivity	Specificity	Precision	Recall	F1	ROC	AUC.
#' @examples
#' \dontrun{
#' statDeepRes('train.10000T.10000F.1.epoch100_on_test.all.csv',
#'              ofile='train.10000T.10000F.1.epoch100_on_test.all.stat.txt')
#' }
#' @export
statDeepRes<-function(deepResCsvs, ofile=NULL, verbose=FALSE, shortFileName=TRUE) {

  rnames=names(deepResCsvs) # name provided
  if (is.null(rnames)) { # not provided, use file name
    rnames=formatPathFileName(basename(deepResCsvs), noExt=TRUE)

    if (length(deepResCsvs)>1 & shortFileName) {
      fname=rnames
      cms=movAPA:::.autoDetectCommonString(fname, sbj=NULL, beAll = TRUE)
      fname=movAPA:::.removeCommonStr(fname, cms, fixed=TRUE)
      rnames=fname
    }

  }

  metrics=c('TP_11','FP_01','TN_00','FN_10', 'Sensitivity', 'Specificity', 'Precision', 'Recall', 'F1', 'ROC', 'AUC')
  stats=data.frame(matrix(data=NA, nrow=length(metrics), ncol=length(deepResCsvs),
                          dimnames = list(metrics,
                                          rnames
                                          )))
  for (f in deepResCsvs) {

    d=read.csv(f)

    if (!all(c('res','predict_label','true_label','score') %in% colnames(d)))
      stop("res, predict_label, true_label, score not all in ",f,"\n")

    tpfp=table(d$res)
    tpfp[c('TP','FP','TN','FN')]=tpfp[c('TP','FP','TN','FN')]
    tpfp[is.na(tpfp)]=0
    stats[c('TP_11','FP_01','TN_00','FN_10'), which(deepResCsvs==f)]=tpfp[c('TP','FP','TN','FN')]

      f1=factor(d$predict_label, levels=c('0','1')) #!add levels
      f2=factor(d$true_label, levels=c('0','1'))

      cm=caret::confusionMatrix(f1, f2, mode = "everything", positive="1")
      cm=cm$byClass
      m=intersect(metrics, names(cm))
      stats[m, which(deepResCsvs==f)]=cm[m]

      # https://topepo.github.io/caret/measuring-performance.html#lift-curves
      d1=d[, c('true_label','predict_label', 'score', 'score')]
      colnames(d1)=c('obs','pred', '1', '0')
      d1$`0`=1-d1$`1`
      d1$obs=factor(d1$obs, levels=c('0','1')); d1$pred=factor(d1$pred, levels=c('0','1'))
      roc=caret::twoClassSummary(d1, lev=levels(d1$obs))['ROC'] # ROC:0.9965681
      auc=caret::prSummary(d1, lev = levels(d1$obs))['AUC'] # AUC: 0.9962009
      stats[c('ROC','AUC'), which(deepResCsvs==f)]=c(roc, auc)
    #}

    # Metrics::accuracy(d$true_label, d$predict_label)
    # Metrics::precision(d$true_label, d$predict_label)
    # Metrics::f1(d$true_label, d$predict_label) # 1 (py=0.9720?)
    # Metrics::recall(d$true_label, d$predict_label)
    #
    # # https://www.statology.org/auc-in-r/
    # Metrics::auc(d$true_label, d$score)  #1=pos, 0=neg; score=large prob for pos; 0.9965681
    # pROC::auc(d$true_label, d$score) #0.9966 (py=0.9966)

    # r=pROC::roc(d$true_label, d$score)
    # plot(r, type = "S",col = "red")

  }

  stats=format(stats, scientific=FALSE)
  stats=as.data.frame(t(stats))
  for (i in 1:ncol(stats)) stats[, i]=as.numeric(stats[,i])

  if (verbose) print(stats)

  if (is.null(ofile)) {
    return(stats)
  }
  if (grepl('csv$',ofile)) {
    write.csv(stats, file=ofile)
    return(stats)
  }
  write.table(stats, file=ofile, quote=F, sep="\t")
  return(stats)
}

#' statDeepResByGrams gets metrics for each gram for a .csv prediction file from DeepIP_test.py.
#' @param deepResCsv one csv files from DeepIP_test.py, with the title column like ">PA615:chr1;+;16972964;AACAAA:1".
#' Or if the `gram` column is already in the csv, then will be used instead of parsing the title column.
#' @param ofile output .csv or other file to store the metrics table. If it is NULL, then not output.
#' @param verbose TRUE to print the metrics table.
#' @param plot TRUE to plot F1 ROC Sn Sp to compare grams.
#' @return A dataframe of the metrics, with each column one gram in file of `deepResCsvs` and each row a metric.
#' Metrics are: TP_11	FP_01	TN_00	FN_10	Sensitivity	Specificity	Precision	Recall	F1	ROC	AUC.
#' grams are obtained from the last chars of the seq title like ">PA615:chr1;+;16972964;AACAAA:1".
#' In the bar plot, the grams are ordered by F1 desc. AATAAA, ATTAAA, and NOPAS are highlighted in darkgrey.
#' @examples
#' \dontrun{
#' stats=statDeepResByGrams('train.10000T.10000F.1.epoch100_ON_test.all.csv',
#'                          ofile=NULL,
#'                          plot=TRUE)
#' }
#' @export
statDeepResByGrams<-function(deepResCsv, ofile=NULL,
                             verbose=FALSE, plot=TRUE) {
  d=read.csv(deepResCsv)

  if (!('gram' %in% colnames(d))) {
    title=strsplit(d$title, split=':|;')
    title=unlist(lapply(title, function(par) {par[length(par)-1]}))
    if (max(nchar(title))>10) stop("the sequence title seem not contain gram (the second last field): ", title[1], '\n')
    d$gram=title
  }

  #d$gram=substr(d$title, nchar(d$title)-7, nchar(d$title)-2)
  ds=split(d, factor(d$gram))
  for (i in 1:length(ds)) {
    g=names(ds)[i]
    write.csv(ds[[i]], file=paste0(g,'.csv'))
  }
  files=paste0(names(ds),'.csv')
  ns=unlist(lapply(ds, nrow))
  files=files[order(ns, decreasing = TRUE)]
  stat=statDeepRes(files, ofile=ofile, verbose=verbose)
  unlink(files)

  if (plot) {
    stats=stat
    stats=stats[, c('Sensitivity','Specificity','F1','ROC'), ]
    stats$gram=rownames(stats)
    stats=tidyr::pivot_longer(stats, cols=!gram, names_to = 'metric', values_to='score')
    stats$score=as.numeric(stats$score)

    # order grams by F1, and highlight AATAAA and ATTAAA, NOPAS
    gorder=stats[stats$metric=='F1', ]
    gorder$gram=forcats::fct_reorder(gorder$gram, gorder$score, .desc = TRUE)

    stats$gram=factor(stats$gram, levels=levels(gorder$gram))

    g=ggplot2::ggplot(data=stats, aes(x=gram, y=score, group=gram, fill=gram)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values=c("AATAAA"='darkgray',"ATTAAA"='darkgray',"NOPAS"='darkgray')) +
      facet_wrap( ~ metric, ncol=2) +
      theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
            strip.text = element_text(size=14, face="bold"),
            legend.position="none")
    print(g)
  }
  return(stat)
}


deepPASTARes2Csv<-function(resFiles) {
  for (f in resFiles) {
    ofile=gsub('txt$', 'csv', f)

    d=read.table(f, header=F)
    if (!ncol(d)==2) stop('DeepPASTA results should have only two columns without header (title, prob)')
    colnames(d)=c('title','score')

    # >PA51504_chr11;-;84093923;84098722;nonPA-pos=3901;TATAAA_0;TATAAA_0	0.63260937
    true_label=substr(d$title, nchar(d$title), nchar(d$title))
    predict_label=ifelse(d$score>0.5, 1, 0)
    res='TP'
    res[true_label==1 & predict_label==0]='FN'
    res[true_label==1 & predict_label==1]='TP'
    res[true_label==0 & predict_label==1]='FP'
    res[true_label==0 & predict_label==0]='TN'

    d$true_label=true_label; d$predict_label=predict_label; d$res=res

    write.csv(d, file=ofile)
  }
  return(gsub('txt$', 'csv', resFiles))
}


#' kcountByPos counts kgrams at each positions in a fasta file.
#' @param fafile a fa file name
#' @param grams vector of grams
#' @param k integer, e.g., 6
#' @param from if NA, then from=1
#' @param to if NA, then to=min length in fafile
#' @param sort if TRUE, then sort grams in output table
#' @param topn if not NULL, then return topn grams with max counts
#' @param perc if TRUE, to output percentage rather than counts
#' @return a count or perc table with rows being grams and cols being positions from [`from` .. (`to`-`K`+1)].
#' @examples
#' \dontrun{
#' cnts=kcountByPos(fafile='test.fa', grams=c('GCCACA','GAAGTC'),
#'                  from=90, to=110, sort=T, topn=50, perc=FALSE)
#' cnts=kcountByPos(fafile='test.fa', k=6, from=90, to=110,
#'                  sort=T, topn=50, perc=FALSE)
#' cnts=kcountByPos(fafile='test.fa', grams=movAPA:::getVarGrams('MM'),
#'                  from=90, to=110, sort=T, topn=50, perc=TRUE)
#' }
#' @export
kcountByPos <- function (fafile=NULL,
                         grams=NULL, k=NULL,
                         from=NA, to=NA,
                         sort=TRUE, topn=50, perc=FALSE) {

  if ( (!is.null(grams) & !is.null(k))  ) stop("kcount: grams or k, only one option")

  grams=movAPA:::getVarGrams(grams)

  if (!is.null(k)) grams=movAPA:::genKgrams(k)
  klen=min(nchar(grams))

  seqs=Biostrings::readDNAStringSet(fafile, format="fasta")

  if (is.na(from)) from=1
  if (is.na(to)) {
    mlen=min(Biostrings::width(seqs))
    to=mlen
  }

  cat('kcount from',from,'to',to,'\n')

  # call movAPA::kcount to count grams in each k-range
  cnts=matrix(data=NA, nrow=length(grams), ncol=to-klen+1-from+1)
  i=1
  for (fi in from:(to-klen+1)) {
    cnt1=movAPA::kcount(seq=seqs, grams=grams, from=fi, to=fi+klen-1, sort=F, topn=NULL, perc=FALSE)
    cnts[, i]=cnt1[, 2]
    i=i+1
  }
  cnts=as.data.frame(cnts)
  rownames(cnts)=grams
  colnames(cnts)=paste0('pos',from:(to-klen+1))

  if (sort | !is.null(topn)) cnts=cnts[order(rowSums(cnts), decreasing = TRUE),]

  if (perc) {
    perc=sapply(cnts, function(par) par/sum(par))
    if (length(perc)==1) {
      names(perc)='perc'
    } else {
      names(perc)=colnames(cnts)
    }
    cnts=perc
  }
  rownames(cnts)=grams

  if (!is.null(topn)) cnts=cnts[1:min(topn, nrow(cnts)), ,drop=F]

  return(cnts)
}



