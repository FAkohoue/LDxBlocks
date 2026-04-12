# ==============================================================================
# haplotypes.R - Comprehensive haplotype analysis for LDxBlocks
# Sections: 1-Phasing  2-Extraction  3-Diversity  4-QTL  5-Matrix  6-Writers
# ==============================================================================
.norm_chr_hap <- function(x) sub("^chr","",as.character(x),ignore.case=TRUE)
`%||%` <- function(x,y) if(is.null(x)) y else x

#' Read Pre-Phased VCF
#' @param vcf_file Path to phased VCF or VCF.gz with 0|1 GT fields.
#' @param min_maf Minimum MAF. Default 0.0.
#' @param verbose Logical. Default TRUE.
#' @return List: hap1, hap2 (SNPs x individuals, 0/1), dosage (0/1/2), snp_info, sample_ids.
#' @export
read_phased_vcf <- function(vcf_file, min_maf=0.0, verbose=TRUE) {
  if (!file.exists(vcf_file)) stop("VCF not found: ", vcf_file, call.=FALSE)
  if (verbose) message("[read_phased_vcf] Reading: ", basename(vcf_file))
  con <- gzfile(vcf_file,"r"); on.exit(try(close(con),silent=TRUE),add=TRUE)
  header <- NULL; data_lines <- character(0)
  repeat {
    ln <- readLines(con, n=1L, warn=FALSE); if (!length(ln)) break
    if (startsWith(ln,"##")) next
    if (startsWith(ln,"#CHROM")) { header <- ln; next }
    data_lines <- c(data_lines, ln)
  }
  if (is.null(header)) stop("No #CHROM header.", call.=FALSE)
  cols <- strsplit(header,"\t",fixed=TRUE)[[1L]]
  sids <- cols[10:length(cols)]; ns <- length(sids); nv <- length(data_lines)
  if (verbose) message("[read_phased_vcf] ",nv," variants x ",ns," samples")
  h1 <- matrix(NA_real_,nv,ns); h2 <- matrix(NA_real_,nv,ns)
  sc <- character(nv); sp <- integer(nv); sid <- character(nv)
  sr <- character(nv); sa <- character(nv)
  for (i in seq_len(nv)) {
    f <- strsplit(data_lines[i],"\t",fixed=TRUE)[[1L]]
    sc[i] <- .norm_chr_hap(f[1]); sp[i] <- as.integer(f[2])
    sid[i] <- f[3]; sr[i] <- f[4]; sa[i] <- f[5]
    fi <- which(strsplit(f[9],":",fixed=TRUE)[[1]]=="GT")
    for (j in seq_len(ns)) {
      gt <- strsplit(f[9+j],":",fixed=TRUE)[[1]][fi]
      al <- strsplit(gt,"[|/]")[[1]]
      if (length(al)==2L && !any(al%in%c(".",NA))) {
        h1[i,j] <- suppressWarnings(as.integer(al[1]))
        h2[i,j] <- suppressWarnings(as.integer(al[2]))
      }
    }
  }
  rownames(h1) <- rownames(h2) <- sid; colnames(h1) <- colnames(h2) <- sids
  dos <- h1+h2
  snp_info <- data.frame(SNP=sid,CHR=sc,POS=sp,REF=sr,ALT=sa,stringsAsFactors=FALSE)
  if (min_maf>0) {
    af <- rowMeans(dos,na.rm=TRUE)/2; maf <- pmin(af,1-af)
    ok <- !is.na(maf)&maf>=min_maf
    h1 <- h1[ok,,drop=FALSE]; h2 <- h2[ok,,drop=FALSE]
    dos <- dos[ok,,drop=FALSE]; snp_info <- snp_info[ok,]
    if (verbose) message("[read_phased_vcf] After MAF>=",min_maf,": ",sum(ok))
  }
  list(hap1=h1, hap2=h2, dosage=dos, snp_info=snp_info, sample_ids=sids, phased=TRUE)
}

#' Statistical Phasing via Beagle 5.x
#' @param input_vcf Unphased VCF path.
#' @param out_prefix Output prefix (Beagle appends .vcf.gz).
#' @param beagle_jar Path to beagle.jar. Default: auto-detected.
#' @param java_path Java executable. Default "java".
#' @param nthreads Threads. Default 1L.
#' @param ref_panel Optional phased reference VCF. Default NULL.
#' @param beagle_args Additional Beagle arguments string.
#' @param verbose Logical. Default TRUE.
#' @return Path to phased VCF.gz. Load with read_phased_vcf().
#' @references Browning et al. (2018) Am J Hum Genet 103:338-348.
#' @export
phase_with_beagle <- function(input_vcf, out_prefix, beagle_jar=NULL,
                              java_path="java", nthreads=1L,
                              ref_panel=NULL, beagle_args="", verbose=TRUE) {
  if (!file.exists(input_vcf)) stop("VCF not found: ",input_vcf,call.=FALSE)
  if (is.null(beagle_jar)) {
    cands <- c(Sys.which("beagle.jar"),"/usr/local/bin/beagle.jar",
               file.path(Sys.getenv("HOME"),"beagle.jar"))
    found <- cands[nzchar(cands)&file.exists(cands)]
    if (!length(found))
      stop("beagle.jar not found.\nDownload: https://faculty.washington.edu/browning/beagle/beagle.html",call.=FALSE)
    beagle_jar <- found[1]
  }
  cmd <- paste(java_path,"-jar",shQuote(beagle_jar),
               paste0("gt=",shQuote(input_vcf)),
               paste0("out=",shQuote(out_prefix)),
               paste0("nthreads=",nthreads),
               if (!is.null(ref_panel)) paste0("ref=",shQuote(ref_panel)) else "",
               beagle_args)
  if (verbose) message("[phase_with_beagle] ",cmd)
  ret <- system(cmd)
  if (ret!=0L) stop("Beagle exited with status ",ret,call.=FALSE)
  out <- paste0(out_prefix,".vcf.gz")
  if (!file.exists(out)) stop("Beagle output not found: ",out,call.=FALSE)
  if (verbose) message("[phase_with_beagle] Done: ",out)
  invisible(out)
}

#' Pedigree-Based Allele Transmission Phasing
#' @description Assigns gametic phase using Mendelian transmission within
#'   parent-offspring trios. Exact when parents are homozygous; ambiguous
#'   positions are NA or randomly resolved when resolve_het=TRUE.
#' @param dosage_mat Numeric matrix (individuals x SNPs, 0/1/2/NA). Rownames must match pedigree$id.
#' @param pedigree Data frame with columns id, sire, dam. Use NA or "0" for unknown parents.
#' @param resolve_het Randomly resolve ambiguous het transmissions. Default FALSE.
#' @param seed RNG seed when resolve_het=TRUE. Default 1L.
#' @param verbose Logical. Default TRUE.
#' @return List: hap1 (sire gamete), hap2 (dam gamete), dosage, n_resolved, n_ambiguous, phased=TRUE.
#' @export
phase_with_pedigree <- function(dosage_mat, pedigree, resolve_het=FALSE,
                                seed=1L, verbose=TRUE) {
  if (!is.matrix(dosage_mat)) dosage_mat <- as.matrix(dosage_mat)
  miss <- setdiff(c("id","sire","dam"),names(pedigree))
  if (length(miss)) stop("pedigree missing: ",paste(miss,collapse=","),call.=FALSE)
  ped <- as.data.frame(lapply(pedigree[c("id","sire","dam")],as.character),stringsAsFactors=FALSE)
  ped[ped%in%c("0","NA","")] <- NA_character_
  ni <- nrow(dosage_mat); nk <- ncol(dosage_mat)
  h1 <- matrix(NA_real_,ni,nk,dimnames=dimnames(dosage_mat))
  h2 <- matrix(NA_real_,ni,nk,dimnames=dimnames(dosage_mat))
  set.seed(seed); nr <- 0L; na_ <- 0L
  for (id in rownames(dosage_mat)) {
    row <- ped[ped$id==id,,drop=FALSE]
    if (!nrow(row)||is.na(row$sire)||is.na(row$dam)||
        !row$sire%in%rownames(dosage_mat)||!row$dam%in%rownames(dosage_mat)) {
      g <- dosage_mat[id,]; h1[id,!is.na(g)] <- ceiling(g[!is.na(g)]/2)
      h2[id,!is.na(g)] <- floor(g[!is.na(g)]/2); next
    }
    go <- dosage_mat[id,]; gs <- dosage_mat[row$sire,]; gd <- dosage_mat[row$dam,]
    for (k in seq_len(nk)) {
      if (is.na(go[k])||is.na(gs[k])||is.na(gd[k])) next
      sa <- unique(if(gs[k]==0)c(0,0) else if(gs[k]==1)c(0,1) else c(1,1))
      da <- unique(if(gd[k]==0)c(0,0) else if(gd[k]==1)c(0,1) else c(1,1))
      vl <- expand.grid(a1=sa,a2=da); vl <- vl[vl$a1+vl$a2==go[k],,drop=FALSE]
      if (nrow(vl)==1L){h1[id,k]<-vl$a1;h2[id,k]<-vl$a2;nr<-nr+1L}
      else if(nrow(vl)>1L&&resolve_het){pk<-vl[sample.int(nrow(vl),1L),];h1[id,k]<-pk$a1;h2[id,k]<-pk$a2;nr<-nr+1L}
      else{na_<-na_+1L}
    }
  }
  if (verbose) message("[phase_with_pedigree] Resolved: ",nr," | Ambiguous: ",na_)
  list(hap1=h1,hap2=h2,dosage=h1+h2,n_resolved=nr,n_ambiguous=na_,phased=TRUE)
}

#' Collapse Phased Gametes to 0/1/2 Dosage
#' @param phased_list List from read_phased_vcf() or phase_with_pedigree().
#' @return Numeric matrix (individuals x SNPs, 0/1/2/NA).
#' @keywords internal
unphase_to_dosage <- function(phased_list) {
  if (!all(c("hap1","hap2")%in%names(phased_list)))
    stop("Need hap1 and hap2 elements.",call.=FALSE)
  phased_list$hap1+phased_list$hap2
}

#' Extract Haplotype Strings Within LD Blocks
#' @description Builds haplotype strings per individual per LD block.
#'   Unphased mode: "012201" (diploid). Phased mode: "011|100" (gametes).
#'   Blocks are always defined within a single chromosome.
#' @param geno Dosage matrix (individuals x SNPs) OR phased list from read_phased_vcf()/phase_with_pedigree().
#' @param snp_info Data frame with SNP, CHR, POS.
#' @param blocks LD block data frame from run_Big_LD_all_chr() with CHR, start.bp, end.bp.
#' @param chr Optional chromosome filter. Default NULL (all).
#' @param min_snps Minimum SNPs per block. Default 3L.
#' @param na_char Missing allele character. Default ".".
#' @return Named list of character vectors (length=n_individuals), one per block.
#'   Carries block_info attribute.
#' @export
extract_haplotypes <- function(geno, snp_info, blocks,
                               chr=NULL, min_snps=3L, na_char=".") {

  # -- Input type detection --------------------------------------------------
  # Three paths:
  #  1. LDxBlocks_backend  - STREAMING: one chromosome at a time, never full
  #     genome in RAM. Each chromosome is extracted, all its blocks processed,
  #     then freed with gc(FALSE) before the next chromosome.
  #  2. Phased list (hap1/hap2) - already in RAM, process as before.
  #  3. Numeric dosage matrix - already in RAM, process as before.
  is_backend <- inherits(geno, "LDxBlocks_backend")
  isp <- !is_backend && is.list(geno) && all(c("hap1","hap2")%in%names(geno))

  if (is_backend) {
    # -- STREAMING PATH: chromosome by chromosome ---------------------------
    iids       <- geno$sample_ids
    si         <- geno$snp_info
    si$CHR     <- .norm_chr_hap(si$CHR)
    blocks$CHR <- .norm_chr_hap(blocks$CHR)
    if (!is.null(chr)) {
      chr    <- .norm_chr_hap(chr)
      blocks <- blocks[blocks$CHR %in% chr, ]
    }
    res <- list()
    bi  <- data.frame(block_id=character(),CHR=character(),start_bp=integer(),
                      end_bp=integer(),n_snps=integer(),phased=logical(),
                      stringsAsFactors=FALSE)
    for (cb in unique(blocks$CHR)) {
      chr_blocks <- blocks[blocks$CHR == cb, ]
      if (!nrow(chr_blocks)) next
      chr_idx <- which(si$CHR == cb)
      if (!length(chr_idx)) next

      # Extract ONLY this chromosome - RAM proportional to one chromosome
      chr_geno  <- read_chunk(geno, chr_idx)   # individuals x chr_SNPs
      chr_si    <- si[chr_idx, ]

      for (b in seq_len(nrow(chr_blocks))) {
        blk <- chr_blocks[b,]
        sb  <- as.numeric(blk$start.bp); eb <- as.numeric(blk$end.bp)
        blk_idx <- which(chr_si$POS >= sb & chr_si$POS <= eb)
        if (length(blk_idx) < min_snps) next
        bid <- paste0("block_",cb,"_",sb,"_",eb)
        hs <- vapply(seq_len(nrow(chr_geno)), function(i) {
          v <- as.character(chr_geno[i, blk_idx])
          v[is.na(chr_geno[i, blk_idx])] <- na_char
          paste(v, collapse="")
        }, character(1L))
        names(hs) <- iids
        res[[bid]] <- hs
        bi <- rbind(bi, data.frame(block_id=bid, CHR=cb,
                                   start_bp=as.integer(sb),
                                   end_bp=as.integer(eb),
                                   n_snps=length(blk_idx), phased=FALSE,
                                   stringsAsFactors=FALSE))
      }
      rm(chr_geno, chr_si); gc(FALSE)  # free chromosome before next
    }
    attr(res,"block_info") <- bi
    return(res)
  }

  if (isp) {
    gm <- t(geno$dosage); h1m <- t(geno$hap1); h2m <- t(geno$hap2)
    iids <- geno$sample_ids%||%colnames(geno$dosage); sg <- rownames(geno$dosage)
  } else {
    if (!is.matrix(geno)) geno <- as.matrix(geno)
    gm <- geno; h1m <- h2m <- NULL; iids <- rownames(geno); sg <- colnames(geno)
  }
  snp_info$CHR <- .norm_chr_hap(snp_info$CHR)
  blocks$CHR   <- .norm_chr_hap(blocks$CHR)
  if (!is.null(chr)) { chr <- .norm_chr_hap(chr); blocks <- blocks[blocks$CHR%in%chr,] }
  res <- list()
  bi  <- data.frame(block_id=character(),CHR=character(),start_bp=integer(),
                    end_bp=integer(),n_snps=integer(),phased=logical(),
                    stringsAsFactors=FALSE)
  for (b in seq_len(nrow(blocks))) {
    blk <- blocks[b,]; cb <- as.character(blk$CHR)
    sb <- as.numeric(blk$start.bp); eb <- as.numeric(blk$end.bp)
    idx <- which(snp_info$CHR==cb&snp_info$POS>=sb&snp_info$POS<=eb)
    if (length(idx)<min_snps) next
    ci <- match(snp_info$SNP[idx],sg); ci <- ci[!is.na(ci)]
    if (length(ci)<min_snps) next
    bid <- paste0("block_",cb,"_",sb,"_",eb)
    if (isp) {
      hs <- vapply(seq_len(nrow(gm)), function(i) {
        v1 <- as.character(h1m[i,ci]); v1[is.na(h1m[i,ci])] <- na_char
        v2 <- as.character(h2m[i,ci]); v2[is.na(h2m[i,ci])] <- na_char
        paste0(paste(v1,collapse=""),"|",paste(v2,collapse=""))
      }, character(1L))
    } else {
      hs <- vapply(seq_len(nrow(gm)), function(i) {
        v <- as.character(gm[i,ci]); v[is.na(gm[i,ci])] <- na_char
        paste(v,collapse="")
      }, character(1L))
    }
    names(hs) <- iids; res[[bid]] <- hs
    bi <- rbind(bi, data.frame(block_id=bid,CHR=cb,start_bp=as.integer(sb),
                               end_bp=as.integer(eb),n_snps=length(ci),
                               phased=isp,stringsAsFactors=FALSE))
  }
  attr(res,"block_info") <- bi; res
}

#' Compute Haplotype Diversity Per Block
#' @description Calculates n_haplotypes, He, Shannon entropy, and
#'   freq_dominant per block. Phased data contributes two gamete observations
#'   per individual.
#' @param haplotypes List from extract_haplotypes().
#' @param missing_string Missing data marker. Default ".".
#' @return Data frame with block_id, CHR, start_bp, end_bp, n_snps, n_ind,
#'   n_haplotypes, He, Shannon, freq_dominant, phased.
#' @export
compute_haplotype_diversity <- function(haplotypes, missing_string=".") {
  bi <- attr(haplotypes,"block_info")
  rows <- lapply(seq_along(haplotypes), function(i) {
    bn <- names(haplotypes)[i]; hap <- haplotypes[[bn]]
    blk <- if(!is.null(bi)) bi[bi$block_id==bn,,drop=FALSE] else NULL
    phased <- isTRUE(blk$phased[1L])
    obs <- if(phased) unlist(strsplit(hap,"|",fixed=TRUE)) else hap
    obs <- obs[!grepl(missing_string,obs,fixed=TRUE)]
    ni  <- if(phased) length(obs)%/%2L else length(obs)
    if (!length(obs)) return(data.frame(block_id=bn,
                                        CHR=if(!is.null(blk))blk$CHR[1] else NA,
                                        start_bp=if(!is.null(blk))blk$start_bp[1] else NA,
                                        end_bp=if(!is.null(blk))blk$end_bp[1] else NA,
                                        n_snps=if(!is.null(blk))blk$n_snps[1] else NA,
                                        n_ind=0L,n_haplotypes=NA,He=NA,Shannon=NA,freq_dominant=NA,
                                        phased=phased,stringsAsFactors=FALSE))
    tbl <- table(obs); freq <- as.numeric(tbl)/sum(tbl)
    He <- 1-sum(freq^2); Sh <- -sum(freq*log(pmax(freq,.Machine$double.eps)))
    data.frame(block_id=bn,
               CHR=if(!is.null(blk))blk$CHR[1] else NA,
               start_bp=if(!is.null(blk))blk$start_bp[1] else NA,
               end_bp=if(!is.null(blk))blk$end_bp[1] else NA,
               n_snps=if(!is.null(blk))blk$n_snps[1] else NA,
               n_ind=ni,n_haplotypes=length(tbl),He=He,Shannon=Sh,freq_dominant=max(freq),
               phased=phased,stringsAsFactors=FALSE)
  })
  do.call(rbind,rows)
}

#' Map GWAS Hits to LD Blocks (Post-GWAS QTL Region Definition)
#' @description Maps significant GWAS markers onto LD blocks to define QTL
#'   regions. Blocks with significant markers from multiple traits are flagged
#'   pleiotropic. Implements the approach of Tong et al. (2024).
#' @param gwas_results Data frame: SNP, CHR, POS. Optional: P, trait.
#' @param blocks LD block data frame from run_Big_LD_all_chr().
#' @param snp_info Full SNP metadata data frame.
#' @param p_threshold Significance threshold. Default 5e-8. NULL = use all markers.
#' @param trait_col Trait column name. Default "trait".
#' @param min_snps Minimum block SNP count. Default 3L.
#' @return Data frame: block_id, CHR, start_bp, end_bp, n_snps_block,
#'   n_sig_markers, lead_snp, lead_p, traits, n_traits, pleiotropic.
#' @references Tong et al. (2024) Theor Appl Genet 137:274.
#' @export
define_qtl_regions <- function(gwas_results, blocks, snp_info,
                               p_threshold=5e-8, trait_col="trait", min_snps=3L) {
  # Accept "Marker" as an alias for "SNP" (OptSLDP / GWAS convention)
  if (!"SNP" %in% names(gwas_results) && "Marker" %in% names(gwas_results))
    gwas_results$SNP <- gwas_results$Marker
  miss <- setdiff(c("SNP","CHR","POS"),names(gwas_results))
  if (length(miss)) stop("gwas_results missing: ",paste(miss,collapse=","),call.=FALSE)
  gwas_results$CHR <- .norm_chr_hap(gwas_results$CHR)
  blocks$CHR <- .norm_chr_hap(blocks$CHR); snp_info$CHR <- .norm_chr_hap(snp_info$CHR)
  if (!is.null(p_threshold)&&"P"%in%names(gwas_results))
    gwas_results <- gwas_results[!is.na(gwas_results$P)&gwas_results$P<=p_threshold,]
  if (!nrow(gwas_results)){message("[define_qtl_regions] No significant markers.");return(data.frame())}
  if (!trait_col%in%names(gwas_results)) gwas_results[[trait_col]] <- "trait"
  rows <- list()
  for (b in seq_len(nrow(blocks))) {
    blk <- blocks[b,]; ch <- as.character(blk$CHR)
    sb <- as.numeric(blk$start.bp); eb <- as.numeric(blk$end.bp)
    hits <- gwas_results[gwas_results$CHR==ch&gwas_results$POS>=sb&gwas_results$POS<=eb,,drop=FALSE]
    if (!nrow(hits)) next
    nb <- sum(snp_info$CHR==ch&snp_info$POS>=sb&snp_info$POS<=eb)
    if (nb<min_snps) next
    traits <- unique(hits[[trait_col]])
    li <- if("P"%in%names(hits)) which.min(hits$P) else 1L
    rows[[length(rows)+1L]] <- data.frame(
      block_id=paste0("block_",ch,"_",sb,"_",eb), CHR=ch,
      start_bp=as.integer(sb), end_bp=as.integer(eb),
      n_snps_block=nb, n_sig_markers=nrow(hits),
      lead_snp=hits$SNP[li],
      lead_p=if("P"%in%names(hits))hits$P[li] else NA_real_,
      traits=paste(sort(traits),collapse=","), n_traits=length(traits),
      pleiotropic=length(traits)>1L, stringsAsFactors=FALSE)
  }
  if (!length(rows)){message("[define_qtl_regions] No overlapping blocks.");return(data.frame())}
  out <- do.call(rbind,rows); out[order(out$CHR,out$start_bp),]
}

#' Build Haplotype Dosage Matrix for Genomic Prediction
#' @description Converts haplotype strings to a numeric matrix for genomic
#'   prediction. Supports phased and unphased input with two encoding schemes.
#'
#'   encoding="additive_012" (default, recommended for GBLUP/rrBLUP/BGLR):
#'     Phased:   0=0 copies, 1=1 copy (het), 2=2 copies (hom)
#'     Unphased: 0=no match, 2=match (1 not identifiable without phase)
#'
#'   encoding="presence_02" (kernel methods, random forest):
#'     Phased:   2=either gamete matches, 0=neither, NA=missing
#'     Unphased: 2=match, 0=no match, NA=missing
#'
#' @param haplotypes List from extract_haplotypes().
#' @param top_n Integer or \code{NULL}. Maximum number of haplotype alleles
#'   to retain per block, ranked by frequency. \code{NULL} (default) retains
#'   all alleles that pass \code{min_freq} -- recommended for most analyses.
#'   Set an integer cap (e.g. \code{top_n = 5L}) only when you need to limit
#'   matrix width for memory reasons on panels with thousands of blocks and
#'   highly diverse haplotypes (many rare alleles above \code{min_freq}).
#' @param encoding "additive_012" (default) or "presence_02".
#' @param missing_string Missing data marker. Default ".".
#' @param scale_features Center and scale columns. Default FALSE.
#' @param min_freq Minimum allele frequency to include. Default 0.01.
#' @return Numeric matrix (individuals x haplotype allele columns).
#' @export
build_haplotype_feature_matrix <- function(haplotypes, top_n=NULL,
                                           encoding=c("additive_012","presence_02"),
                                           missing_string=".", scale_features=FALSE,
                                           min_freq=0.01) {
  encoding <- match.arg(encoding)
  bi       <- attr(haplotypes,"block_info")
  inames   <- names(haplotypes[[1L]])
  mats     <- vector("list",length(haplotypes))
  for (bk in seq_along(haplotypes)) {
    bn <- names(haplotypes)[bk]; hap <- haplotypes[[bn]]
    phased <- isTRUE(bi$phased[bi$block_id==bn][1L])
    if (phased) {
      parts <- strsplit(hap,"|",fixed=TRUE)
      g1 <- vapply(parts,`[`,character(1L),1L); g2 <- vapply(parts,`[`,character(1L),2L)
      miss <- grepl(missing_string,g1,fixed=TRUE)|grepl(missing_string,g2,fixed=TRUE)
      gam  <- c(g1[!grepl(missing_string,g1,fixed=TRUE)],g2[!grepl(missing_string,g2,fixed=TRUE)])
    } else {
      miss <- grepl(missing_string,hap,fixed=TRUE); gam <- hap[!miss]
    }
    tbl <- sort(table(gam),decreasing=TRUE); fv <- as.numeric(tbl)/sum(tbl)
    tbl  <- tbl[fv >= min_freq]
    tops <- if (is.null(top_n)) names(tbl) else names(tbl)[seq_len(min(as.integer(top_n), length(tbl)))]
    if (!length(tops)){mats[[bk]]<-NULL;next}
    mat <- matrix(NA_real_,nrow=length(inames),ncol=length(tops),
                  dimnames=list(inames,paste0(bn,"__hap",seq_along(tops))))
    for (j in seq_along(tops)) {
      ref <- tops[j]
      d <- if (identical(encoding,"additive_012")) {
        if (phased) ifelse(miss,NA_real_,(g1==ref)*1+(g2==ref)*1)
        else ifelse(miss,NA_real_,ifelse(hap==ref,2,0))
      } else {
        if (phased) ifelse(miss,NA_real_,ifelse(g1==ref|g2==ref,2,0))
        else ifelse(miss,NA_real_,ifelse(hap==ref,2,0))
      }
      mat[,j] <- d
    }
    mats[[bk]] <- mat
  }
  mats <- Filter(Negate(is.null),mats)
  if (!length(mats)) stop("No columns produced. Lower min_freq.",call.=FALSE)
  feat <- do.call(cbind,mats)
  if (scale_features) { feat <- scale(feat); feat[is.nan(feat)] <- 0 }
  feat
}

#' Write Haplotype Character (Nucleotide) Matrix
#'
#' Writes a matrix where each cell contains the nucleotide sequence of the
#' haplotype allele carried by each individual. Rows are haplotype alleles,
#' columns are individuals. This is the most interpretable format: you can
#' read directly which nucleotides define each haplotype allele and which
#' individuals carry it.
#'
#' @details
#' The cell value for individual i at haplotype allele h is:
#' \itemize{
#'   \item The nucleotide sequence (e.g. \code{"AGTTA"}) if the individual
#'     carries that allele (dosage = 2 for unphased, or present in either
#'     gamete for phased).
#'   \item \code{"-"} if the individual does not carry that allele.
#'   \item \code{"."} if the individual has missing data in that block.
#' }
#'
#' @param haplotypes List from \code{\link{extract_haplotypes}}.
#' @param snp_info   Data frame with \code{SNP}, \code{CHR}, \code{POS},
#'   \code{REF}, \code{ALT}.
#' @param out_file   Output file path (tab-delimited).
#' @param min_freq   Minimum haplotype frequency. Default \code{0.01}.
#' @param top_n      Integer or \code{NULL}. Cap alleles per block.
#'   \code{NULL} (default) keeps all above \code{min_freq}.
#' @param missing_string Missing genotype marker. Default \code{"."}.
#' @param verbose    Logical. Default \code{TRUE}.
#' @return Invisibly returns \code{out_file}.
#' @export
write_haplotype_character <- function(haplotypes, snp_info, out_file,
                                      min_freq       = 0.01,
                                      top_n          = NULL,
                                      missing_string = ".",
                                      verbose        = TRUE) {

  if (!all(c("CHR","POS","REF","ALT") %in% names(snp_info)))
    stop("snp_info must have columns CHR, POS, REF, ALT.", call. = FALSE)

  snp_info$CHR <- .norm_chr_hap(snp_info$CHR)
  bi           <- attr(haplotypes, "block_info")

  # Decode one dosage string to nucleotide sequence
  decode_str <- function(dstr, ref_v, alt_v) {
    chars <- strsplit(dstr, "", fixed = TRUE)[[1L]]
    n     <- min(length(chars), length(ref_v))
    vapply(seq_len(n), function(i) {
      switch(chars[i],
             "0" = ref_v[i],
             "2" = alt_v[i],
             "1" = paste0(ref_v[i], "/", alt_v[i]),
             "N")
    }, character(1L)) |> paste(collapse = "")
  }

  rows_all <- vector("list", length(haplotypes))

  for (bk in seq_along(haplotypes)) {
    bn   <- names(haplotypes)[bk]
    hap  <- haplotypes[[bn]]
    brow <- if (!is.null(bi)) bi[bi$block_id == bn, , drop = FALSE] else NULL

    # Parse CHR, start_bp from block name
    parts_bn <- strsplit(bn, "_", fixed = TRUE)[[1L]]
    chr  <- if (!is.null(brow) && nrow(brow)) brow$CHR[1L] else parts_bn[2L]
    sb   <- if (!is.null(brow) && nrow(brow)) brow$start_bp[1L] else
      suppressWarnings(as.integer(parts_bn[3L]))
    eb   <- if (!is.null(brow) && nrow(brow)) brow$end_bp[1L] else
      suppressWarnings(as.integer(parts_bn[4L]))

    # SNPs in block
    blk_snps <- snp_info[snp_info$CHR == chr &
                           snp_info$POS >= sb &
                           snp_info$POS <= eb, , drop = FALSE]
    blk_snps <- blk_snps[order(blk_snps$POS), ]
    if (!nrow(blk_snps)) next

    ref_v <- as.character(blk_snps$REF)
    alt_v <- as.character(blk_snps$ALT)

    # Alleles string for metadata column: REF1/ALT1;REF2/ALT2;...
    alleles_str <- paste(paste0(ref_v, "/", alt_v), collapse = ";")

    # Frequency table
    miss <- grepl(missing_string, hap, fixed = TRUE)
    gam  <- hap[!miss]
    if (!length(gam)) next
    tbl  <- sort(table(gam), decreasing = TRUE)
    fv   <- as.numeric(tbl) / sum(tbl)
    tbl  <- tbl[fv >= min_freq]
    tops <- if (is.null(top_n)) names(tbl)
    else names(tbl)[seq_len(min(as.integer(top_n), length(tbl)))]
    if (!length(tops)) next

    # For each top haplotype allele, build a row:
    # cols = hap_id | CHR | start_bp | end_bp | n_snps | Alleles | ind1 | ind2 | ...
    # cell = nucleotide_sequence if individual carries it, "-" if not, "." if missing
    for (r in seq_along(tops)) {
      dstr   <- tops[r]
      nuc_seq <- decode_str(dstr, ref_v, alt_v)
      hap_id  <- paste0(bn, "__hap", r)

      cell_vals <- ifelse(
        miss, missing_string,
        ifelse(hap == dstr, nuc_seq, "-")
      )

      # Build row: metadata cols + one col per individual
      ind_df <- as.data.frame(
        matrix(cell_vals, nrow = 1L,
               dimnames = list(NULL, names(hap))),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      row_df <- data.frame(
        hap_id   = hap_id,
        CHR      = chr,
        start_bp = as.integer(sb),
        end_bp   = as.integer(eb),
        n_snps   = nrow(blk_snps),
        Alleles  = alleles_str,
        ind_df,
        check.names = FALSE, stringsAsFactors = FALSE,
        row.names = NULL
      )
      rows_all[[length(rows_all) + 1L]] <- row_df
    }
  }

  if (!length(rows_all)) {
    message("[write_haplotype_character] No haplotypes passed filters. File not written.")
    return(invisible(out_file))
  }

  out <- do.call(rbind, rows_all)
  data.table::fwrite(out, out_file, sep = "\t", quote = FALSE, na = ".")
  if (verbose) message("[write_haplotype_character] ", out_file,
                       " (", nrow(out), " haplotypes x ",
                       ncol(out) - 6L, " individuals)")
  invisible(out_file)
}

#' Write Haplotype Feature Matrix as Numeric Dosage Table
#'
#' @description
#' Writes the haplotype dosage matrix in a tab-delimited format with
#' haplotype alleles as rows and individuals as columns. Metadata columns
#' (\code{hap_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#' \code{n_snps}, \code{alleles}, \code{frequency}) precede the individual
#' columns. Individual cells contain 0/1/2/NA dosage values.
#'
#' @param hap_matrix Numeric matrix (individuals x haplotype alleles) from
#'   \code{\link{build_haplotype_feature_matrix}}.
#' @param out_file   Output file path.
#' @param haplotypes List from \code{\link{extract_haplotypes}}. When supplied
#'   together with \code{snp_info}, the \code{alleles} and \code{frequency}
#'   metadata columns are populated.
#' @param snp_info   Data frame with \code{CHR}, \code{POS}, \code{REF},
#'   \code{ALT}. Required for \code{alleles} column.
#' @param sep        Field separator. Default \code{","}.
#' @param na_str     NA string. Default \code{"NA"}.
#' @param min_freq   Minimum frequency used when computing \code{alleles}.
#'   Default \code{0.01}.
#' @param missing_string Missing genotype marker. Default \code{"."}.
#' @param verbose    Logical. Default \code{TRUE}.
#' @return Invisibly returns \code{out_file}.
#' @export
write_haplotype_numeric <- function(hap_matrix, out_file,
                                    haplotypes = NULL,
                                    snp_info   = NULL,
                                    sep        = "\t",
                                    na_str     = "NA",
                                    min_freq   = 0.01,
                                    missing_string = ".",
                                    verbose    = TRUE) {
  # Output orientation: haplotype alleles as ROWS, individuals as COLUMNS.
  #
  # Metadata columns (before individual columns):
  #   hap_id        - haplotype allele identifier
  #   CHR           - chromosome
  #   start_bp      - block start position (bp)
  #   end_bp        - block end position (bp)
  #   n_snps        - number of SNPs in the block
  #   alleles       - nucleotide sequence of this specific haplotype allele,
  #                   decoded from the dosage string using SNP REF/ALT.
  #                   Unphased: gametic sequence (e.g. "AGTTA" for 5 SNPs).
  #                   Phased: gametic sequence of one strand (haploid).
  #                   Dosage = 2: both gametes carry this allele (unphased: present).
  #                   Dosage = 1: one gamete carries it (phased data only).
  #                   Dosage = 0: neither gamete carries this allele.
  #   frequency     - observed frequency in the panel
  #
  # Individual columns:
  #   Unphased: 0 = does not carry this allele, 2 = carries it, NA = missing.
  #   Phased:   0 = neither gamete carries it, 1 = one gamete carries it,
  #             2 = both gametes carry it, NA = missing.

  nhap    <- ncol(hap_matrix)
  hnm     <- colnames(hap_matrix)
  ind_ids <- rownames(hap_matrix)
  if (is.null(hnm))     hnm     <- paste0("hap", seq_len(nhap))
  if (is.null(ind_ids)) ind_ids <- paste0("Ind", seq_len(nrow(hap_matrix)))

  # Parse CHR, start_bp, end_bp from column names:
  # block_{CHR}_{start_bp}_{end_bp}__{hapN}
  parts    <- strsplit(hnm, "_", fixed = TRUE)
  chr_col  <- vapply(parts, function(x) if (length(x)>=2L) x[2L] else "NA", character(1L))
  spos_col <- suppressWarnings(
    as.integer(vapply(parts, function(x) if (length(x)>=3L) x[3L] else "0", character(1L))))
  epos_col <- suppressWarnings(
    as.integer(vapply(parts, function(x) if (length(x)>=4L) x[4L] else "0", character(1L))))
  spos_col[is.na(spos_col)] <- 0L
  epos_col[is.na(epos_col)] <- 0L

  # block_info for n_snps
  bi <- if (!is.null(haplotypes)) attr(haplotypes, "block_info") else NULL

  # Nucleotide decoder: converts a dosage string to a nucleotide sequence.
  # Works for both unphased ("02110") and phased gametic strings ("0110")
  # because extract_haplotypes() already splits phased strings into individual
  # gametes before tabling frequencies. So tops[] always contains single-strand
  # gametic dosage strings regardless of phasing.
  # "0" -> REF, "2" -> ALT, "1" -> REF/ALT (het, phased), "." -> N
  decode_str <- function(dstr, ref_v, alt_v) {
    chars <- strsplit(dstr, "", fixed = TRUE)[[1L]]
    n     <- min(length(chars), length(ref_v))
    paste(vapply(seq_len(n), function(i)
      switch(chars[i], "0"=ref_v[i], "2"=alt_v[i],
             "1"=paste0(ref_v[i],"/",alt_v[i]), "N"), character(1L)),
      collapse="")
  }

  # Compute metadata per haplotype allele
  alt_seq_col  <- rep(NA_character_, nhap)
  freq_col     <- rep(NA_real_, nhap)
  n_snps_col   <- rep(NA_integer_, nhap)

  if (!is.null(haplotypes) && !is.null(snp_info) &&
      all(c("CHR","POS","REF","ALT") %in% names(snp_info))) {
    snp_info$CHR <- .norm_chr_hap(snp_info$CHR)
    for (j in seq_len(nhap)) {
      bn <- sub("__hap[0-9]+$", "", hnm[j])
      hap_rank <- suppressWarnings(as.integer(sub(".*__hap","",hnm[j])))
      if (is.na(hap_rank)) next

      hap_strs <- haplotypes[[bn]]
      if (is.null(hap_strs)) next

      # Block SNPs
      blk_snps <- snp_info[snp_info$CHR == chr_col[j] &
                             snp_info$POS >= spos_col[j] &
                             snp_info$POS <= epos_col[j], , drop=FALSE]
      blk_snps <- blk_snps[order(blk_snps$POS), ]
      if (!nrow(blk_snps)) next

      ref_v <- as.character(blk_snps$REF)
      alt_v <- as.character(blk_snps$ALT)
      n_snps_col[j] <- nrow(blk_snps)

      # Frequency table of this block
      miss <- grepl(missing_string, hap_strs, fixed=TRUE)
      gam  <- hap_strs[!miss]
      if (!length(gam)) next
      tbl  <- sort(table(gam), decreasing=TRUE)
      fv   <- as.numeric(tbl)/sum(tbl)
      tbl  <- tbl[fv >= min_freq]
      tops <- names(tbl)

      # This specific haplotype allele (ALT sequence)
      if (hap_rank <= length(tops)) {
        alt_hapstr     <- tops[hap_rank]
        alt_seq_col[j] <- decode_str(alt_hapstr, ref_v, alt_v)  # stored as "alleles" column
        freq_col[j]    <- round(as.numeric(tbl[hap_rank])/sum(tbl), 4)
      }
    }
  }

  # Block n_snps from block_info if not computed above
  if (!is.null(bi)) {
    for (j in seq_len(nhap)) {
      if (!is.na(n_snps_col[j])) next
      bn <- sub("__hap[0-9]+$", "", hnm[j])
      brow <- bi[bi$block_id == bn, , drop=FALSE]
      if (nrow(brow)) n_snps_col[j] <- brow$n_snps[1L]
    }
  }

  # Transpose dosage matrix: haplotypes as rows
  t_mat <- t(hap_matrix)

  out <- data.frame(
    hap_id       = hnm,
    CHR          = chr_col,
    start_bp     = spos_col,
    end_bp       = epos_col,
    n_snps       = n_snps_col,
    alleles      = alt_seq_col,
    frequency    = freq_col,
    as.data.frame(t_mat, check.names=FALSE),
    check.names=FALSE, stringsAsFactors=FALSE, row.names=NULL
  )
  data.table::fwrite(out, out_file, sep=sep, quote=FALSE, na=na_str)
  if (verbose) message("[write_haplotype_numeric] ", out_file,
                       " (", nhap, " haplotypes x ", length(ind_ids), " individuals)")
  invisible(out_file)
}

#' Decode Haplotype Strings to Nucleotide Sequences
#'
#' Converts the dosage-encoded haplotype strings produced by
#' \code{\link{extract_haplotypes}} (e.g. \code{"02110"}) into
#' nucleotide sequences (e.g. \code{"AGTTА"}) using the REF and ALT
#' alleles of each SNP in the block.
#'
#' @details
#' Each character in a haplotype string is the dosage at one SNP in the block:
#' \itemize{
#'   \item \code{"0"} = homozygous REF  -> REF nucleotide (e.g. \code{A})
#'   \item \code{"1"} = heterozygous    -> REF/ALT (e.g. \code{A/G})
#'   \item \code{"2"} = homozygous ALT  -> ALT nucleotide (e.g. \code{G})
#'   \item \code{"."} = missing         -> \code{N}
#' }
#'
#' The result is a data frame with one row per unique haplotype allele per
#' block, showing its nucleotide sequence, frequency, and the REF/ALT at each
#' SNP position. This is the most interpretable representation of what each
#' haplotype allele actually encodes biologically.
#'
#' @param haplotypes List from \code{\link{extract_haplotypes}}.
#' @param snp_info   Data frame with columns \code{SNP}, \code{CHR},
#'   \code{POS}, \code{REF}, \code{ALT}. Must contain all SNPs in the blocks.
#' @param min_freq   Minimum haplotype frequency to include. Default \code{0.01}.
#' @param top_n      Integer or \code{NULL}. Maximum alleles per block.
#'   \code{NULL} (default) retains all above \code{min_freq}.
#' @param missing_string Missing genotype marker. Default \code{"."}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{block_id}{Block identifier.}
#'   \item{CHR}{Chromosome.}
#'   \item{start_bp, end_bp}{Block boundaries.}
#'   \item{hap_rank}{Rank by frequency (1 = most common).}
#'   \item{hap_id}{Column name as it appears in the feature matrix.}
#'   \item{dosage_string}{Raw dosage string e.g. \code{"02110"}.}
#'   \item{nucleotide_sequence}{Decoded nucleotide string e.g. \code{"AGTTА"}.}
#'   \item{frequency}{Observed frequency across non-missing individuals.}
#'   \item{n_carriers}{Number of individuals carrying this haplotype (dosage > 0).}
#'   \item{snp_positions}{Semicolon-separated CHR:POS of each SNP in the block.}
#'   \item{snp_alleles}{Semicolon-separated REF/ALT for each SNP.}
#' }
#'
#' @examples
#' data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
#' haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 3)
#' decoded <- decode_haplotype_strings(haps, ldx_snp_info)
#' head(decoded[, c("block_id","hap_rank","dosage_string",
#'                   "nucleotide_sequence","frequency")])
#'
#' @export
decode_haplotype_strings <- function(haplotypes, snp_info,
                                     min_freq       = 0.01,
                                     top_n          = NULL,
                                     missing_string = ".") {

  bi <- attr(haplotypes, "block_info")
  if (is.null(bi))
    stop("haplotypes must carry a block_info attribute from extract_haplotypes().",
         call. = FALSE)
  if (!all(c("CHR","POS","REF","ALT") %in% names(snp_info)))
    stop("snp_info must have columns CHR, POS, REF, ALT.", call. = FALSE)

  snp_info$CHR <- .norm_chr_hap(snp_info$CHR)
  rows_out     <- vector("list", length(haplotypes))

  for (bk in seq_along(haplotypes)) {

    bn  <- names(haplotypes)[bk]
    hap <- haplotypes[[bn]]
    brow <- bi[bi$block_id == bn, , drop = FALSE]
    if (!nrow(brow)) next

    chr <- brow$CHR[1L]
    sb  <- brow$start_bp[1L]
    eb  <- brow$end_bp[1L]

    # SNPs in this block from snp_info (ordered by position)
    blk_snps <- snp_info[snp_info$CHR == chr &
                           snp_info$POS >= sb &
                           snp_info$POS <= eb, , drop = FALSE]
    blk_snps <- blk_snps[order(blk_snps$POS), ]
    n_snps   <- nrow(blk_snps)

    if (!n_snps) next

    ref_v <- as.character(blk_snps$REF)
    alt_v <- as.character(blk_snps$ALT)
    pos_v <- blk_snps$POS

    # Frequency table of haplotype strings
    miss <- grepl(missing_string, hap, fixed = TRUE)
    gam  <- hap[!miss]
    if (!length(gam)) next

    tbl  <- sort(table(gam), decreasing = TRUE)
    fv   <- as.numeric(tbl) / sum(tbl)
    tbl  <- tbl[fv >= min_freq]
    if (!length(tbl)) next
    tops <- if (is.null(top_n)) names(tbl)
    else names(tbl)[seq_len(min(as.integer(top_n), length(tbl)))]

    # Decode each top haplotype string to nucleotides
    decode_one <- function(dstr) {
      chars <- strsplit(dstr, "", fixed = TRUE)[[1L]]
      # Pad or trim if lengths differ (safety guard)
      n <- min(length(chars), n_snps)
      nuc <- character(n)
      for (i in seq_len(n)) {
        nuc[i] <- switch(chars[i],
                         "0" = ref_v[i],
                         "2" = alt_v[i],
                         "1" = paste0(ref_v[i], "/", alt_v[i]),
                         "." = "N",
                         "N"
        )
      }
      paste(nuc, collapse = "")
    }

    n_carriers_fn <- function(dstr) sum(hap == dstr, na.rm = TRUE)

    for (r in seq_along(tops)) {
      dstr <- tops[r]
      rows_out[[length(rows_out) + 1L]] <- data.frame(
        block_id           = bn,
        CHR                = chr,
        start_bp           = as.integer(sb),
        end_bp             = as.integer(eb),
        hap_rank           = r,
        hap_id             = paste0(bn, "__hap", r),
        dosage_string      = dstr,
        nucleotide_sequence = decode_one(dstr),
        frequency          = round(as.numeric(tbl[dstr]) / sum(tbl), 4),
        n_carriers         = n_carriers_fn(dstr),
        snp_positions      = paste(paste0(chr, ":", pos_v), collapse = ";"),
        snp_alleles        = paste(paste0(ref_v, "/", alt_v), collapse = ";"),
        stringsAsFactors   = FALSE
      )
    }
  }

  if (!length(rows_out))
    return(data.frame())

  do.call(rbind, rows_out)
}


#' Write Haplotype Diversity Table
#' @param diversity Data frame from compute_haplotype_diversity().
#' @param out_file Output CSV path.
#' @param append_summary Append genome-wide mean row. Default TRUE.
#' @param verbose Logical. Default TRUE.
#' @return Invisibly returns out_file.
#' @export
write_haplotype_diversity <- function(diversity, out_file,
                                      append_summary=TRUE, verbose=TRUE) {
  if (append_summary) {
    sr <- data.frame(block_id="GENOME",CHR="ALL",start_bp=NA_integer_,
                     end_bp=NA_integer_,n_snps=NA_integer_,n_ind=NA_integer_,
                     n_haplotypes=round(mean(diversity$n_haplotypes,na.rm=TRUE),1),
                     He=mean(diversity$He,na.rm=TRUE),
                     Shannon=mean(diversity$Shannon,na.rm=TRUE),
                     freq_dominant=mean(diversity$freq_dominant,na.rm=TRUE),
                     phased=NA,stringsAsFactors=FALSE)
    diversity <- rbind(diversity,sr)
  }
  data.table::fwrite(diversity,out_file,sep=",",quote=FALSE,na="NA")
  if(verbose) message("[write_haplotype_diversity] ",out_file)
  invisible(out_file)
}
