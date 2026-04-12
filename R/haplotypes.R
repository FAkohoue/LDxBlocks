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

  # ── Input type detection ──────────────────────────────────────────────────
  # Three paths:
  #  1. LDxBlocks_backend  — STREAMING: one chromosome at a time, never full
  #     genome in RAM. Each chromosome is extracted, all its blocks processed,
  #     then freed with gc(FALSE) before the next chromosome.
  #  2. Phased list (hap1/hap2) — already in RAM, process as before.
  #  3. Numeric dosage matrix — already in RAM, process as before.
  is_backend <- inherits(geno, "LDxBlocks_backend")
  isp <- !is_backend && is.list(geno) && all(c("hap1","hap2")%in%names(geno))

  if (is_backend) {
    # ── STREAMING PATH: chromosome by chromosome ───────────────────────────
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

      # Extract ONLY this chromosome — RAM proportional to one chromosome
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
#' @param top_n Top haplotype alleles per block. Default 5L.
#' @param encoding "additive_012" (default) or "presence_02".
#' @param missing_string Missing data marker. Default ".".
#' @param scale_features Center and scale columns. Default FALSE.
#' @param min_freq Minimum allele frequency to include. Default 0.01.
#' @return Numeric matrix (individuals x haplotype allele columns).
#' @export
build_haplotype_feature_matrix <- function(haplotypes, top_n=5L,
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
    tbl <- tbl[fv>=min_freq]; tops <- names(tbl)[seq_len(min(top_n,length(tbl)))]
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
  if (!length(mats)) stop("No columns produced. Lower min_freq or top_n.",call.=FALSE)
  feat <- do.call(cbind,mats)
  if (scale_features) { feat <- scale(feat); feat[is.nan(feat)] <- 0 }
  feat
}

#' Write Haplotype Matrix as Numeric CSV
#' @description Rows=individuals, columns=haplotype alleles (0/1/2 or 0/2).
#'   Compatible with rrBLUP, BGLR, ASReml-R.
#' @param hap_matrix Numeric matrix (individuals x haplotype alleles).
#' @param out_file Output CSV path.
#' @param sep Separator. Default ",".
#' @param na_str NA string. Default "NA".
#' @param verbose Logical. Default TRUE.
#' @return Invisibly returns out_file.
#' @export
write_haplotype_numeric <- function(hap_matrix, out_file,
                                    sep="," ,na_str="NA", verbose=TRUE) {
  out <- data.frame(Sample=rownames(hap_matrix),
                    as.data.frame(hap_matrix,check.names=FALSE),
                    check.names=FALSE,stringsAsFactors=FALSE)
  data.table::fwrite(out,out_file,sep=sep,quote=FALSE,na=na_str)
  if(verbose) message("[write_haplotype_numeric] ",out_file,
                      " (",nrow(hap_matrix)," ind x ",ncol(hap_matrix)," cols)")
  invisible(out_file)
}

#' Write Haplotype Matrix in HapMap Format
#' @description Rows=haplotype alleles, columns=individuals.
#'   Encoding: 0->HH, 1->HA, 2->AA, NA->NN (H=ref token, A=alt token).
#'   Compatible with TASSEL and GAPIT.
#' @param hap_matrix Numeric matrix (individuals x haplotype alleles).
#' @param out_file Output path (.hmp.txt recommended).
#' @param ref_token Reference allele symbol. Default "H".
#' @param alt_token Alternate allele symbol. Default "A".
#' @param verbose Logical. Default TRUE.
#' @return Invisibly returns out_file.
#' @export
write_haplotype_hapmap <- function(hap_matrix, out_file,
                                   ref_token="H", alt_token="A", verbose=TRUE) {
  nhap    <- ncol(hap_matrix)
  hnm     <- colnames(hap_matrix)
  ind_ids <- rownames(hap_matrix)
  if (is.null(hnm))     hnm     <- paste0("hap", seq_len(nhap))
  if (is.null(ind_ids)) ind_ids <- paste0("Ind", seq_len(nrow(hap_matrix)))

  blks <- sub("__hap[0-9]+$","",hnm)
  rnks <- suppressWarnings(as.integer(sub(".*__hap","",hnm)))
  rnks[is.na(rnks)] <- seq_len(sum(is.na(rnks)))

  enc <- function(v,r,a) {
    o <- character(length(v))
    o[!is.na(v)&v==0] <- paste0(r,r); o[!is.na(v)&v==1] <- paste0(r,a)
    o[!is.na(v)&v==2] <- paste0(a,a); o[is.na(v)] <- "NN"; o
  }

  # Write line by line with explicit tab joins — avoids any fwrite quoting or
  # special-character issues with column names like "rs#" and "assembly#".
  header_line <- paste(c("rs#","alleles","chrom","pos","strand",
                         "assembly#","center","protLSID","assayLSID",
                         "panelLSID","QCcode", ind_ids), collapse="\t")

  data_lines <- vapply(seq_len(nhap), function(j) {
    gt <- enc(hap_matrix[, j], r=ref_token, a=alt_token)
    paste(c(hnm[j], paste0(ref_token,"/",alt_token), blks[j], rnks[j],
            "+","NA","NA","NA","NA","NA","NA", gt), collapse="\t")
  }, character(1L))

  writeLines(c(header_line, data_lines), con=out_file)
  if (verbose) message("[write_haplotype_hapmap] ",out_file,
                       " (",nhap," alleles x ",nrow(hap_matrix)," ind)")
  invisible(out_file)
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
