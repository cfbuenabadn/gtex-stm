library(dplyr)
library(tidyverse)
library(bedr)
library(readr)
library(ebpmf.alpha)
dir.create(tempdir())

merge_isoforms <- function(isoforms){

    isoform_list <- list(isoform_1 =  list(df = isoforms$isoform_1, factors = c('isoform_1')))
    i = 2
    
    for (name_iso in names(isoforms)){
    
        iso <- isoforms[[name_iso]]
        
        found_match <- FALSE
        for (name_iso2 in names(isoform_list)){
            iso2 <- isoform_list[[name_iso2]][['df']]
            # iso2 <- iso2$df
            if (found_match) {
                next
            } else {
            
                if (dim(iso)[1] == dim(iso2)[1]){
    
    
                    length_iso <- dim(iso)[1]
        
                    print(length_iso)
                    
                    exon_len_iso <- iso$end - iso$start
                    exon_len_iso2 <- iso2$end - iso2$start
                    
                    start_diff = (iso$start - iso2$start) %>% abs()
                    end_diff = (iso$end - iso2$end) %>% abs()
    

                    if (length_iso == 1) {
                        suma_start <- 0
                        suma_end <- 0
                    } else {
                        suma_start <- sum(start_diff[2:length_iso])
                        suma_end <- sum(end_diff[1:((length_iso)-1)])
                    }
        
                    first_exon_diff = start_diff[1]
                    last_exon_diff = end_diff[length_iso]
                    
                    frac_start_1 = first_exon_diff/(exon_len_iso[1])
                    frac_start_2 = first_exon_diff/(exon_len_iso2[1])
        
                    frac_end_1 = last_exon_diff/(exon_len_iso[length_iso])
                    frac_end_2 = last_exon_diff/(exon_len_iso2[length_iso])

                    first_exon_is_diff = ((first_exon_diff > 200) & ((frac_start_1 > 0.25) | (frac_start_2 > 0.25)))
                    last_exon_is_diff = ((last_exon_diff > 200) & ((frac_end_1 > 0.25) | (frac_end_2 > 0.25)))
    
               
                    if ((suma_start == 0) & (suma_end == 0) & (! first_exon_is_diff) & (! last_exon_is_diff)){
                        isoform_list[[name_iso2]]$factors <- c(isoform_list[[name_iso2]]$factors, name_iso)
                        found_match <- TRUE
                    } 
                
            }
        }
        }
        if (! found_match) {
    
            new_name_iso <- paste0('isoform_', i) 
            isoform_list <- append(isoform_list, list(list(df = iso, factors = c(name_iso))))
    
            names(isoform_list)[length(isoform_list)] <- new_name_iso
    
            i <- i+1
        }
    }
    return (isoform_list)
}


smooth_factor <- function(factor, cutoff_ = 0.25, smooth_fraction = 0.25, step_fraction = 10, exon_quant = 0.99, 
                          cutoff_strict = 0.1, #cutoff_strict = 0.25,#cutoff_strict = 0.1, 
                          pos_window = 100) {
  
  cutoff <- cutoff_
  factor <- factor / quantile(factor, exon_quant)
  step <- length(factor) / step_fraction
  
  factor_smooth <- numeric(length(factor))
    
  last_pos <- NULL
  
  for (idx in seq_along(factor)) {
    f <- factor[idx]
    
    if (f <= cutoff_strict) {
      factor_smooth[idx] <- 0
      last_pos <- NULL
      next
    }

    
    if (idx == 1) {
      if (f > cutoff) {
        f_smooth <- 1
        last_pos <- 0
      } else {
        f_smooth <- 0
        last_pos <- NULL
      }
    } else {
      if (factor_smooth[idx - 1] == 0) {
        if (f < cutoff) {
          f_smooth <- 0
          last_pos <- NULL
        } else {
          f_smooth <- 1
          last_pos <- idx
        }
      } else {
        last_pos <- max(last_pos, (idx - pos_window))
        f_exon_peak <- quantile(factor[last_pos:idx], exon_quant)
        
        if (f < (f_exon_peak * smooth_fraction) || f < cutoff_strict) {
            
          f_smooth <- 0
          cutoff <- max(f, cutoff_, cutoff_strict)
          last_pos <- NULL
        } else {
          f_smooth <- 1
        }
      }
    }
    
    factor_smooth[idx] <- f_smooth
  }

  first_smooth <- factor_smooth
  
  last_pos <- NULL
  
  for (idx in rev(seq_along(factor))) {
    f <- factor[idx]
    
    if (idx == length(factor)) {
      if (factor_smooth[idx] == 1) {
        last_pos <- length(factor) - 1
      } else {
        next
      }
    } else {
      previous_f <- factor_smooth[idx + 1]
      current_f <- factor_smooth[idx]
      
      if (previous_f == 0 && current_f == 1) {
        last_pos <- idx
        next
      } else if (previous_f == 1) {
        last_pos <- min(last_pos, (idx + pos_window))
        
        if (current_f == 0 && (last_pos - idx) < 10) {
          factor_smooth[idx:(last_pos + 1)] <- 0
          next
        }
        
        f_exon_peak <- quantile(factor[(idx + 1):(last_pos + 1)], exon_quant)
        
        if (f > (f_exon_peak * smooth_fraction) && f > cutoff_strict) {
          factor_smooth[idx] <- 1
        } else {
          next
        }
      }
    }
  }
  
  return(factor_smooth)
}


factor_lm <- function(factor, strand=NULL){

    print('HELLOWALLS')
    print(strand)
    x <- (1:length(factor))
    y <- as.vector(factor)
    y <- y/quantile(y,0.99)

    y_ <- y[y > 0.25]
    x_ <- x[y > 0.25]
    model <- lm(y_~x_) %>% summary
    intercept <- model$coefficients[1,1]
    coef <- model$coefficients[2,1]

    print(coef)

    if (is.null(strand)) {

        y_model <- (x*(coef))+intercept
        
        y_corrected <- y/pmax(y_model, 0.1)
        y_corrected <- ifelse(y > 0.1, y_corrected, y)

    } else {if (((strand == 'plus') & (coef > 0)) | ((strand=='minus') & (coef < 0))) {
        print('HELLO WALLS')
        print(coef)
        y_model <- (x*(coef))+intercept
        y_corrected <- y/pmax(y_model, 0.1)
        y_corrected <- ifelse(y > 0.1, y_corrected, y)

    } else {y_corrected = factor}

    return(y_corrected)}
}


correct_factor <- function(factor, junctions_bed, coords, correct_bias=TRUE, smooth_fraction = 0.25, strand=NULL){
    if (correct_bias) {y <- factor_lm(factor, strand)} else {y <- factor}
    y <- y/quantile(y, 0.99)
    binary_y <- smooth_factor(pmin(y, 1), exon_quant=1, smooth_fraction = smooth_fraction) #smooth_factor(y_corrected)
    segments_df <- find_continuous_segments(coords, binary_y)

    if ((dim(segments_df)[1] < 2) | (dim(junctions_bed)[1] == 0)) { 
        return(segments_df)
    } else {

    factor_junctions_bed <- get_factor_gaps(segments_df)
    bed_intersection <- bedr(input = list(a = factor_junctions_bed, b = junctions_bed), 
          method = "intersect", 
          params = "-wao"
          )

    corrected_exons <- get_corrected_exons(segments_df, factor_junctions_bed, bed_intersection, junctions_bed)

    print(segments_df)
    return(corrected_exons)
}
}


get_corrected_exons <- function(segments_df, factor_junctions_bed, bed_intersection, junctions_bed){
    corrected_exons <- data.frame(chrom = character(0), start = integer(0), end = integer(0))
    
    start_list <- c(segments_df$start[1])
    end_list <- c()
    for (i in 1:length(factor_junctions_bed$start)){
        # print("")
        # print(i)
        factor_junc <- factor_junctions_bed$names[i]
        best_junction <- bed_intersection %>% filter(names == factor_junc) %>% find_best_junction()
    
        if (is.null(best_junction)){
            next
        } else if (best_junction == 'gap') {
             end_list <- c(end_list, as.integer(factor_junctions_bed$start[i]))
            start_list <- c(start_list, as.integer(factor_junctions_bed$end[i]))
        } else{
        
        end_list <- c(end_list, (junctions_bed %>% filter(junc_names == best_junction) %>% pull(start)))
        start_list <- c(start_list, (junctions_bed %>% filter(junc_names == best_junction) %>% pull(end)))
        }
        }
    end_list <- c(end_list, segments_df$end[length(segments_df$end)])

    for (i in 1:length(start_list)){
        df_row <- data.frame(chrom = segments_df$chrom[1], start = start_list[i], end = end_list[i], stringsAsFactors = FALSE)
        corrected_exons <- rbind(corrected_exons, df_row)
        }
    return(corrected_exons)
}


get_isoforms <- function(EF, junctions_bed, coordinates, correct_bias=TRUE, smooth_fraction = 0.25, strand=NULL){
    K <- (EF %>% dim())[2]

    EF <- EF[2:((EF %>% dim())[1]-1),] 

    stopifnot(dim(EF)[1] == length(coordinates))
    
    isoforms <- list()
    for (k in 1:K){
        isoform_k <- correct_factor(EF[,k], junctions_bed, coordinates, correct_bias=correct_bias, smooth_fraction = smooth_fraction, strand=strand)
        isoforms[[paste0('isoform_', k)]] <- isoform_k
        }
    return(isoforms)
}

find_continuous_segments <- function(idx, x) {
  segments_df <- data.frame(chrom = character(0), start = integer(0), end = integer(0))
  in_segment <- FALSE
  start_idx <- NULL
  current_chrom <- NULL
  
  for (i in seq_along(x)) {
    if (x[i] == 1 && !in_segment) {
      in_segment <- TRUE
      start_idx <- strsplit(idx[i], ':')[[1]][2] %>% as.integer()
    } else if (x[i] == 0 && in_segment) {
      in_segment <- FALSE
      end_idx <- strsplit(idx[i - 1], ':')[[1]][2] %>% as.integer()
      segment <- data.frame(chrom = current_chrom, start = start_idx, end = end_idx, stringsAsFactors = FALSE)
      segments_df <- rbind(segments_df, segment)
    }
    
    if (in_segment && is.null(current_chrom)) {
      current_chrom <- strsplit(idx[i], ":")[[1]][1]
    }
  }
  
  # Check if the last segment is still continuing till the end
  if (in_segment) {
    end_idx <- strsplit(idx[length(idx)], ':')[[1]][2] %>% as.integer()
    segment <- data.frame(chrom = current_chrom, start = start_idx, end = end_idx, stringsAsFactors = FALSE)
    segments_df <- rbind(segments_df, segment)
  }
  
  return(segments_df)
}

get_factor_gaps <- function(df) {
  # Create a new dataframe with one less row
  new_df <- data.frame(chrom = df$chrom[-1],
                       start = df$end[-nrow(df)],
                       end = df$start[-1])

  new_df$names <- paste0('factor_junction', 1:dim(new_df)[1])
  factor_junctions_bed <- bedr::convert2bed(new_df)
  return(factor_junctions_bed)
}


find_best_junction <- function(bed_intersection_slice, gene_size=10000){
    gap_start <- as.integer(bed_intersection_slice$start)[1]
    gap_end <- as.integer(bed_intersection_slice$end)[1]
    gap_size <- (gap_end - gap_start)
    best_junc <- NULL
    any_close <- FALSE
    any_start_close <- FALSE
    any_end_close <- FALSE
    

    if (dim(bed_intersection_slice)[1] == 0) {
        if (gap_size >= gap_min_size){
            best_junc <- "gap"
            return(best_junc)
        } else {return (NULL)}
    } else if (bed_intersection_slice$V9[1] == '0'){
        gap_min_size <- gene_size/10
        if (gap_size >= gap_min_size){
            best_junc <- "gap"
        }
    } else {
        
        best_distance <- 1000
        best_percent <- 0
        
        overlapping_junctions <- bed_intersection_slice$junc_names %>% unique()
        for (junc in overlapping_junctions){
            bed_intersection_junc <- bed_intersection_slice %>% filter(junc_names == junc)
            junc_start <- bed_intersection_junc$start.b %>% as.integer()
            junc_end <- bed_intersection_junc$end.b %>% as.integer()
            junc_size <- bed_intersection_junc$V9 %>% as.integer()

            start_close <- abs(gap_start - junc_start)
            end_close <- abs(gap_end - junc_end)
            gap_percent <- junc_size/gap_size
            current_distance <- start_close + end_close

            # print(current_distance)

            # print(gap_percent)

            if (start_close <= 100) {any_start_close <- TRUE}
            if (end_close <= 100) {any_end_close <- TRUE}

            if ((start_close <= 100) || (end_close <= 100)) {any_close <- TRUE}

            if ((start_close <= 100) && (end_close <= 100) && (gap_percent >= 0.75)){

                if ((current_distance < best_distance) && (gap_percent > best_percent)){
                    best_junc <- junc
                    best_distance <- current_distance
                    best_percent <- gap_percent
                }
                }
            
            }
    }

    if (is.null(best_junc)){
            if ((gap_size > 1000) && (any_close)) {best_junc <- 'gap'}
            else if ((any_start_close) && (any_end_close)) {best_junc <- 'gap'}
        }
    
    return (best_junc)
    }


run_tabix <- function(coords, junc_file, gene_id){
    coords_end <- (coords[length(coords)] %>% str_split(':'))[[1]][2]
    coords_tabix <- paste(coords[1], coords_end, sep='-')
    
    junctions_bed <- tabix(coords_tabix, junc_file)

    if (is.null(junctions_bed)) {
        junctions_bed <- data.frame(chrom = character(0), start = integer(0), end = integer(0), junc_names = character(0))
        return(junctions_bed)
        }
        
    junctions_bed  <- junctions_bed %>% select(c('chrom', 'start', 'end', 'gene'))

    if (dim(junctions_bed)[1] == 0) { return(junctions_bed) }
    
    junctions_bed <- junctions_bed %>% filter(startsWith(gene, gene_id))

    if (dim(junctions_bed)[1] == 0) { return(junctions_bed) }
    
    junctions_bed$start <- junctions_bed %>% pull(start) %>% as.integer()
    junctions_bed$end <- junctions_bed %>% pull(end) %>% as.integer()
    
    junctions_bed$junc_names <- paste0('junction', 1:dim(junctions_bed)[1])

    junctions_bed <- junctions_bed %>% select(c('chrom', 'start', 'end', 'junc_names'))

    junctions_bed <- bedr::convert2bed(junctions_bed)

    return(junctions_bed)
}
