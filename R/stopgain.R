#' Return the upstream and downstream stop gains for a mutation
#'
#' @description This function is for the neoantigens project.
#' We call it with a mutation, gene and sequence.
#'
#' @param mutation string
#' @param gene string
#' @param sequence string
#' @return An object with various info for the stop gain
#' @export
#' @examples
#' stop_gain_factors("ins2AA","atagtagtagtagct")
#' @references
#' \url{https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html}
#' @author Rachel Alcraft, \url{https://github.com/instituteofcancerresearch/PkgStopGain}

stop_gain_factors <- function(mutation, sequence) {  
  aa_orig <- get_protein_from_dna(sequence)
  mut_pos <- apply_mutation(mutation, sequence)  
  if (mut_pos$ok == FALSE) {
    return(list(ok = FALSE, msg = mut_pos$msg))
  } else {
    down_stop_gain <- get_down_stop_gain(sequence, mut_pos$new_seq, aa_orig, 0)
    up_stop_gain <- get_up_stop_gain(sequence, mut_pos$new_seq,aa_orig, mut_pos$offset)
    return(list(
      # from the mutation
      ok = TRUE,
      dna_orig = sequence,
      aa_orig = aa_orig,
      type = mut_pos$type,
      new_dna = mut_pos$new_seq,
      # from the stop-gain
      down_all = down_stop_gain$aa_all,
      down_stop_gain = down_stop_gain$aa_diff,
      down_dna_stop_gain = down_stop_gain$dna_diff,
      up_all = up_stop_gain$aa_all,
      up_stop_gain = up_stop_gain$aa_diff,
      up_dna_stop_gain = up_stop_gain$dna_diff,
      version=1.00
    ))
  }
}

get_protein_from_dna <- function(dna, shift = 0) {  
  protein <- ""
  vals <- seq(1 + shift, nchar(dna), by = 3)
  for (i in vals) {
    codon <- substring(dna, i,i + 2)
    aa <- get_aa_from_codon(codon)    
    protein <- paste(protein, aa, sep = "")
  }
  return(protein)
}

get_aa_from_codon <- function(codon) {
  codon <- paste(",", toupper(codon), ",", sep = "")
  if (grepl(codon,",TAA,TAG,TGA,")) {
    return("*")
  } else if (grepl(codon, ",GCT,GCC,GCA,GCG,")) {
    return("A")
  } else if (grepl(codon, ",TGT,TGC,")) {
    return("C")
  } else if (grepl(codon, ",GAT,GAC,")) {
    return("D")
  } else if (grepl(codon, ",GAA,GAG,")) {
    return("E")
  } else if (grepl(codon, ",TTT,TTC,")) {
    return("F")
  } else if (grepl(codon, ",GGT,GGC,GGA,GGG,")) {
    return("G")
  } else if (grepl(codon, ",CAT,CAC,")) {
    return("H")
  } else if (grepl(codon, ",ATT,ATC,ATA,")) {
    return("I")
  } else if (grepl(codon, ",AAA,AAG,")) {
    return("K")
  } else if (grepl(codon, ",TTA,TTG,CTT,CTC,CTA,CTG,")) {
    return("L")
  } else if (grepl(codon, ",ATG,")) {
    return("M")
  } else if (grepl(codon, ",AAT,AAC,")) {
    return("N")
  } else if (grepl(codon, ",CCT,CCC,CCA,CCG,")) {
    return("P")
  } else if (grepl(codon, ",CAA,CAG,")) {
    return("Q")
  } else if (grepl(codon, ",CGT,CGC,CGA,CGG,AGA,AGG,")) {
    return("R")
  } else if (grepl(codon, ",TCT,TCC,TCA,TCG,AGT,AGC,")) {
    return("S")
  } else if (grepl(codon, ",ACT,ACC,ACA,ACG,")) {
    return("T")
  } else if (grepl(codon, ",GTT,GTC,GTA,GTG,")) {
    return("V")
  } else if (grepl(codon, ",TGG,")) {
    return("W")
  } else if (grepl(codon, ",TAT,TAC,")) {
    return("Y")
  }
}

get_down_stop_gain <- function(dna_seq, new_seq, aa_seq, shift = 0) {
  # The downstream stop gain is first stop codon after the mutation
  # (downstream of the mutation, to the right)
  # The sequence will shift from the wildtype downstream until there is a stop
  # WWWWWWWWWXXXXXXXXXX*YYYYYYYYYYYYY
  # where W is the wiltype sequence,
  # X is shifted sequence after the mutation,
  # Y is the shifted sequence after the stop codon
  # We want to return the xxxxxxxxx*
  # And the equivalent dna sequence
  aa_with_down_shift <- get_protein_from_dna(new_seq, shift = shift)
  min_len <- min(nchar(aa_with_down_shift), nchar(aa_seq))  
  char_diffs <- ""
  dna_diffs <- ""
  has_changed <- FALSE
  for (i in 1:min_len) {
    char_orig <- substring(aa_seq, i, i)
    char_new <- substring(aa_with_down_shift, i, i)
    codon <- substring(new_seq, i*3 - 2, i*3)        
    if (char_orig != char_new || has_changed) {
      has_changed <- TRUE
      char_diffs <- paste(char_diffs, char_new, sep = "")
      dna_diffs <- paste(dna_diffs, codon, sep = "")
      if (char_new == "*") {
        break
      }
    }
  }
  return(list(
    aa_diff = char_diffs,
    dna_diff = dna_diffs,
    aa_all = aa_with_down_shift
  ))
}

get_up_stop_gain <- function(dna_seq, new_seq, aa_seq, shift = 0) {
  # The upstream stop gain is first stop codon before the mutation
  # (upstream of the mutation, to the left)
  # This assumes that the protein sequence is read from right to left
  # The sequence will shift from the wildtype upstream until there is a stop
  # YYYYYYYYY*XXXXXXXXXXWWWWWWWWWW
  # where W is the wiltype sequence,
  # X is shifted sequence after the mutation (right to left),
  # Y is the shifted sequence after the stop codon (right to left)
  # We want to returnn the *xxxxxxxxx
  # And the equivalent dna sequence
  aa_with_up_shift <- get_protein_from_dna(new_seq, shift = shift)
  min_len <- min(nchar(aa_with_up_shift), nchar(aa_seq))
  char_diffs <- ""
  dna_diffs <- ""
  has_changed <- FALSE
  # reverse backwards down the sequences
  rev_aa_seq <- stringi::stri_reverse(aa_seq)
  rev_aa_with_up_shift <- stringi::stri_reverse(aa_with_up_shift)  
  codon_pos <- nchar(new_seq) - 2
  for (i in 1:min_len) {
    char_orig <- substring(rev_aa_seq, i, i)
    char_new <- substring(rev_aa_with_up_shift, i, i)
    codon <- substring(new_seq, codon_pos, codon_pos + 2)
    codon_pos <- codon_pos - 3
    #print(paste(i, char_orig, char_new, codon_pos, codon_pos + 2, codon)) DEBUG
    if (char_orig != char_new || has_changed) {
      has_changed <- TRUE
      char_diffs <- paste(char_new, char_diffs, sep = "")
      dna_diffs <- paste(codon, dna_diffs, sep = "")
      if (char_new == "*") {
        break
      }
    }
  }
  return(list(
    aa_diff = char_diffs,
    dna_diff = dna_diffs,
    aa_all = aa_with_up_shift
   ))
}