#' Return the upstream and downstream neoantigen information along with distance and sequence
#'
#' @description This function is for the neoantigens project.
#' We call it with a mutation, gene and sequence.
#'
#' @param mutation string
#' @param gene string
#' @param sequence string
#' @return An object with various info
#' @export
#' @examples
#' ff_rtest_add_basic(1, 1)
#' ff_rtest_add_basic(c(10,10), c(1,3))
#' @references
#' \url{https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html}
#' @author Rachel Alcraft, \url{https://github.com/instituteofcancerresearch/r-reversions}

up_down_factors <- function(mutation, sequence) {  
  aa_orig <- get_protein_from_dna(sequence)
  mut_pos <- apply_mutation(mutation, sequence)
  if (mut_pos$ok == FALSE) {
    return(list(ok = FALSE, msg = mut_pos$msg))
  } else {    
    distance <- mut_pos$distance    
    return(list(
      ok = TRUE,
      dna_orig = sequence,
      aa_orig = aa_orig,
      type = mut_pos$type,      
      distance = distance,
      ins_piece = mut_pos$ins,
      del_piece = mut_pos$del,
      offset = mut_pos$offset,
      new_dna = mut_pos$new_seq,
      down_aa = mut_pos$down_aa,
      up_aa = mut_pos$up_aa,
      up_stop = mut_pos$up_stop,
      down_stop = mut_pos$down_stop,
      up_stub = mut_pos$up_stub,
      down_stub = mut_pos$down_stub,
      down_non_wild = mut_pos$down_non_wild,
      len_orig = nchar(sequence)
    ))
  }
}

########################################################################################
apply_mutation <- function(mutation, sequence) {          

  # A bit of basic string prep
  mutation <- gsub(" ", "", mutation)
  mutation <- gsub("c.", "", mutation)
  sequence <- tolower(sequence)

  if (grepl("delins",mutation)) {
    mut_bits <- get_delins_mut_bits(mutation, sequence)
  } else if (grepl("ins", mutation)) {
    mut_bits <- get_ins_mut_bits(mutation, sequence)
  } else if (grepl("del", mutation)) {
    mut_bits <- get_del_mut_bits(mutation, sequence)
  } else if (grepl("dup", mutation)) {
    mut_bits <- get_dup_mut_bits(mutation, sequence)
  } else {
    return(list(start = 0, end = 0, ok = FALSE, msg = "Not handled mut type"))
  }

  # Preparing the variables we need to apply the mutation
  type <- mut_bits$type
  del_piece <- mut_bits$del
  ins_piece <- mut_bits$ins
  pos_start <- mut_bits$start
  pos_end <- mut_bits$end
  offset <- mut_bits$offset  
  modulo <- offset %% 3
  if (modulo == 0) {
    return(list(type = type, start = -1, end = -1, ok = FALSE, 
    msg = "Mutation does not cause frameshift"))
      }
  
  new_seq <- paste(
                    substring(sequence, 1, pos_start-1),
                    toupper(ins_piece),
                    substring(sequence, pos_end+1, nchar(sequence)),sep= "")
  
  aa_with_up_shift <- get_protein_from_dna(new_seq, shift = modulo)
  aa_with_down_shift <- get_protein_from_dna(new_seq, shift = 0)
    
  # We want the last stop in the upstream shift
  has_stop <- FALSE
  last_stop <- nchar(aa_with_up_shift)
  for (i in seq(nchar(aa_with_up_shift), by=-1)) {
    if (substring(aa_with_up_shift,i,i) == "*") {
      last_stop <- i
      has_stop <- TRUE
      break
    }
  }
  if (!has_stop) {
    last_stop <- 1
  }
  # up_stub <- substring(aa_with_up_shift, last_stop,nchar(aa_with_up_shift))
  

  # we want the first stopin the downstream shift
  has_stop <- FALSE
  first_stop <- 1
  for (i in seq(1, nchar(aa_with_down_shift))) {
    if (substring(aa_with_down_shift, i, i) == "*") {
      first_stop <- i
      has_stop <- TRUE
      break
    }
  }
  if (!has_stop) {
    first_stop <- nchar(aa_with_down_shift)
  }
  # down_stub <- substring(aa_with_down_shift, 1, first_stop)
  
  up_plus_10 <- substring(aa_with_up_shift, last_stop, (pos_start + offset)/3 + 10)
  down_plus_10 <- substring(aa_with_down_shift, pos_start/3 - 10, first_stop)

  # Downstream diffs, start at the beginning of original sequence and downstream and start adding when there is a diff
  pos_down <- 0
  for (i in seq(1, nchar(sequence))) {
    if (substring(sequence, i, i) != substring(aa_with_down_shift, i, i)) {
      pos_down <- i
      break
    }
  }

  down_non_wild <- substring(aa_with_down_shift, pos_down, nchar(aa_with_down_shift))
  
  

      
  distance <- first_stop - last_stop
  return(list(
    type = type,
    start = pos_start,
    end = pos_end,
    offset = offset,
    modulo = modulo,
    ins = ins_piece,
    del = del_piece,
    new_seq = new_seq,
    down_aa = aa_with_down_shift,
    up_aa = aa_with_up_shift,
    up_stop = last_stop,
    down_stop = first_stop,
    up_plus_10 = up_plus_10,
    down_plus_10 = down_plus_10,
    down_non_wild = down_non_wild,
    distance = distance,
    ok = TRUE
  ))

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


get_delins_mut_bits <- function(mutation, sequence) {
  type <- "delins"
  mut_pieces <- unlist(strsplit(mutation, "delins"))
  ins_piece <- mut_pieces[2]
  if (grepl("[0-9]", ins_piece)) {
    return(list(
                type = "delins",
                start = -1,
                end = -1,
                ok = FALSE,
                msg = "Contains numbers not dna"))
  }
  del_pieces <- mut_pieces[1]
  if (grepl("_", del_pieces)) {
    del_pieces <- unlist(strsplit(del_pieces, "_"))
    pos_start <- as.numeric(del_pieces[1])
    pos_end <- as.numeric(del_pieces[2])
  } else {
    pos_start <- as.numeric(del_pieces)
    pos_end <- pos_start + 1
  }
  del_piece <- substring(sequence, pos_start, pos_end)
  offset <- nchar(ins_piece) - nchar(del_piece)
  return(list(
    type = type,
    start = pos_start,
    end = pos_end,
    offset = offset,
    ins = ins_piece,
    del = del_piece
  ))
}

get_del_mut_bits <- function(mutation, sequence) {
  type <- "del"
  mut_pieces <- unlist(strsplit(mutation, "del"))
  pos_pieces <- mut_pieces[1]
  if (grepl("_", pos_pieces)) {
    pos_pieces <- unlist(strsplit(pos_pieces, "_"))
    pos_start <- as.numeric(pos_pieces[1])
    pos_end <- as.numeric(pos_pieces[2])
  } else {
    pos_start <- as.numeric(pos_pieces)
    pos_end <- pos_start
  }
  del_piece <- substring(sequence, pos_start, pos_end)
  offset <- - nchar(del_piece)
  return(list(
    type = type,
    start = pos_start,
    end = pos_end,
    offset = offset,
    ins = "",
    del = del_piece
  ))
}

get_ins_mut_bits <- function(mutation,sequence) {
  type <- "ins"
  mut_pieces <- unlist(strsplit(mutation, "ins"))  
  ins_piece <- mut_pieces[2]
  if (grepl("[0-9]", ins_piece)) {
    return(list(type = "ins", start = -1, end = -1, ok = FALSE,
                msg = "Contains numbers not dna"))
  }
  pos_pieces <- mut_pieces[1]
  if (grepl("_", pos_pieces)) {
    pos_pieces <- unlist(strsplit(pos_pieces, "_"))
    pos_start <- as.numeric(pos_pieces[1])
    pos_end <- as.numeric(pos_pieces[2])
  } else {
    pos_start <- as.numeric(pos_pieces)
    pos_end <- pos_start + 1
  }
  
  if (pos_end != pos_start + 1) {
    return(list(type = "ins", start = -1, end = -1, ok = FALSE,
                msg = "Insertion needs to be between positions"))
  }
  pos_end <- pos_start + 1
  offset <- nchar(ins_piece)
  return(list(
    type = type,
    start = pos_start + 1,
    end = pos_end - 1,
    offset = offset,
    ins = ins_piece,
    del = ""
  ))
  
}

get_dup_mut_bits <- function(mutation,sequence) {
  type <- "dup"
  mut_pieces <- unlist(strsplit(mutation, "dup"))
  pos_pieces <- mut_pieces[1]
  if (grepl("_", pos_pieces)) {
    pos_pieces <- unlist(strsplit(pos_pieces, "_"))
    pos_start <- as.numeric(pos_pieces[1])
    pos_end <- as.numeric(pos_pieces[2])
  } else {
    pos_start <- as.numeric(pos_pieces)
    pos_end <- pos_start
  }
  ins_piece <- toupper(substring(sequence, pos_start, pos_end))  
  pos_end <- pos_start + 1
  offset <- nchar(ins_piece)
  return(list(
    type = type,
    start = pos_start + 1,
    end = pos_end - 1,
    offset = offset,
    ins = ins_piece,
    del = ""
  ))
}