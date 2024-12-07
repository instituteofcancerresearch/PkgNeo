#' Processes the mutation object
#'
#' @description This function is for the neoantigens project.
#' We call it with a mutation and sequence.
#'
#' @param mutation string
#' @param sequence string
#' @return An object with various info for mutation
#' @export
#' @examples
#' apply_mutation("ins2AA","atagtagtagtagct")
#' @references
#' \url{https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html}
#' @author Rachel Alcraft, \url{https://github.com/instituteofcancerresearch/PkgStopGain}


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
    
  return(list(
    type = type,
    start = pos_start,
    end = pos_end,
    offset = offset,
    modulo = modulo,
    ins = ins_piece,
    del = del_piece,
    new_seq = new_seq,
    ok = TRUE
  ))

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


