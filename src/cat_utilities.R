# Cat utilities ----
# ||||||||||||||||||

cat.script_title <- function(title) {
  len <- nchar(title)
  spacer <- strrep("%", len)
  
  if(VERBOSE){
    cat(paste("\n%%%", spacer, "%%%", sep = ""))
    cat(paste("\n%% ", title,  " %%", sep = ""))
    cat(paste("\n%%%", spacer, "%%%\n\n", sep = ""))
  }
}

cat.section_title <- function(title) {
  len <- nchar(title)
  spacer <- strrep("|", len)
  
  if(VERBOSE){
    cat(paste("\n|||", spacer, "|||", sep = ""))
    cat(paste("\n|| ", title,  " ||", sep = ""))
    cat(paste("\n|||", spacer, "|||\n\n", sep = ""))
  }
}

cat.subsection_title <- function(title) {
  if(VERBOSE){
    cat(paste("\n- ", title,  "\n\n", sep = ""))
  }
}