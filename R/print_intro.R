# If not seen as to annoying it will be called at the start misipi_rna()
.print_intro <- function(roi, bam, genome, method) {
  # Widest line here is 81
  cat("\n")
  cat(r"(__/\\\\____________/\\\\___________/\\\\\\\\\\\__________/\\\\\\\\\\\\\_________)", "\n")
  cat(r"( _\/\\\\\\________/\\\\\\_________/\\\/////////\\\_______\/\\\/////////\\\_______)", "\n")
  cat(r"(  _\/\\\//\\\____/\\\//\\\__/\\\__\//\\\______\///___/\\\_\/\\\_______\/\\\__/\\\_)", "\n")
  cat(r"(   _\/\\\\///\\\/\\\/_\/\\\_\///____\////\\\_________\///__\/\\\\\\\\\\\\\/__\///__)", "\n")
  cat(r"(    _\/\\\__\///\\\/___\/\\\__/\\\______\////\\\_______/\\\_\/\\\/////////_____/\\\_)", "\n")
  cat(r"(     _\/\\\____\///_____\/\\\_\/\\\_________\////\\\___\/\\\_\/\\\_____________\/\\\_)", "\n")
  cat(r"(      _\/\\\_____________\/\\\_\/\\\__/\\\______\//\\\__\/\\\_\/\\\_____________\/\\\_)", "\n")
  cat(r"(       _\/\\\_____________\/\\\_\/\\\_\///\\\\\\\\\\\/___\/\\\_\/\\\_____________\/\\\_)", "\n")
  cat(r"(        _\///______________\///__\///____\///////////_____\///__\///______________\///__)", "\n")
  cat("\n")

  msg1 <- paste("Processing method:", method)
  msg2 <- paste("BED:", vars$roi, "| BAM:", vars$bam_file, "| GENOME:", vars$genome)
  
  cli::cli_inform(c(msg1, msg2, ""))
}
