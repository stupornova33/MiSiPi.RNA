.log <- function(msg, wkdir, level = c("info", "warn", "error", "debug"), func = NULL, include_timestamp = FALSE, add_newline = FALSE) {
  if (is.null(msg)) {
    return(NULL)
  }
  
  level <- match.arg(level)
  
  # Adds an additional new line creating a blank line in the log file
  if (add_newline) {
    msg <- paste0(msg, "\n")
  }
  
  # If func (calling function) is provided, add it to front of msg
  if (!is.null(func)) {
    msg_prefix <- paste0()
    
    msg <- paste0("Function: ", func, ": ", msg)
  }
  
  # Add timestamp to msg if provided
  if (include_timestamp) {
    now <- lubridate::now()
    
    msg <- paste0("[", now, "] ", msg)
  }
  
  # In the case of warnings and errors, append their tag to the front of the message
  log_line_prefix <- switch(
    level,
    "info" = "",
    "warn" = "[WARNING]",
    "error" = "[ERROR]",
    "debug" = "[DEBUG]"
  )
  
  if (level != "info") {
    msg <- paste(log_line_prefix, msg)
  }
  
  filepath <- file.path(wkdir, "log.txt")
  
  con <- file(filepath, "a")
  writeLines(msg, con)
  close(con)
}
