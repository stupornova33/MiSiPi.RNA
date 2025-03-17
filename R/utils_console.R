.overwrite_warning <- function(results_dir) {
  no_slash <- gsub("/", "", results_dir)
  inform <- glue::glue("It appears you have previous results stored in /{results_dir}")
  user_prompt <- "Is it okay to overwrite your old results?"
  choices <- c("No!" = "Save them!", "Sure!" = "Overwrite them!", "No way!" = "Keep them!")
  
  user_option <- utils::menu(
    choices = glue::glue(
      "{format(names(choices), justify = 'right')} --> {choices}"
    ),
    title = glue::glue(
      "\n{inform}\n",
      "{user_prompt} (0 to exit)"
    )
  )
  
  if (user_option == 0) {
    cli::cli_inform("Stopping at user request. Come back anytime!")
    # Store old settings, change them, then change them back when exiting
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  } else if (user_option == 1 | user_option == 3) {
    backup_dir <- paste0(no_slash, "_", format(lubridate::now(), "%Y%m%d_%H%M%S/"))
    cli::cli_inform(c(
      "Backing up old results to the following location:",
      backup_dir,
      ""))
    file.rename(results_dir, backup_dir)
  } else if (user_option == 2) {
    cli::cli_inform(c(
      "Overwriting old results. Writing new results to:",
      glue::glue("{results_dir}/"),
      ""))
    unlink(results_dir, recursive = TRUE)
  } else {
    cli::cli_warn(glue::glue("Unknown option {user_option} passed to backup prompt."))
    backup_dir <- paste0(no_slash, "_", format(lubridate::now(), "%Y%m%d_%H%M%S/"))
    cli::cli_inform(c(
      "Backing up old results just in case to the following location:",
      backup_dir,
      ""))
    file.rename(results_dir, backup_dir)
  }
}

.inform_complete <- function(results_dir) {
  cli::cli_inform(c(
    "",
    "Processing complete!",
    glue::glue("Results are stored in: --> {results_dir}/")
  ))
}

.inform_iteration <- function(i, i_max, chrom_name, strand = NULL) {
  if (is.null(strand)) {
    msg <- glue::glue("{i} out of {i_max} | {chrom_name}")
  } else {
    msg <- glue::glue("{i} out of {i_max} | {chrom_name} | {strand} strand")
  }
  cli::cli_inform(msg)
}