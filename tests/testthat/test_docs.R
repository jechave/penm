# The help pages, as a user reads them.
#
# Roxygen comments are written in markdown. If roxygen's markdown support is
# off, the markdown is copied into the help pages literally, and R CMD check
# does not notice: to it, it is text. Here every page is rendered to text, the
# way ?topic shows it, and searched for markdown that was not converted.

# the man pages: from the source tree when it is there (devtools::test()), from
# the installed package otherwise (R CMD check)
rendered_help <- function() {
  man_dir <- test_path("..", "..", "man")
  if (dir.exists(man_dir)) {
    files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
    rds <- lapply(files, tools::parse_Rd)
    names(rds) <- basename(files)
  } else {
    rds <- tools::Rd_db("penm")
  }
  vapply(rds, function(rd) {
    paste(utils::capture.output(tools::Rd2txt(rd, options = list(underline_titles = FALSE))),
          collapse = "\n")
  }, character(1))
}

test_that("no help page shows markdown or LaTeX literally", {
  pages <- rendered_help()
  expect_gt(length(pages), 20)

  literal <- c(
    backtick = "`",
    bold     = "\\*\\*",
    link     = "\\[[A-Za-z_.]+(\\(\\))?\\]",
    # a formula without a plain-text form shows its LaTeX in text help
    latex    = "\\\\(mathbf|delta|sigma|lambda|sqrt|frac|sum)\\b"
  )
  for (kind in names(literal)) {
    offending <- names(pages)[grepl(literal[[kind]], pages)]
    expect(length(offending) == 0,
           paste0("literal ", kind, " in: ", paste(offending, collapse = ", ")))
  }
})
