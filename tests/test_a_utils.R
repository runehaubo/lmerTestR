# test_a_utils.R

library(lmerTest)

# test safeDeparse() - equivalence and differences to deparse():
deparse_args <- formals(deparse)
safeDeparse_args <- formals(lmerTest:::safeDeparse)
stopifnot(
  all.equal(names(deparse_args), names(safeDeparse_args)),
  all.equal(deparse_args[!names(deparse_args) %in% "width.cutoff"],
            safeDeparse_args[!names(safeDeparse_args) %in% "width.cutoff"]),
  all.equal(deparse_args[["width.cutoff"]], 60L),
  all.equal(safeDeparse_args[["width.cutoff"]], 500L)
)

