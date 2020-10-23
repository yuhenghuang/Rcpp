library(Rcpp)

sourceCpp("Modules.cpp")

unif <- new(Uniform, 0, 10)

unif$draw(10)

unif$min

loadModule("unif_module", TRUE)

example_int <- new (Example, 2L)

example_int$add()


example_int_int <- new(Example, 1L, 3L) 

example_int_int$add()


example_str <- new(Example, "abc", "de")

example_str$add()

example_str_int <- new(Example, "abc", 4L)

example_str_int$add()

example_iii <- new(Example, 1L, 4L, 3L)

example_iii$add()

example_def <- new(Example)

example_def$add()


setMethod("show", Example, function(object) {
  msg <- paste("x + y =", object$add())
  writeLines(msg)
})

Example$new(1L, 2L)
