library(Rcpp)

sourceCpp("Modules.cpp")

unif <- new(Uniform, 0, 10)

unif$draw(10)

unif$min

loadModule("unif_module", TRUE)


example_int <- new(Example, 1L, 3L) 

example_int$add()


example_str <- new(Example, "abc", "de")

example_str$add()
