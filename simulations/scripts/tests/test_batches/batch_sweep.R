f = glue::glue

DATA = "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tests/batch_test"

B = c(1, 2, 5, 10, 20, 50, 100, 200)

res = map(B, function(x) {
    list(
        C = RcppCNPy::npyLoad(f("C_nb{x}.npy")),
        kdist = RcppCNPy::npyLoad(f("kdist_nb{x}.npy")),
        trep = RcppCNPy::npyLoad(f("trep_nb{x}.npy"))
    )
})


map(res, function(x) {
    list(
        mean_total = mean(rowSums(x$C)), ## average number of singletons per individual
        t0 = colSums(x$C) * x$trep
    )
})
