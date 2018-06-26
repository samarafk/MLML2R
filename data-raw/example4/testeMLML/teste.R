library(MLML2R)

data(MethylatedBS_sim)
data(MethylatedOxBS_sim)
data(UnMethylatedBS_sim)
data(UnMethylatedOxBS_sim)
data(MethylatedTAB_sim)
data(UnMethylatedTAB_sim)

source("MLML2R.R")

MethylatedBS_sim[1:10,1] <- NA
MethylatedBS_sim[11:20,1] <- 0
UnMethylatedBS_sim[11:20,1] <- 0
UnMethylatedOxBS_sim[11:20,1] <-0
MethylatedOxBS_sim[11:20,1] <-0
UnMethylatedTAB_sim[11:20,1] <-0
MethylatedTAB_sim[11:20,1] <- 0

results_em <- MLML(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                   L.matrix = UnMethylatedOxBS_sim, M.matrix = MethylatedOxBS_sim,iterative=TRUE)
results_exact <- MLML(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                      L.matrix = UnMethylatedOxBS_sim, M.matrix = MethylatedOxBS_sim)

results_em1 <- MLML2(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                   L.matrix = UnMethylatedOxBS_sim, M.matrix = MethylatedOxBS_sim,iterative=TRUE)



results_em <- MLML(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                   L.matrix = UnMethylatedOxBS_sim, M.matrix = MethylatedOxBS_sim,
                   G.matrix = UnMethylatedTAB_sim, H.matrix = MethylatedTAB_sim,iterative=TRUE)

results_exact <- MLML(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                      L.matrix = UnMethylatedOxBS_sim, M.matrix = MethylatedOxBS_sim,
                      G.matrix = UnMethylatedTAB_sim, H.matrix = MethylatedTAB_sim)


results_em2 <- MLML2(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                   L.matrix = UnMethylatedOxBS_sim, M.matrix = MethylatedOxBS_sim,
                   G.matrix = UnMethylatedTAB_sim, H.matrix = MethylatedTAB_sim,iterative=TRUE)

results_exact2 <- MLML2(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                      L.matrix = UnMethylatedOxBS_sim, M.matrix = MethylatedOxBS_sim,
                      G.matrix = UnMethylatedTAB_sim, H.matrix = MethylatedTAB_sim)



data(MethylatedBS_sim)
data(MethylatedOxBS_sim)
data(UnMethylatedBS_sim)
data(UnMethylatedOxBS_sim)
data(MethylatedTAB_sim)
data(UnMethylatedTAB_sim)


MethylatedBS_sim[1:10,1] <- NA
UnMethylatedBS_sim[1:10,1] <- NA
UnMethylatedOxBS_sim[1:10,1] <- NA
MethylatedOxBS_sim[1:10,1] <- NA
UnMethylatedTAB_sim[1:10,1] <- NA
MethylatedTAB_sim[1:10,1] <- NA

results_em <- MLML(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                   L.matrix = UnMethylatedOxBS_sim, M.matrix = MethylatedOxBS_sim,
                   G.matrix = UnMethylatedTAB_sim, H.matrix = MethylatedTAB_sim,iterative=TRUE)

results_exact <- MLML(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                      L.matrix = UnMethylatedOxBS_sim, M.matrix = MethylatedOxBS_sim,
                      G.matrix = UnMethylatedTAB_sim, H.matrix = MethylatedTAB_sim)


results_em2 <- MLML2(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                     L.matrix = UnMethylatedOxBS_sim, M.matrix = MethylatedOxBS_sim,
                     G.matrix = UnMethylatedTAB_sim, H.matrix = MethylatedTAB_sim,iterative=TRUE)

results_exact2 <- MLML2(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                        L.matrix = UnMethylatedOxBS_sim, M.matrix = MethylatedOxBS_sim,
                        G.matrix = UnMethylatedTAB_sim, H.matrix = MethylatedTAB_sim)

# obtain MLE via EM-algorithm for BS+TAB:
results_em <- MLML(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                   G.matrix = UnMethylatedTAB_sim, H.matrix = MethylatedTAB_sim,iterative=TRUE)

# obtain constrained exact MLE for BS+TAB:
results_exact <- MLML(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                      G.matrix = UnMethylatedTAB_sim, H.matrix = MethylatedTAB_sim)

# obtain MLE via EM-algorithm for BS+TAB:
results_em2 <- MLML2(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                   G.matrix = UnMethylatedTAB_sim, H.matrix = MethylatedTAB_sim,iterative=TRUE)

# obtain constrained exact MLE for BS+TAB:
results_exact2 <- MLML2(T.matrix = MethylatedBS_sim , U.matrix = UnMethylatedBS_sim,
                      G.matrix = UnMethylatedTAB_sim, H.matrix = MethylatedTAB_sim)

