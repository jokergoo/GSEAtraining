library(org.Hs.eg.db)
org.Hs.egSYMBOL2EG[["TP53"]]

lt = as.list(org.Hs.egGO2ALLEGS)
l = sapply(lt, function(x) {
    "7157" %in% x
})

gs1 = names(lt)[l]

##
lt = as.list(org.Hs.egGO2EG)
l = sapply(lt, function(x) {
    "7157" %in% x
})

gs2 = names(lt)[l]

library(GO.db)

onto = Ontology(gs2)
all_gs = NULL

ind = which(onto == "CC")
all_gs = c(all_gs, gs2[ind], unlist(as.list(GOCCANCESTOR[ gs2[ind] ])))

ind = which(onto == "MF")
all_gs = c(all_gs, gs2[ind], unlist(as.list(GOMFANCESTOR[ gs2[ind] ])))

ind = which(onto == "BP")
all_gs = c(all_gs, gs2[ind], unlist(as.list(GOBPANCESTOR[ gs2[ind] ])))

gs3 = unique(all_gs)
