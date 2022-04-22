library(stringdist)
library(DNABarcodes)
library(data.table)


BC_candidates_dist3 <- create.dnabarcodes(n=7, dist=3, metric = "hamming", heuristic = "ashlock", cores=45 , filter.gc = F)
cand.dt <- data.table::data.table(bc = BC_candidates_dist3)
cand.dt[,paste0("base_",1:7) := tstrsplit(bc, "")]
cand.dt[, id := paste0("BC_",1:.N)]
cand.dt.l <- melt(cand.dt, id.var = 'id')
cand.gc <- cand.dt.l[,.(gc = sum(value %in% c("G","C"))/length(unique(variable))), by = id]
cand.dt <- merge(cand.dt, cand.gc, by = "id")
cand.dt[,gc_penalty := abs(0.5-gc)]
setorder(cand.dt, gc_penalty)

bcs <- cand.dt[1:(13*24)] #13 spike backbones * 24 each
setorder(bcs,id)
bcs[, spike := rep(paste0("spike_",1:13),24)]
setorder(bcs, spike)
fwrite(bcs[,c("bc","spike","gc"), with = F],  file = "/mnt/storage3/chrisz/proj/UMIspike/newSpikes/molspike_barcodes_v2.txt", sep = "\t", quote = F)
