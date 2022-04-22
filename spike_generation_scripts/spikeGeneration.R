library(data.table)
setDTthreads(12)
library(Biostrings)
library(ggplot2)



# optional code: check distribution of GC & lenght of typical huma --------

#load cDNA sequences of the Human genome
cdna_seq <- Biostrings::readDNAStringSet("~/resources/genomes/Human/Homo_sapiens.GRCh38.95.cDNA.fa")
#load APPRIS data on principal isoforms 
appris <- fread("~/resources/genomes/Human/appris_data.principal.txt", header = F)
#load Ensembl transcript biotype information
biotypes <- fread("~/resources/genomes/Human/Homo_sapiens.GRCh38.95.chr.tx_biotypes.txt", header =F )
#load UMI count tables for UHRR bulk sequencing
uhrr <- read.table("~/Downloads/UHRR_UMIcounts_1mio.txt")
uhrr_exprs <- data.table(GeneID = row.names(uhrr), mu = rowMeans(uhrr))

principal_tx <- cdna_seq[substr(names(cdna_seq),1,15) %in% appris$V3]

principal_tx_dat <- data.frame(Biostrings::alphabetFrequency(principal_tx)[,1:4])
setDT(principal_tx_dat)
principal_tx_dat[,txID := names(principal_tx)][
                 ,txID := substr(txID,1,15)][
                 ,txLength := A+C+G+T][
                 ,txGC := ((G+C)/txLength)*100]
principal_tx_dat <- merge(principal_tx_dat, appris[,c("V2","V3"), with = F], by.x = "txID", by.y = "V3", all.x = TRUE)
principal_tx_dat <- merge(principal_tx_dat, biotypes[,c("V1","V3"), with = F], by.x = "txID", by.y = "V1", all.x = TRUE)
principal_tx_dat <- merge(principal_tx_dat, uhrr_exprs, by.x = "V2", by.y = "GeneID", all.x = TRUE)

table(principal_tx_dat$V3)


summary(principal_tx_dat[mu > 1]$txLength)

ggplot(principal_tx_dat,aes(txLength)) + geom_density() +  scale_x_log10(labels = scales::label_comma()) + annotation_logticks(sides = 'b') + cowplot::theme_cowplot() + geom_vline(xintercept = 2505, lty = 'dashed')
ggplot(principal_tx_dat,aes(txGC)) + geom_density() +  scale_x_continuous(labels = scales::label_comma())  + cowplot::theme_cowplot() + geom_vline(xintercept = c(40,50,60), lty = 'dashed')

# functions ---------------------------------------------------------------

check_max_constant_bases <- function(dna_vector){
  for(i in 1:length(dna_vector)){
    if(i == 1) {
      count <- 0
    }else{
      if(dna_vector[i] != last_base) { count <- c(count,0) } else {count <- c(count,count[i-1]+1)}
    }
    last_base <- dna_vector[i]
  }
  return(max(count))
}
sequence_generator <- function(n_bases, max_n_stretch = 4, gc_content = 50, gc_window = 1, iterations = 1000000, threads = 8){
  seqs <- parallel::mclapply(1:iterations, FUN = function(x) {
    tmp <- sample(x = c("G","C","T","A"), replace = T, size = n_bases, prob = c(gc_content/200,gc_content/200,(100-gc_content)/200,(100-gc_content)/200))

    #remove strings that do not have good GC content
    gc_obs <- (sum(table(tmp)[c("G","C")])/length(tmp))*100
    if(!is.null(gc_window)){
      if(abs(gc_obs-gc_content) > gc_window){
        return(NA)
      }
    }
    
    #check for constant bases
    if(!is.null(max_n_stretch)){
      if(check_max_constant_bases(tmp) > max_n_stretch){
        return(NA)
      }
    }
    
    #make output string
    return_val <- paste(tmp, collapse = '')
    return(return_val)
  }, mc.cores = threads)
  
  seqs <- unlist(seqs)
  seqs <- seqs[!is.na(seqs)]
  
  seqs <- DNAStringSet(seqs)
  
  #names(seqs) <- paste0("sequence_",1:length(seqs))
  return(seqs)
}



# generate randomized candidate sequences ---------------------------------


seqs_gc40_3000bp <- sequence_generator(n_bases = 3000, max_n_stretch = 4, gc_content = 40, gc_window = 1, threads = 20)
names(seqs_gc40_3000bp) <- paste0("gc40_sequence_",1:length(seqs_gc40_3000bp))
seqs_gc50_3000bp <- sequence_generator(n_bases = 3000, max_n_stretch = 4, gc_content = 50, gc_window = 1, threads = 20)
names(seqs_gc50_3000bp) <- paste0("gc50_sequence_",1:length(seqs_gc50_3000bp))
seqs_gc60_3000bp <- sequence_generator(n_bases = 3000, max_n_stretch = 4, gc_content = 60, gc_window = 1, threads = 20)
names(seqs_gc60_3000bp) <- paste0("gc60_sequence_",1:length(seqs_gc60_3000bp))

all_seqs <- c(seqs_gc40_3000bp,seqs_gc50_3000bp,seqs_gc60_3000bp)


# simulate reads & check for mappable fragments against human geno --------

### generate reads for mapping
#BiocManager::install("polyester")
library(polyester)
readspertx <- round(100 * width(all_seqs) / 100) #100x coverage
fold_changes = matrix(rep(1,2*length(all_seqs)), nrow=length(all_seqs))
writeXStringSet(all_seqs, '~/projects/spike/polyester_reads/generated_seqs.fa')

simulate_experiment('~/projects/spike/polyester_reads/generated_seqs.fa', reads_per_transcript=readspertx, outdir='~/projects/spike/polyester_reads/', fold_changes = fold_changes, num_reps = c(1,1), paired = FALSE, gzip = TRUE, readlen = 50) 
system('~/programs/STAR-2.7.3a/bin/Linux_x86_64_static/STAR --readNameSeparator - --outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5 --outFilterMismatchNmax 20 --outSAMmultNmax 20 --outFilterMultimapNmax 1000 --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --runThreadN 24 --sjdbOverhang 49 --sjdbGTFfile ~/resources/genomes/Human/Homo_sapiens.GRCh38.95.chr.gtf --genomeDir ~/resources/genomes/Human/STAR7idx_primary_noGTF --readFilesCommand zcat --readFilesIn ~/projects/spike/polyester_reads/sample_01.fasta.gz,~/projects/spike/polyester_reads/sample_02.fasta.gz --outFileNamePrefix ~/projects/spike/polyester_reads/starmapping. ')

mapping_reads <- fread(cmd = "~/programs/samtools-1.12/samtools view -F4 ~/projects/spike/polyester_reads/starmapping.Aligned.sortedByCoord.out.bam | cut -f1 ", header = F)
mapping_reads[, seq_id := tstrsplit(V1, "[/]", keep = 2)]
mappable_seqs <- unique(mapping_reads$seq_id)

filtered_seqs <- all_seqs[which(!names(all_seqs) %in% mappable_seqs)]




# GC sliding window checking -------------------------------------------------------

window <- 100
# compute the GC content in a sliding window (as a fraction) for a sequence no. 364
gc_windows <- lapply(1:length(filtered_seqs),
            function(x) {
              out <- data.table(gc = rowSums(letterFrequencyInSlidingView(filtered_seqs[[x]], window, c("G", "C")))/window)
              out[,idx := seq(.N)]
              return(out)
              })
names(gc_windows) <- names(filtered_seqs)
gc_windows <- rbindlist(gc_windows, idcol = "seqID")
gc_windows[, gc_content := tstrsplit(seqID,"_", keep = 1)]

ggplot(gc_windows[seqID %in% sample(unique(gc_windows$seqID),300, replace = F)], aes(x = idx, y = gc, group = seqID)) +
  geom_smooth(method = 'loess', alpha = 0.1, se = F, span = 0.5) +
  cowplot::theme_cowplot() +
  facet_wrap(~gc_content)

ggplot(gc_windows[seqID %in% sample(unique(gc_windows$seqID),300, replace = F)], aes(x = seqID, y = gc)) +
  geom_boxplot() +
  cowplot::theme_cowplot() +
  facet_wrap(~gc_content, scales = "free_x")


gc_window_summary <- gc_windows[, .(max_gc = max(gc), min_gc = min(gc), mean_gc = mean(gc)), by = c("seqID","gc_content")]
gc_window_summary[, minmax_diff := max_gc-min_gc]
setorder(gc_window_summary, minmax_diff)


# restriction analysis ----------------------------------------------------

library(DECIPHER)
data(RESTRICTION_ENZYMES)

d <- DigestDNA(RESTRICTION_ENZYMES[c("HindIII")], filtered_seqs, type = "positions")
d <- DigestDNA("TAACTATAACG/GTCCTAAGGTAGCGAA", filtered_seqs, type = "positions")
