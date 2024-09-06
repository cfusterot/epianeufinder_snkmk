suppressMessages(library("GenomicAlignments"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("plyranges"))
suppressMessages(library("Rsamtools"))
suppressMessages(library("GenomeInfoDb"))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("Biostrings"))
suppressMessages(library("BiocGenerics"))
suppressMessages(library("S4Vectors"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("devtools"))
suppressMessages(library("Signac"))
library("epiAneufinder")
suppressMessages(library("R.utils"))


# -- Analysis set-up -- # 
x = snakemake@params[['sample_id']]
path_to_files <- snakemake@input[['path']]
path_to_out <- snakemake@params[['dir']]
windowSize=snakemake@params[['windowSize']]
blacklist=snakemake@params[['blacklist']]
genome=snakemake@params[['genome']]
reuse.existing=snakemake@params[['reuse']]

# -- Setting directories and paths -- # 
dir.create(file.path(path_to_files,"cnv_custom"))
path_to_fragments <- file.path(path_to_files, "cellranger_count", "outs", "fragments.tsv.gz")
filtered_experiment <- list.files(file.path(path_to_files, "RDS_Files"), full.names = TRUE)
filtered_fragments_path <- file.path(path_to_files,"cnv_custom")

# -- Filter cells in original fragment file -- #
message(paste0("Reading sample ", x, " single-cell RDS"))
atac <- readRDS(filtered_experiment)
cell_names <- unlist(sapply(strsplit(colnames(atac), "_"), `[`, 2, simplify=FALSE))
message("Generating a filtered fragment file")
FilterCells(fragments = path_to_fragments,
            cells = cell_names,
            outfile = file.path(filtered_fragments_path, "filtered_fragments.tsv.gz")) 
message("Decompressing filtered fragment file")
gunzip(file.path(filtered_fragments_path, "filtered_fragments.tsv.gz"), remove=FALSE, overwrite=TRUE)

# -- Setting epianefinder parameters -- # 
# minFrags=2000 was used for the graft samples, while 20000 was used for
# the rest of the samples  
message(paste0("Sample ", x, " contains: ", length(colnames(atac)), " cells"))
message(paste0("Sample ", x, " contains a total of ", sum(atac@meta.data$TSS_fragments), " fragments"))
minFrags = round(sum(atac@meta.data$TSS_fragments)/length(colnames(atac)), digits = -3)
message(paste0("Setting minFrags parameter to ", minFrags))
windowSize= 1e5 #1e5
message(paste0("Setting windowSize parameter to ",windowSize)) 
 
# -- Running epianeufinder -- #
message(paste0("Running epianeufinder for sample ", x))
message(paste0("Analysing ", " fragment file", path_to_fragments))
epiAneufinder(input = file.path(filtered_fragments_path, "filtered_fragments.tsv"), #Enter path to your fragments.tsv file or the folder containing bam files
        outdir = file.path(path_to_out, x, "cnv_custom"), #Path to the directory where results should be written 
        blacklist = blacklist, #Path to .BED containing the blacklisted regions of your genome
        windowSize = windowSize, 
        genome = genome, #Substitute with relevant BSgenome
        #exclude = c('chrX','chrY', 'chrM'),
        reuse.existing = reuse.existing,
        title_karyo = x, 
        ncores = 1,
        minFrags = minFrags)   
