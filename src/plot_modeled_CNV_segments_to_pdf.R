#!/usr/bin/env Rscript
# Austin Southard-Smith
# Coming back after debugging to the point where stuff is working: The way some of this script works is a crime. I am sorry to anybody who tries to update/modify it.
# I found the code from GATK used for plotting the modeled segments and am going to try to modify it to use ggplot so I can save it as an editable PDF with certain elements rasterized instead of the full thing.
# https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/copynumber/plotting/PlotModeledSegments.java
# https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/copynumber/plotting/PlottingUtils.java
# https://github.com/broadinstitute/gatk/blob/master/src/main/resources/org/broadinstitute/hellbender/tools/copynumber/plotting/CNVPlottingLibrary.R
# https://github.com/broadinstitute/gatk/blob/master/src/main/resources/org/broadinstitute/hellbender/tools/copynumber/plotting/PlotModeledSegments.R
library(optparse)
library(data.table)
# libraries I am loading to plot myself.
library(tidyverse)
library(ggrastr)
library(ggpubr)
library(ggrepel)
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
#library(ggnewscale)
# library(ggrepel) I don't plan on labeling but am including this in case.
option_list = list(
    make_option(c("--sample_name"), dest="sample_name", type="character", metavar="character"),
    make_option(c("--denoised_copy_ratios_file"), dest="denoised_copy_ratios_file", type="character", metavar="character"),
    make_option(c("--allelic_counts_file"), dest="allelic_counts_file", type="character", metavar="character"),
    make_option(c("--copy_ratio_call_file"), dest="copy_ratio_call_file", type="character", metavar="character"),
    make_option(c("--modeled_segments_file"), dest="modeled_segments_file", type="character", metavar="character"),
    make_option(c("--sequence_dictionary"), dest="sequence_dictionary", default=NA, type="character", metavar="character"),    # the same sequence dictionary file as what is used by GATK4SCNA
    make_option(c("--contig_names"), dest="contig_names", default=NA, type="character", metavar="character"),      # string with elements separated by ",". only used if sequence dictionary is not provided. 
    make_option(c("--contig_lengths"), dest="contig_lengths", default=NA, type="character", metavar="character"),  # string with elements separated by ",". only used if sequence dictionary is not provided.
    make_option(c("--cytoband_file"), dest="cytoband_file", type="character", metavar="character"),
    make_option(c("--gene_location_file"), dest="gene_location_file", type="character", metavar="character"),
    make_option(c("--maximum_copy_ratio"), dest="maximum_copy_ratio", default = 4, type="double"),
    make_option(c("--point_size_copy_ratio"), dest="point_size_copy_ratio", default=0.4, type="double"),
    make_option(c("--point_size_allele_fraction"), dest="point_size_allele_fraction", default=0.4, type="double"),
    make_option(c("--output_folder"), dest="output_folder", type="character", metavar="character"),
    make_option(c("--retain_contigs"), dest="retain_contigs", default=NA, type="character", metavar="character"),
    make_option(c("--plot_genes"), dest="plot_genes", default='', type="character", metavar="character"))

opt = parse_args(OptionParser(option_list=option_list))

# 
# sample_name = "HT065B1-S1H1A3Y3D1_1"
# denoised_copy_ratios_file = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/HTAN/freeze_v2_2024-07-23/HT065B1-S1H1A3Y3D1_1/HT065B1-S1H1A3Y3D1_1.T.denoisedCR.tsv"
# allelic_counts_file = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/HTAN/freeze_v2_2024-07-23/HT065B1-S1H1A3Y3D1_1/HT065B1-S1H1A3Y3D1_1.T.hets.tsv"
# copy_ratio_call_file = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/HTAN/freeze_v2_2024-07-23/HT065B1-S1H1A3Y3D1_1/HT065B1-S1H1A3Y3D1_1.T.called.seg"
# contig_names_string  = paste(sequence_dict_df_sub$contig, collapse="CONTIG_DELIMITER")
# contig_lengths_string = paste(as.character(sequence_dict_df_sub$contig_length), collapse="CONTIG_DELIMITER")
# maximum_copy_ratio = 4
# point_size_copy_ratio = 0.2
# point_size_allele_fraction = 0.4
# output_file = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/loss_of_heterozygosity/freeze_v2_2024-07-22/HT065B1-S1H1A3N1/HT065B1-S1H1A3Y3D1_1.T.PlotModeledSegments.modeled.pdf"
# cytoband_file = "/diskmnt/Projects/Users/austins2/tools/inferCNV/cytoBand.txt"
# gene_location_file = "/diskmnt/Projects/Users/austins2/tools/inferCNV/genes.Cellranger-2020-A.refdata-gex-GRCh38-2020-A.gene_filterd.remove_excess_columns.trimmed.no_duplicates.txt"


ReadTSV = function(tsv_file) {
    # We need to filter out header lines beginning with '@';
    # however, the standard 'fread("grep ...")' causes issues with the default Docker container, so we use a temporary file.
    # See https://github.com/broadinstitute/gatk/issues/4140.
    temp_file = tempfile()
    system(sprintf('grep -v ^@ "%s" > %s', tsv_file, temp_file))
    return(suppressWarnings(fread(temp_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, data.table=FALSE, showProgress=FALSE, verbose=FALSE)))
}

modeled_segments_file = opt[["modeled_segments_file"]]
modeled_segments_df = ReadTSV(modeled_segments_file)
if (is.na(opt[["retain_contigs"]])) {
    retained_contigs = unique(modeled_segments_df$CONTIG)
} else {
    retained_contigs = unlist(strsplit(opt[["retain_contigs"]], ","))
    if (sum(unique(modeled_segments_df$CONTIG) %in% retained_contigs) > 0) {
        modeled_segments_df = modeled_segments_df[(modeled_segments_df$CONTIG %in% retained_contigs),]
    } else {
        retained_contigs = unique(modeled_segments_df$CONTIG)
    }
}
# print(opt[["contig_names"]])
# print(opt[["contig_lengths"]])
if (is.na(opt[["contig_names"]]) & is.na(opt[["contig_lengths"]])) {
    if (is.na(opt[["sequence_dictionary"]])) {
        stop("you must provide either a sequence dictionary file (e.g. --sequence_dictionary /diskmnt/Datasets/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.dict) or both of the following: a comma separated list of contig/chromosome names and a comma separated list of contig lengths (e.g. --contig_names 'chr1,chr2,chr17' --contig_lengths '248956422,242193529,83257441'")
    }
    sequence_dict_df = read.table(opt[["sequence_dictionary"]], skip = 1, header = F)
    colnames(sequence_dict_df) = c("boilerplate","contig_orig","contig_length_orig","m5_orig","url")
    sequence_dict_df$contig = gsub("^.{0,3}", "", sequence_dict_df$contig_orig)
    sequence_dict_df$contig_length = gsub("^.{0,3}", "", sequence_dict_df$contig_length_orig)
    sequence_dict_df_sub = sequence_dict_df[(sequence_dict_df$contig %in% retained_contigs),]
    
    contig_names_string  = paste(sequence_dict_df_sub$contig, collapse="CONTIG_DELIMITER")
    contig_lengths_string = paste(as.character(sequence_dict_df_sub$contig_length), collapse="CONTIG_DELIMITER")
    contig_names = as.list(strsplit(contig_names_string, "CONTIG_DELIMITER")[[1]]) 
    contig_lengths = as.list(strsplit(contig_lengths_string, "CONTIG_DELIMITER")[[1]])
} else {
    if (is.na(opt[["contig_names"]]) | is.na(opt[["contig_lengths"]])) {
        stop("you must provide either a sequence dictionary file (e.g. -- sequence_dictionary /diskmnt/Datasets/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.dict) or both of the following: a comma separated list of contig/chromosome names and a comma separated list of contig lengths (e.g. --contig_names 'chr1,chr2,chr17' --contig_lengths '248956422,242193529,83257441'")
    }
    contig_names_string = opt[["contig_names"]] #this is a string of of contig names compiled from the sequence dictionary where each name is separated by the string "CONTIG_DELIMITER". It is the same order as the contig_lengths_string. It will be turned into a list of contig names below. I need to make this string and save it to a file.
    contig_lengths_string = opt[["contig_lengths"]] #this is a string of of contig names compiled from the sequence dictionary where each name is separated by the string "CONTIG_DELIMITER". It is the same order as the contig_names_string. It will be turned into a list of contig names below. I need to make this string and save it to a file.
    contig_names = as.list(strsplit(contig_names_string, ",")[[1]]) 
    contig_lengths = as.list(strsplit(contig_lengths_string, ",")[[1]])
    retained_contig_indices = match(retained_contigs, unlist(contig_names))
    contig_names = contig_names[retained_contig_indices]
    contig_lengths = contig_lengths[retained_contig_indices]
}
sample_name = opt[["sample_name"]]
denoised_copy_ratios_file = opt[["denoised_copy_ratios_file"]]
allelic_counts_file = opt[["allelic_counts_file"]]
cytoband_file = opt[["cytoband_file"]]
maximum_copy_ratio = opt[["maximum_copy_ratio"]]
point_size_copy_ratio = opt[["point_size_copy_ratio"]]
point_size_allele_fraction = opt[["point_size_allele_fraction"]]
output_folder = opt[["output_folder"]]
output_prefix = paste0(output_folder,"/",sample_name,"_PlotModeledSegments")
copy_ratio_call_file = opt[["copy_ratio_call_file"]]
gene_location_file = opt[["gene_location_file"]]
# contig_ends = cumsum(contig_lengths)
# contig_starts = c(0, head(contig_ends, -1))
# contig_middles = round((contig_starts + contig_ends) / 2)
if (nchar(opt[["plot_genes"]]) > 0) {
    genes_to_mark = unlist(strsplit(opt[["plot_genes"]], ","))
} else {
    genes_to_mark = c()
}
#print(modeled_segments_df)
# print(contig_names)
# print(contig_lengths)

# modeled_segments_file = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/HTAN/freeze_v2_2024-07-23/HT065B1-S1H1A3Y3D1_1/HT065B1-S1H1A3Y3D1_1.T.modelFinal.seg"
# sample_name = "HT065B1-S1H1A3Y3D1_1"
# denoised_copy_ratios_file = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/HTAN/freeze_v2_2024-07-23/HT065B1-S1H1A3Y3D1_1/HT065B1-S1H1A3Y3D1_1.T.denoisedCR.tsv"
# allelic_counts_file = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/HTAN/freeze_v2_2024-07-23/HT065B1-S1H1A3Y3D1_1/HT065B1-S1H1A3Y3D1_1.T.hets.tsv"
# copy_ratio_call_file = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/HTAN/freeze_v2_2024-07-23/HT065B1-S1H1A3Y3D1_1/HT065B1-S1H1A3Y3D1_1.T.called.seg"
# maximum_copy_ratio = 4
# point_size_copy_ratio = 0.2
# point_size_allele_fraction = 0.4
# output_file = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/loss_of_heterozygosity/freeze_v2_2024-07-22/HT065B1-S1H1A3N1/HT065B1-S1H1A3Y3D1_1.T.PlotModeledSegments.modeled.pdf"
# cytoband_file = "/diskmnt/Projects/Users/austins2/tools/inferCNV/cytoBand.txt" # must have the following columns: c("contig", "start", "end", "name","Giemsa_stain"). Columns can have different names, but must have 5 columns. the names of the bands in the 4th column must contain either a "p" or a "q" somewhere in the string of characters. The chromosomes in the contig column (first column) must be identical format of the contig names in the sequence dictionary file.
# gene_location_file = "/diskmnt/Projects/Users/austins2/tools/inferCNV/genes.Cellranger-2020-A.refdata-gex-GRCh38-2020-A.gene_filterd.remove_excess_columns.trimmed.no_duplicates.txt" #The table must have the following 4 columns c("gene","chromosome","start","stop"). There can be more but the first columns must be those previously listed. The chromosomes in the chromosome column (second column) must be identical format of the contig names in the sequence dictionary file.
# sequence_dictionary = "/diskmnt/Datasets/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.dict" # This must be the same sequence dictionary that was used for initial CNV calling.
# contig_names_string  = paste(sequence_dict_df_sub$contig, collapse="CONTIG_DELIMITER")
# contig_lengths_string = paste(as.character(sequence_dict_df_sub$contig_length), collapse="CONTIG_DELIMITER")

# ReadTSV = function(tsv_file) {
#     # We need to filter out header lines beginning with '@';
#     # however, the standard 'fread("grep ...")' causes issues with the default Docker container, so we use a temporary file.
#     # See https://github.com/broadinstitute/gatk/issues/4140.
#     temp_file = tempfile()
#     system(sprintf('grep -v ^@ "%s" > %s', tsv_file, temp_file))
#     return(suppressWarnings(fread(temp_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, data.table=FALSE, showProgress=FALSE, verbose=FALSE)))
# }
# 
# modeled_segments_df = ReadTSV(modeled_segments_file)
# 
# retained_contigs <- unique(modeled_segments_df$CONTIG)
# sequence_dict_df = read.table(sequence_dictionary, skip = 1, header = F)
# colnames(sequence_dict_df) = c("boilerplate","contig_orig","contig_length_orig","m5_orig","url")
# sequence_dict_df$contig = gsub("^.{0,3}", "", sequence_dict_df$contig_orig)
# sequence_dict_df$contig_length = gsub("^.{0,3}", "", sequence_dict_df$contig_length_orig)
# sequence_dict_df_sub <- sequence_dict_df[(sequence_dict_df$contig %in% retained_contigs),]
# 
# contig_names_string  = paste(sequence_dict_df_sub$contig, collapse="CONTIG_DELIMITER")
# contig_lengths_string = paste(as.character(sequence_dict_df_sub$contig_length), collapse="CONTIG_DELIMITER")
# contig_names = as.list(strsplit(contig_names_string, "CONTIG_DELIMITER")[[1]]) 
# contig_lengths = as.list(strsplit(contig_lengths_string, "CONTIG_DELIMITER")[[1]])
# contig_ends = cumsum(contig_lengths)
# contig_starts = c(0, head(contig_ends, -1))
# contig_middles = round((contig_starts + contig_ends) / 2)


# SetUpPlot(sample_name, "denoised copy ratio", 0, maximum_denoised_copy_ratio, "contig", contig_names, contig_starts, contig_ends, TRUE)
SetUpPlot = function(sample_name, y.lab, y.min, y.max, x.lab, contig_names, contig_starts, contig_ends, contig_middles, do_label_contigs) {
    num_contigs = length(contig_names)
    contig_centers = (contig_starts + contig_ends) / 2
    genome_length = contig_ends[num_contigs]
    use.col = rep(c("grey90", "white"), length.out =num_contigs)
    ggplot_object = ggplot() + 
        geom_rect(aes(xmin = contig_starts, xmax = contig_ends, ymin = floor(y.min), ymax = ceiling(y.max)), fill = NA, color = "black") +
        scale_x_continuous(breaks = contig_middles, labels = contig_names, guide = guide_axis(n.dodge=2)) +
        scale_y_continuous(limits = c(y.min, y.max), breaks = seq(from = floor(y.min), to = ceiling(y.max), by = 1)) + theme_classic() +
        theme(text = element_text(color = "black", size = 8),
              axis.text.x = element_text(color = "black", size = 8),
              axis.text.y = element_text(color = "black", size = 8),
              #legend.position="none",
              legend.text = element_text(color = "black", size = 8),
              axis.title=element_text(size=12),
              axis.ticks = element_line(color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank()) + # element_rect(colour = "black", fill=NA, linewidth=2)) +
        labs(y = paste0(y.lab),
             x = paste0(x.lab))
    pdf(paste0(output_prefix,"_setup.pdf"), width = 12, height = 3.5 * 1, useDingbats=F)
    print(ggplot_object)
    dev.off()
    return(ggplot_object)
}
PlotCopyRatiosWithModeledSegments_genes = function(ggplot_object, denoised_copy_ratios_df, modeled_segments_df, gene_location_file, contig_names, contig_starts, genes_highlight, point_size=0.2) {
    points_start_index = 1
    mean_copy_ratio_vec = c()
    genomic_coordinates_vec = c()
    denoised_copy_ratios_vec = c()
    calls_vec = c()
    segment_coordinates_vec = c()
    segment_vec = c()
    plot_genes = F
    if (length(genes_highlight) > 0) {
        gene_location_df = read.table(gene_location_file,sep='\t',header = F)
        colnames(gene_location_df) = c("gene","chromosome","start","stop")
        gene_location_df$middle = floor((gene_location_df$start + gene_location_df$stop) / 2)
        modeled_segments_df$genes = NA
        modeled_segments_df$gene_middles = NA
        modeled_segments_df$gene_maxes = NA
        for (gene in genes_highlight) {
            if (gene %in% gene_location_df$gene) {
                gene_max = layer_scales(ggplot_object)$y$range$range[2]*0.9
                gene_middle = gene_location_df$middle[gene_location_df$gene == gene]
                gene_contig = gene_location_df$chromosome[gene_location_df$gene == gene]
                gene_offset = contig_starts[match(gene_contig, contig_names)]
                gene_middle_offset = gene_offset + gene_middle
                is_gene_in_any_segment = sum((modeled_segments_df$START < gene_middle) & (modeled_segments_df$END > gene_middle) & (modeled_segments_df$CONTIG == gene_contig)) > 0
                if (is_gene_in_any_segment) {
                    plot_genes = T
                    first_gene_in_segment = is.na(modeled_segments_df$genes[(modeled_segments_df$START < gene_middle) & (modeled_segments_df$END > gene_middle) & (modeled_segments_df$CONTIG == gene_contig)])
                    if (first_gene_in_segment) {
                        modeled_segments_df$genes[(modeled_segments_df$START < gene_middle) & (modeled_segments_df$END > gene_middle) & (modeled_segments_df$CONTIG == gene_contig)] = gene
                        modeled_segments_df$gene_middles[(modeled_segments_df$START < gene_middle) & (modeled_segments_df$END > gene_middle) & (modeled_segments_df$CONTIG == gene_contig)] = as.character(gene_middle_offset)
                        modeled_segments_df$gene_maxes[(modeled_segments_df$START < gene_middle) & (modeled_segments_df$END > gene_middle) & (modeled_segments_df$CONTIG == gene_contig)] = gene_max
                    } else {
                        prior_genes = modeled_segments_df$genes[(modeled_segments_df$START < gene_middle) & (modeled_segments_df$END > gene_middle) & (modeled_segments_df$CONTIG == gene_contig)]
                        prior_middle = modeled_segments_df$gene_middles[(modeled_segments_df$START < gene_middle) & (modeled_segments_df$END > gene_middle) & (modeled_segments_df$CONTIG == gene_contig)]
                        modeled_segments_df$genes[(modeled_segments_df$START < gene_middle) & (modeled_segments_df$END > gene_middle) & (modeled_segments_df$CONTIG == gene_contig)] = paste(prior_genes,gene,sep=",")
                        modeled_segments_df$gene_middles[(modeled_segments_df$START < gene_middle) & (modeled_segments_df$END > gene_middle) & (modeled_segments_df$CONTIG == gene_contig)] = paste(prior_middle,as.character(gene_middle_offset),sep=",")
                    }
                }
            }
        }
    }
    gene_x = c()
    gene_y = c()
    gene_labels = c()
    for (s in 1:nrow(modeled_segments_df)) {
        #skip segments with no points
        num_points = modeled_segments_df[s, "NUM_POINTS_COPY_RATIO"]
        if (num_points == 0) {
            next
        }
        points_end_index = points_start_index + num_points
        
        contig = modeled_segments_df[s, "CONTIG"]
        offset = contig_starts[match(contig, contig_names)]
        
        segment_start = offset + modeled_segments_df[s, "START"]
        segment_end = offset + modeled_segments_df[s, "END"]
        genomic_coordinates = offset + denoised_copy_ratios_df[points_start_index:points_end_index, "MIDDLE"]
        genomic_coordinates_vec = c(genomic_coordinates_vec, genomic_coordinates)
        segment_coordinates_vec = c(segment_coordinates_vec,segment_start,segment_end)
        segment_vec = c(segment_vec,s,s)
        denoised_copy_ratios = denoised_copy_ratios_df[points_start_index:points_end_index, "COPY_RATIO"]
        denoised_copy_ratios_vec = c(denoised_copy_ratios_vec, denoised_copy_ratios)
        denoised_copy_ratio_calls = denoised_copy_ratios_df[points_start_index:points_end_index, "segment_call"]
        calls_vec = c(calls_vec, denoised_copy_ratio_calls)
        mean_copy_ratio = 2^modeled_segments_df[s, "MEAN_LOG2_COPY_RATIO"]
        mean_copy_ratio_vec = c(mean_copy_ratio_vec, mean_copy_ratio, mean_copy_ratio)
        points_start_index = points_start_index + num_points
        if (plot_genes) {
            if (!(is.na(modeled_segments_df[s,"genes"]))) {
                #print(s)
                gene_csv = modeled_segments_df[s,"genes"]
                middles_csv = modeled_segments_df[s,"gene_middles"]
                #print(paste(gene_csv, middles_csv))
                genes_list = unlist(strsplit(gene_csv, ","))
                middle_list = as.numeric(unlist(strsplit(middles_csv, ",")))
                #print(middle_list)
                double_middles = rep(middle_list, each=2)
                #print(double_middles)
                double_genes = rep(genes_list, each=2)
                gene_labels = c(gene_labels, double_genes)
                gene_x = c(gene_x, double_middles)
                #print(gene_x)
                #print(paste(mean_copy_ratio, gene_max))
                gene_y = c(gene_y, rep(c(mean_copy_ratio, gene_max), length(genes_list)))
            }
        }
    }
    new_segment_vec = c()
    new_mean_cr_vec = c()
    for (i in 1:length(segment_coordinates_vec)) {
        if (i%%3 == 0) {
            new_segment_vec = c(new_segment_vec, segment_coordinates_vec[i-1]+1, segment_coordinates_vec[i])
            new_mean_cr_vec = c(new_mean_cr_vec, NA, mean_copy_ratio_vec[i])
        } else {
            new_segment_vec = c(new_segment_vec, segment_coordinates_vec[i])
            new_mean_cr_vec = c(new_mean_cr_vec, mean_copy_ratio_vec[i])
        }
    }
    calls_vec = factor(calls_vec, levels = c('amplified','neutral','deleted'))
    segment_vec = as.character(segment_vec)
    segment_vec = factor(segment_vec)
    # make a factor of fake alpha values with all of the levels you will need for the additional legend. This is so that we can highjack the alpha aesthetic legend to make a legend for the GATK4 MCMC posterior 10%,50%,90% distribution boxes. I cannot for the life of me figure out how else to do this without using cowplot+grobs and I'm trying to stick strictly to ggplot2 right now.
    # this modifies the example of https://stackoverflow.com/questions/57421042/ggplot2-split-one-legend-two-color-scales-and-delete-another 
# make a factor of fake alpha values with all of the levels you will need for the additional legend. This is so that we can highjack the alpha aesthetic legend to make a legend for the GATK4 MCMC posterior 10%,50%,90% distribution boxes. I cannot for the life of me figure out how else to do this without using cowplot+grobs and I'm trying to stick strictly to ggplot2 right now.
    # this modifies the example of https://stackoverflow.com/questions/57421042/ggplot2-split-one-legend-two-color-scales-and-delete-another 
    alphas_vec = factor(rep((c("test")),length.out=length(segment_vec)))
    #print(head(calls_vec))
    #ggplot_object = ggplot_object +
    #    geom_point(aes(x = genomic_coordinates_vec, y = denoised_copy_ratios_vec, colour = calls_vec), size = point_size) +
    #    scale_colour_manual(name = "segment call", guide = "legend", breaks = c("amplified","neutral","deleted"),values =c("amplified" = "#b21f2c","neutral" = "#878787","deleted" = "#2167ac"), drop = FALSE) + 
    #    theme(legend.position="bottom",
    #          legend.text = element_text(color = "black", size = 8)) +
    #    labs(title = paste0(sample_name, " copy ratios and modeled segments")) +
    #    guides(colour = guide_legend(override.aes = list(size = 4))) +
    #    geom_line(aes(x = segment_coordinates_vec, y = mean_copy_ratio_vec, group = segment_vec, colour = alphas_vec), linewidth=1, show.legend = T) +
#scale_colour_manual(name = " ", breaks=c("test"), values=c("modeled segment mean"="black"), drop = F, guide="legend") +
#        geom_line(aes(x = c(min(segment_coordinates_vec),max(segment_coordinates_vec)), y = c(1,1)), color = "#02c0c5", alpha=0.8, linetype=2, linewidth=1, show.legend = F)
    #calls_vec = factor(calls_vec, levels = c('amplified','neutral','deleted'))
    #segment_vec = as.character(segment_vec)
    #segment_vec = factor(segment_vec)
    #
    # below works but has a line through the points of the geom_point colors legend for segment call
    #
    #calls_vec = factor(calls_vec, levels = c('amplified','neutral','deleted'))
    #segment_vec = as.character(segment_vec)
    #segment_vec = factor(segment_vec)
    #alphas_vec = factor(rep(as.double(c(0.999999989999999,1)),length.out=length(segment_vec)), levels = c("0.999999989999999","1"))
    #ggplot_object = ggplot_object +
    #    geom_point(aes(x = genomic_coordinates_vec, y = denoised_copy_ratios_vec, colour = calls_vec), size = point_size) +
    #    scale_colour_manual(values = c(
    #        "amplified" = "#b21f2c",
    #        "neutral" = "#878787",
    #        "deleted" = "#2167ac")) +
    #    theme(legend.position="bottom",
    #          legend.text = element_text(color = "black", size = 8)) +
    #    labs(title = paste0(sample_name, " copy ratios and modeled segments")) +
    #    guides(colour = guide_legend(override.aes = list(size = 4))) +
    #    geom_line(aes(x = segment_coordinates_vec, y = mean_copy_ratio_vec,alpha = alphas_vec, group = segment_vec), linewidth=1, show.legend = T) +
    #    geom_line(aes(x = c(min(segment_coordinates_vec),max(segment_coordinates_vec)), y = c(1,1)), color = "#02c0c5", alpha=0.8, linetype=2, linewidth=0.5, show.legend = F) +
    #    scale_alpha_manual(values=c(1, 1), labels = c(" ","modeled segment mean"), guide=F) +
    
    #    guides(alpha = guide_legend(
    #            #title="modeled segment mean", # label comes from https://gatk.broadinstitute.org/hc/en-us/articles/360035890011--How-to-part-II-Sensitively-detect-copy-ratio-alterations-and-allelic-segments#8.1
    #            override.aes = list(colour = c("white", "black"), 
    #                                size=c(5,5), 
    #                                drop = F,
    #                                shape = c(0,0),
    #                                alpha = c(1,1),
    #                                linetype = c(1,1))))
    #
    # below this might work if you are able to get ggnewscale to install and run correctly. The colors were right and the line was separate, however the line was still overlapping the points. 
    #
    # alphas_vec = factor(rep((c("test")),length.out=length(segment_vec)))
    # ggplot_object = ggplot_object +
    #     geom_point(aes(x = genomic_coordinates_vec, y = denoised_copy_ratios_vec, colour = calls_vec), size = point_size) +
    #     scale_colour_manual(values = c(
    #         "amplified" = "#b21f2c",
    #         "neutral" = "#878787",
    #         "deleted" = "#2167ac")) +
    #     theme(legend.position="bottom",
    #           legend.text = element_text(color = "black", size = 8)) +
    #     labs(title = paste0(sample_name, " copy ratios and modeled segments")) +
    #     guides(colour = guide_legend(override.aes = list(size = 4))) +
    #     new_scale_color() +
    #     geom_line(aes(x = segment_coordinates_vec, y = mean_copy_ratio_vec, group = segment_vec, color = alphas_vec), linewidth=1, show.legend = T) +
    #     scale_color_manual(name = " ", values = c("test" = "black"),labels = c("modeled_segment_mean")) +
    #     geom_line(aes(x = c(min(segment_coordinates_vec),max(segment_coordinates_vec)), y = c(1,1)), color = "#02c0c5", alpha=0.8, linetype=2, linewidth=0.5, show.legend = F)
    # trying scale_linetype_manual
    #
    alphas_vec = factor(rep((c("test")),length.out=length(segment_vec)))
    ggplot_object = ggplot_object +
        geom_point(aes(x = genomic_coordinates_vec, y = denoised_copy_ratios_vec, colour = calls_vec), size = point_size) +
        scale_colour_manual(name ="segment call",
            values = c(
            "amplified" = "#b21f2c",
            "neutral" = "#878787",
            "deleted" = "#2167ac"))# +
    num_colors = length(unique(ggplot_build(ggplot_object)$data[[2]]["colour"])$colour)
    #print(unique(ggplot_build(ggplot_object)$data[[2]]["colour"])$colour)
    #print(num_colors)
    #print(paste0("num_colors ",num_colors))
    ggplot_object = ggplot_object + theme(legend.position="bottom",
              legend.text = element_text(color = "black", size = 8)) +
        labs(title = paste0(sample_name, " copy ratios and modeled segments")) +
        guides(colour = guide_legend(override.aes = list(size = 4))) +
        geom_line(aes(x = segment_coordinates_vec, y = mean_copy_ratio_vec, group = segment_vec, linetype = alphas_vec), color = "black", linewidth=1, show.legend = T) +
        scale_linetype_manual(name = " ", values = c("test" = "solid"), labels = c("modeled segment mean")) +
        geom_line(aes(x = c(min(segment_coordinates_vec),max(segment_coordinates_vec)), y = c(1,1)), color = "#02c0c5", alpha=0.8, linetype=2, linewidth=0.5, show.legend = F) +
    guides(color = guide_legend(override.aes = list(linetype = NA,size = rep(5,num_colors)))) #e(5,5,5))))
    #
    # below this works bud doesn't have a legend for the segment.
    #
    # calls_vec = factor(calls_vec, levels = c('amplified','neutral','deleted'))
    # segment_vec = as.character(segment_vec)
    # segment_vec = factor(segment_vec)
    # ggplot_object = ggplot_object +
    #    geom_point(aes(x = genomic_coordinates_vec, y = denoised_copy_ratios_vec, colour = calls_vec), size = point_size) +
    #    scale_colour_manual(values = c(
    #        "amplified" = "#b21f2c",
    #        "neutral" = "#878787",
    #        "deleted" = "#2167ac")) +
    #    theme(legend.position="bottom",
    #          legend.text = element_text(color = "black", size = 8)) +
    #    labs(title = paste0(sample_name, " copy ratios and modeled segments")) +
    #    guides(colour = guide_legend(override.aes = list(size = 4))) +
    #    geom_line(aes(x = segment_coordinates_vec, y = mean_copy_ratio_vec, group = segment_vec), linewidth=1, show.legend = F) +
    #    geom_line(aes(x = c(min(segment_coordinates_vec),max(segment_coordinates_vec)), y = c(1,1)), color = "#02c0c5", alpha=0.8, linetype=2, linewidth=0.5, show.legend = F)
    if (plot_genes) {
        gene_x_start = as.numeric(gene_x[seq(1,length(gene_x), 2)])
        gene_x_stop = as.numeric(gene_x[seq(2,length(gene_x), 2)])
        gene_y_start = as.numeric(gene_y[seq(1,length(gene_y), 2)])
        gene_y_stop = as.numeric(gene_y[seq(2,length(gene_y), 2)])
        gene_symbol = factor(gene_labels[seq(1, length(gene_labels), 2)])
        gene_df = data.frame(gene_x_start = gene_x_start,
                             gene_x_stop = gene_x_stop,
                             gene_y_start = gene_y_start,
                             gene_y_stop = gene_y_stop,
                             gene_symbol = gene_symbol)
        #print(gene_df)
        ggplot_object = ggplot_object +
            #coord_cartesian(ylim = c(min(gene_y_start)*0.1, max(gene_y_stop)*1.1), clip = 'off') +
            #geom_segment(data = gene_df, aes(x = gene_x_start, xend= gene_x_stop, y = gene_y_start, yend = gene_y_stop, group = gene_symbol), color = "#4baf4d", linetype=1, linewidth=0.5, show.legend = F) +
            #geom_text(aes(x = gene_x_start, y = gene_y_stop + gene_y_stop*0.1, label = gene_symbol), color = "#4baf4d", size = 4)
            #+
            geom_text_repel(size = 3, show.legend = FALSE, min.segment.length = unit(0, 'lines'), direction = "y", box.padding = 1, max.overlaps = Inf, segment.size = 0.5, aes( x=gene_x_start, y=gene_y_start, label=gene_symbol), color = "black") #+
            #theme(plot.margin = margin(1,0,0,0, "in"))
    }
        #guides(group = "none")
    return(ggplot_object)
}
#source("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/loss_of_heterozygosity/HT065B1-S1H1A3N1/PlotCopyRatiosWithModeledSegments_genes.R")

PlotAlternateAlleleFractionsWithModeledSegments = function(ggplot_object, allelic_counts_df, modeled_segments_df, contig_names, contig_starts, point_size=0.2) {
    points_start_index = 1
    box_minor_allele_10_vec = c()
    minor_allele_50_vec = c()
    box_minor_allele_90_vec = c()
    box_segment_start_vec = c()
    box_segment_stop_vec = c()
    box_major_allele_10_vec = c()
    major_allele_50_vec = c()
    box_major_allele_90_vec = c()
    genomic_coordinates_vec = c()
    alternate_allele_fractions_vec = c()
    calls_vec = c()
    segment_coordinates_vec = c()
    segment_vec = c()
    for (s in 1:nrow(modeled_segments_df)) {
        #skip segments with no points
        num_points = modeled_segments_df[s, "NUM_POINTS_ALLELE_FRACTION"]
        if (num_points == 0) {
            next
        }
        points_end_index = points_start_index + num_points
        contig = modeled_segments_df[s, "CONTIG"]
        offset = contig_starts[match(contig, contig_names)]
        segment_start = offset + modeled_segments_df[s, "START"]
        segment_end = offset + modeled_segments_df[s, "END"]
        genomic_coordinates = offset + allelic_counts_df[points_start_index:points_end_index, "POSITION"]
        genomic_coordinates_vec = c(genomic_coordinates_vec, genomic_coordinates)
        segment_coordinates_vec = c(segment_coordinates_vec,segment_start,segment_end)
        segment_vec = c(segment_vec,s,s)
        ref_counts = allelic_counts_df[points_start_index:points_end_index, "REF_COUNT"]
        alt_counts = allelic_counts_df[points_start_index:points_end_index, "ALT_COUNT"]
        alternate_allele_fractions = alt_counts / (alt_counts + ref_counts)
        alternate_allele_fractions_vec = c(alternate_allele_fractions_vec, alternate_allele_fractions)
        allelic_counts_calls = allelic_counts_df[points_start_index:points_end_index, "segment_call"]
        calls_vec = c(calls_vec, allelic_counts_calls)
        minor_allele_fraction_posterior_10 = modeled_segments_df[s, "MINOR_ALLELE_FRACTION_POSTERIOR_10"]
        minor_allele_fraction_posterior_50 = modeled_segments_df[s, "MINOR_ALLELE_FRACTION_POSTERIOR_50"]
        minor_allele_fraction_posterior_90 = modeled_segments_df[s, "MINOR_ALLELE_FRACTION_POSTERIOR_90"]
        box_minor_allele_10_vec = c(box_minor_allele_10_vec, minor_allele_fraction_posterior_10)
        minor_allele_50_vec = c(minor_allele_50_vec, minor_allele_fraction_posterior_50, minor_allele_fraction_posterior_50)
        box_minor_allele_90_vec = c(box_minor_allele_90_vec, minor_allele_fraction_posterior_90)
        box_segment_start_vec = c(box_segment_start_vec, segment_start)
        box_segment_stop_vec = c(box_segment_stop_vec, segment_end)
        major_allele_fraction_posterior_10 = 1 - minor_allele_fraction_posterior_10
        major_allele_fraction_posterior_50 = 1 - minor_allele_fraction_posterior_50
        major_allele_fraction_posterior_90 = 1 - minor_allele_fraction_posterior_90
        box_major_allele_10_vec = c(box_major_allele_10_vec, major_allele_fraction_posterior_10)
        major_allele_50_vec = c(major_allele_50_vec, major_allele_fraction_posterior_50, major_allele_fraction_posterior_50)
        box_major_allele_90_vec = c(box_major_allele_90_vec, major_allele_fraction_posterior_90)
        points_start_index = points_start_index + num_points
    }
    calls_vec = factor(calls_vec, levels = c('amplified','neutral','deleted'))
    segment_vec = as.character(segment_vec)
    segment_vec = factor(segment_vec)
    # make a factor of fake alpha values with all of the levels you will need for the additional legend. This is so that we can highjack the alpha aesthetic legend to make a legend for the GATK4 MCMC posterior 10%,50%,90% distribution boxes. I cannot for the life of me figure out how else to do this without using cowplot+grobs and I'm trying to stick strictly to ggplot2 right now.
    alphas_vec = factor(rep(as.double(c(0.999999989999999,1)),length.out=length(calls_vec)))
    ggplot_object = ggplot_object +
        geom_point(aes(x = genomic_coordinates_vec, y = alternate_allele_fractions_vec, colour = calls_vec, alpha = alphas_vec), size = point_size) +
        theme(legend.position="bottom",
              legend.text = element_text(color = "black", size = 8)) +
        labs(title = paste0(sample_name, " allele frequencies")) +
        guides(colour = guide_legend(override.aes = list(size = 4))) +
        geom_line(aes(x = segment_coordinates_vec, y = minor_allele_50_vec, group = segment_vec), color = "#4baf4d", linewidth=0.5, show.legend = F) + #minor allele line
        geom_rect(aes(xmin = box_segment_start_vec, xmax = box_segment_stop_vec, ymin = box_minor_allele_10_vec, ymax = box_minor_allele_90_vec), fill = NA, color = "#4baf4d") + # minor allele box show.legend=T just resulted in the current legend being removed
        geom_line(aes(x = segment_coordinates_vec, y = major_allele_50_vec, group = segment_vec), color = "#984f9d", linewidth=0.5, show.legend = F) + #major allele line
        geom_rect(aes(xmin = box_segment_start_vec, xmax = box_segment_stop_vec, ymin = box_major_allele_90_vec, ymax = box_major_allele_10_vec), fill = NA, color = "#984f9d") + # major allele box
        geom_line(aes(x = c(min(segment_coordinates_vec),max(segment_coordinates_vec)), y = c(0.5,0.5)), color = "#02c0c5", alpha=0.8, linetype=2, linewidth=0.5, show.legend = T) +
        scale_colour_manual(name = "segment call", breaks = c('amplified','neutral','deleted'), guide = "legend", values =c("amplified" = "#b21f2c","neutral" = "#878787","deleted" = "#2167ac")) #+
    #print(unique(ggplot_build(ggplot_object)$data[[2]]["colour"]))
    num_colors = length(unique(ggplot_build(ggplot_object)$data[[2]]["colour"])$colour)
    #print(num_colors)
    #print(paste0("num_colors ",num_colors))
    ggplot_object = ggplot_object + 
        scale_alpha_manual(values=c(1, 1),labels = c("minor allele", "major allele"), guide=T) +
        guides(alpha = guide_legend(
            title="GATK4 MCMC posterior distribution", # label comes from https://gatk.broadinstitute.org/hc/en-us/articles/360035890011--How-to-part-II-Sensitively-detect-copy-ratio-alterations-and-allelic-segments#8.1
            override.aes = list(colour = c("#4baf4d", "#984f9d"), 
                                #fill = c(NA,NA),
                                size=c(5,5), 
                                drop = F,
                                shape = c(0,0), # to get the lines to show up in the boxes I would need to have a second alpha guide that also 
                                alpha = c(1,1),
                                #linewidth=c(1,1),
                                linetype = c(1,1))))+
        guides(color = guide_legend(override.aes = list(linetype = NA, size = rep(c(5),num_colors-1))))
                                #label = c("minor allele", "major allele"))))
    # Figure out how to add the line+box to the legend (that is going to be a pain in the ass since it isn't a box plot).
    return(ggplot_object)
    # the below works but it only gives a single combined legend for the segment calls and the GATK4 MCMC posterior distributions (major and minor alleles) which is less than ideal
    # calls_vec = factor(calls_vec, levels = c('amplified','neutral','deleted',"Minor allele GATK4 MCMC posterior distribution","Major allele GATK4 MCMC posterior distribution"))
    # segment_vec = as.character(segment_vec)
    # segment_vec = factor(segment_vec)
    # ggplot_object = ggplot_object +
    #     geom_point(aes(x = genomic_coordinates_vec, y = alternate_allele_fractions_vec, colour = calls_vec), size = point_size) +
    #     theme(legend.position="bottom",
    #           legend.text = element_text(color = "black", size = 8)) +
    #     labs(title = paste0(sample_name, " allele frequencies")) +
    #     guides(colour = guide_legend(override.aes = list(size = 4))) +
    #     geom_line(aes(x = segment_coordinates_vec, y = minor_allele_50_vec, group = segment_vec), color = "#4baf4d", linewidth=0.5, show.legend = F) + #minor allele line
    #     geom_rect(aes(xmin = box_segment_start_vec, xmax = box_segment_stop_vec, ymin = box_minor_allele_10_vec, ymax = box_minor_allele_90_vec), fill = NA, color = "#4baf4d") + # minor allele box show.legend=T just resulted in the current legend being removed
    #     #scale_colour_manual(values="#4baf4d") +
    #     #scale_colour_manual(values = c("GATK4 MCMC minor allele posterior summary" = "#4baf4d")) + # label comes from https://gatk.broadinstitute.org/hc/en-us/articles/360035890011--How-to-part-II-Sensitively-detect-copy-ratio-alterations-and-allelic-segments#8.1
    #     geom_line(aes(x = segment_coordinates_vec, y = major_allele_50_vec, group = segment_vec), color = "#984f9d", linewidth=0.5, show.legend = F) + #major allele line
    #     geom_rect(aes(xmin = box_segment_start_vec, xmax = box_segment_stop_vec, ymin = box_major_allele_90_vec, ymax = box_major_allele_10_vec), fill = NA, color = "#984f9d") + # major allele box
    #     #scale_colour_manual(values = c("GATK4 MCMC major allele posterior summary" = "#984f9d")) + # label comes from https://gatk.broadinstitute.org/hc/en-us/articles/360035890011--How-to-part-II-Sensitively-detect-copy-ratio-alterations-and-allelic-segments#8.1
    #     geom_line(aes(x = c(min(segment_coordinates_vec),max(segment_coordinates_vec)), y = c(0.5,0.5)), color = "#02c0c5", alpha=0.8, linetype=2, linewidth=0.5, show.legend = F) +
    #     scale_colour_manual(name = "segment call", breaks = c('amplified','neutral','deleted',"Minor allele GATK4 MCMC posterior distribution","Major allele GATK4 MCMC posterior distribution"), guide = "legend", values =c("amplified" = "#b21f2c","neutral" = "#878787","deleted" = "#2167ac", "Minor allele GATK4 MCMC posterior distribution"="#4baf4d", "Major allele GATK4 MCMC posterior distribution"="#984f9d"), drop = FALSE) +
    #     guides(colour = guide_legend(
    #         title="segment call",
    #         override.aes = list(colour = c('#b21f2c','#878787','#2167ac',"#4baf4d", "#984f9d"), 
    #                             size=c(5,5,5,5,5), 
    #                             drop = F,
    #                             shape = c(16,16,16,0,0),
    #                             linetype = c(0,0,0,1,1),
    #                             label = c('amplified','neutral','deleted',"Minor allele GATK4 MCMC posterior distribution", "Major allele GATK4 MCMC posterior distribution"))))
    # # Figure out how to add the line+box to the legend (that is going to be a pain in the ass since it isn't a box plot).
    # return(ggplot_object)
}
#source("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/loss_of_heterozygosity/HT065B1-S1H1A3N1/PlotAlternateAlleleFractionsWithModeledSegments.R")
#plotting is extracted to a function for debugging purposes
WriteModeledSegmentsPlot = function(sample_name, allelic_counts_file, denoised_copy_ratios_file, modeled_segments_df, copy_ratio_call_file, cytoband_file, contig_names, contig_lengths, genes_to_mark, gene_location_file, output_dir, output_prefix) {
    #modeled_segments_df = ReadTSV(modeled_segments_file)
    #print(modeled_segments_df)
    #num_plots = ifelse(all(file.exists(c(denoised_copy_ratios_file, allelic_counts_file, cytoband_file))), 3, ifelse(all(file.exists(c(denoised_copy_ratios_file, cytoband_file))), 2, ifelse(all(file.exists(c(cytoband_file, allelic_counts_file))), 2, ifelse(all(file.exists(c(denoised_copy_ratios_file, allelic_counts_file))), 2, 1))))
    # this is for testing purposes
    called_segments_df = ReadTSV(copy_ratio_call_file)
    num_contigs = length(contig_names)
    called_segments_df = called_segments_df[(called_segments_df$CONTIG %in% contig_names),]
    # png(output_file, 12, 3.5 * num_plots, units="in", type="cairo", res=300, bg="white")
    # par(mfrow=c(num_plots, 1), cex=0.75, las=1)
    # sets graphical parameters mfrow will be replaced by ggarrange(). cex is used to format text and will be replaced by theme. las is used to format axis labels as always horizontal. this will be handled by theme.
    #par(mfrow=c(num_plots, 1), cex=0.75, las=1)
    # making the karyotype part of the plot
    contig_ends = cumsum(contig_lengths)
    contig_starts = c(0, head(contig_ends, -1))
    contig_middles = round((contig_starts + contig_ends) / 2)
    #print(cytoband_file)
    if (file.exists(cytoband_file) && cytoband_file != "null") {
        cytoband_df = read.table(cytoband_file,sep='\t',header=F)
        colnames(cytoband_df) = c("contig", "start", "end", "name","Giemsa_stain")
        cytoband_df = cytoband_df[(cytoband_df$contig %in% contig_names),]
        cytoband_df$arm = NA
        cytoband_df$arm[(grepl("p",cytoband_df$name))] <- "p"
        cytoband_df$arm[(grepl("q",cytoband_df$name))] <- "q"
        cytoband_df$arm[cytoband_df$Giemsa_stain == "acen"] <- "centromere"
        
        karyotype_df = data.frame(contig = rep(unlist(contig_names),3),
                                  contig.start = rep(contig_starts, 3),
                                  contig.end = rep(contig_ends,3),
                                  contig.length = rep(unlist(contig_lengths),3))
        karyotype_df = karyotype_df[order(karyotype_df$contig.start),]
        row.names(karyotype_df) = c(1:length(karyotype_df$contig.start))
        karyotype_df$arm <- rep(c("p","centromere","q"),length(unlist(contig_names)))
        karyotype_df$cytoband.start <- NA
        karyotype_df$cytoband.end <- NA
        karyotype_df$cytoband.length = NA
        num_contigs = length(contig_names)
        for (i in 1:num_contigs) {
            contig.name = contig_names[[i]][1]
            #print(contig.name)
            for (arm.name in c("p","centromere","q")) {
                #print(arm.name)
                contig_arm_max = max(cytoband_df$end[((cytoband_df$contig == contig.name) & (cytoband_df$arm == arm.name)) ])
                contig_arm_min = min(cytoband_df$start[(cytoband_df$contig == contig.name) & (cytoband_df$arm == arm.name) ])
                karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == arm.name)),"cytoband.start"] = contig_arm_min
                karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == arm.name)),"cytoband.end"] = contig_arm_max
                karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == arm.name)),"cytoband.length"] = (contig_arm_max - contig_arm_min)
            }
            if (arm.name == "centromere") {
                if (!(karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == "p")), "cytoband.end"] == karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == arm.name)), "cytoband.start"])) {
                    middle = (karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == "p")), "cytoband.end"] + karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == arm.name)), "cytoband.start"]) / 2
                    karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == "p")), "cytoband.end"] = middle
                    karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == arm.name)), "cytoband.start"] = middle
                }
            } else if (arm.name == "q") {
                if (!(karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == "centromere")), "cytoband.end"] == karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == arm.name)), "cytoband.start"])) {
                    middle = (karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == "centromere")), "cytoband.end"] + karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == arm.name)), "cytoband.start"]) / 2
                    karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == "centromere")), "cytoband.end"] = middle
                    karyotype_df[((karyotype_df$contig == contig.name) & (karyotype_df$arm == arm.name)), "cytoband.start"] = middle
                }         
            }
        }
        karyotype_df$x.cytoband.start = NA
        karyotype_df$x.cytoband.end = NA
        karyotype_df$x.cytoband.start[karyotype_df$arm == "p"] = karyotype_df$contig.start[karyotype_df$arm == "p"]
        karyotype_df$x.cytoband.end[karyotype_df$arm == "p"] = karyotype_df$contig.start[karyotype_df$arm == "p"] + karyotype_df$cytoband.end[karyotype_df$arm == "p"]
        karyotype_df$x.cytoband.start[karyotype_df$arm == "centromere"] = karyotype_df$contig.start[karyotype_df$arm == "centromere"] + karyotype_df$cytoband.start[karyotype_df$arm == "centromere"]
        karyotype_df$x.cytoband.end[karyotype_df$arm == "centromere"] = karyotype_df$contig.start[karyotype_df$arm == "centromere"] + karyotype_df$cytoband.end[karyotype_df$arm == "centromere"]
        karyotype_df$x.cytoband.start[karyotype_df$arm == "q"] = karyotype_df$contig.start[karyotype_df$arm == "q"] + karyotype_df$cytoband.start[karyotype_df$arm == "q"]
        karyotype_df$x.cytoband.end[karyotype_df$arm == "q"] = karyotype_df$contig.start[karyotype_df$arm == "q"] + karyotype_df$cytoband.end[karyotype_df$arm == "q"]
        orig = karyotype_df
        for (contig in unique(karyotype_df$contig)) {
            if (!(karyotype_df$contig.start[(karyotype_df$contig == contig) & (karyotype_df$arm == "p")] == karyotype_df$x.cytoband.start[(karyotype_df$contig == contig) & (karyotype_df$arm == "p")])) {
                karyotype_df$x.cytoband.start[(karyotype_df$contig == contig) & (karyotype_df$arm == "p")] = karyotype_df$contig.start[(karyotype_df$contig == contig) & (karyotype_df$arm == "p")]
            }
            if (!(karyotype_df$x.cytoband.end[(karyotype_df$contig == contig) & (karyotype_df$arm == "p")] == karyotype_df$x.cytoband.start[(karyotype_df$contig == contig) & (karyotype_df$arm == "centromere")])) {
                karyotype_df$x.cytoband.start[(karyotype_df$contig == contig) & (karyotype_df$arm == "centromere")] = karyotype_df$x.cytoband.end[(karyotype_df$contig == contig) & (karyotype_df$arm == "p")]
            }
            if (!(karyotype_df$x.cytoband.end[(karyotype_df$contig == contig) & (karyotype_df$arm == "centromere")] == karyotype_df$x.cytoband.start[(karyotype_df$contig == contig) & (karyotype_df$arm == "q")])) {
                karyotype_df$x.cytoband.end[(karyotype_df$contig == contig) & (karyotype_df$arm == "centromere")] = karyotype_df$x.cytoband.start[(karyotype_df$contig == contig) & (karyotype_df$arm == "q")]
            }
            if (!(karyotype_df$x.cytoband.end[(karyotype_df$contig == contig) & (karyotype_df$arm == "q")] == karyotype_df$contig.end[(karyotype_df$contig == contig) & (karyotype_df$arm == "q")])) {
                karyotype_df$x.cytoband.end[(karyotype_df$contig == contig) & (karyotype_df$arm == "q")] = karyotype_df$contig.end[(karyotype_df$contig == contig) & (karyotype_df$arm == "q")]
            }
        }
        if (!(unique(orig$x.cytoband.end == karyotype_df$x.cytoband.end))) {
            print("modified the end coordinates of at least one of the cytobands to match the contig coordinates on the plot")
        }
        if (!(unique(orig$x.cytoband.start == karyotype_df$x.cytoband.start))) {
            print("modified the start coordinates of at least one of the cytobands to match the contig coordinates on the plot")
        }
        karyotype_df$arm = factor(karyotype_df$arm, levels = c("p","centromere","q"))
        #print(contig_names)
        #print(contig_starts)
        #print(contig_ends)
        #print(contig_middles)
        print("making p3")
        p3 = SetUpPlot(sample_name, "denoised copy ratio", 0, 1, "contig", contig_names, contig_starts, contig_ends, contig_middles, TRUE)
        p3 = p3 + 
            geom_rect(data = karyotype_df, aes(xmin = x.cytoband.start, xmax = x.cytoband.end, ymin = 0, ymax = 1, fill = arm)) +
            scale_fill_manual(values = c(
                "p" = "#7f7f7f",
                "centromere" = "#2efce7",
                "q" = "#cccccc"), ) +
            geom_rect(aes(xmin = contig_starts, xmax = contig_ends, ymin = 0, ymax = 1), fill = NA, color = "black") +
            #scale_x_continuous(breaks = contig_middles, labels = contig_names, guide = guide_axis(n.dodge=2)) +
            #scale_y_continuous(limits = c(y.min, y.max), breaks = seq(from = floor(y.min), to = ceiling(y.max), by = 1)) +theme_classic() +
            theme(text = element_text(color = "black", size = 8),
                  #axis.text.x = element_text(color = "black", size = 8),
                  #axis.text.y = element_text(color = "black", size = 8),
                  legend.position="bottom",
                  legend.text = element_text(color = "black", size = 8),
                  axis.title=element_text(size=12),
                  axis.ticks = element_line(color = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.line.x = element_blank(),
                  axis.title.x = element_blank(),
                  #axis.title.y = element_blank(),
                  #axis.ticks.y = element_blank(),
                  #axis.line.y = element_blank(),
                  title=element_blank(),
                  axis.text.x = element_blank()) + # element_rect(colour = "black", fill=NA, linewidth=2)) +
            labs(title = paste0(sample_name, " karyotype"),
                 y = paste0("denoised copy ratio"),
                 x = paste0("karyotype"))
        
        pdf(paste0(output_prefix,"_karyotype.pdf"), width = 12, height = 3.5 * 1, useDingbats=F)
        print(p3)
        dev.off()
    }
    if (file.exists(denoised_copy_ratios_file) && denoised_copy_ratios_file != "null") {
        print("preparing p1 input")
        print(modeled_segments_df)
        denoised_copy_ratios_df = ReadTSV(denoised_copy_ratios_file)
        denoised_copy_ratios_df = denoised_copy_ratios_df[(denoised_copy_ratios_df$CONTIG %in% contig_names),]
        #transform to linear copy ratio
        denoised_copy_ratios_df[["COPY_RATIO"]] = 2^denoised_copy_ratios_df[["LOG2_COPY_RATIO"]]
        
        #determine copy-ratio midpoints
        denoised_copy_ratios_df[["MIDDLE"]] = round((denoised_copy_ratios_df[["START"]] + denoised_copy_ratios_df[["END"]]) / 2)
        denoised_copy_ratios_df$segment_call = NA
        modeled_segments_df$MEAN_LOG2_COPY_RATIO = NA
        for (seg in 1:nrow(called_segments_df)) {
            #print(seg)
            seg_contig = called_segments_df[seg,"CONTIG"]
            seg_start = called_segments_df[seg,"START"]
            seg_stop = called_segments_df[seg,"END"]
            seg_call = ifelse(called_segments_df[seg,"CALL"] == "0", "neutral", ifelse(called_segments_df[seg,"CALL"] == "+","amplified","deleted"))
            seg_log2_cr = called_segments_df[seg,"MEAN_LOG2_COPY_RATIO"]
            denoised_copy_ratios_df$segment_call[(denoised_copy_ratios_df$CONTIG == seg_contig) & (denoised_copy_ratios_df$START >= seg_start) & (denoised_copy_ratios_df$END <= seg_stop)] = seg_call
            modeled_segments_df$MEAN_LOG2_COPY_RATIO[(modeled_segments_df$CONTIG == seg_contig) & (modeled_segments_df$START >= seg_start) & (modeled_segments_df$END <= seg_stop)] = seg_log2_cr
        }
        #print(modeled_segments_df)
        # plot up to maximum_copy_ratio (or full range, if maximum_copy_ratio = Infinity)
        maximum_denoised_copy_ratio = if(is.finite(maximum_copy_ratio)) maximum_copy_ratio else 1.05 * max(denoised_copy_ratios_df[["COPY_RATIO"]])
        print("making p1")
        p1 <- SetUpPlot(sample_name, "denoised copy ratio", 0, maximum_denoised_copy_ratio, "contig", contig_names, contig_starts, contig_ends, contig_middles, TRUE)
        p1 <- PlotCopyRatiosWithModeledSegments_genes(p1, denoised_copy_ratios_df, modeled_segments_df, gene_location_file, contig_names, contig_starts, genes_to_mark, point_size_copy_ratio)
        p1 <- rasterize(p1, layers='Point', dpi=600)
        pdf(paste0(output_prefix,"_segments_genes.pdf"), width = 12, height = 3.5 * 1, useDingbats=F)
        print(p1)
        dev.off()
    }
    
    if (file.exists(allelic_counts_file) && allelic_counts_file != "null") {
        print("preparing p2 input")
	allelic_counts_df = ReadTSV(allelic_counts_file)
        allelic_counts_df = allelic_counts_df[(allelic_counts_df$CONTIG %in% contig_names),]
        allelic_counts_df$segment_call = NA
        for (seg in 1:nrow(called_segments_df)) {
            #print(seg)
            seg_contig = called_segments_df[seg,"CONTIG"]
            seg_start = called_segments_df[seg,"START"]
            seg_stop = called_segments_df[seg,"END"]
            seg_call = ifelse(called_segments_df[seg,"CALL"] == "0", "neutral", ifelse(called_segments_df[seg,"CALL"] == "+","amplified","deleted"))
            seg_log2_cr = called_segments_df[seg,"MEAN_LOG2_COPY_RATIO"]
            allelic_counts_df$segment_call[(allelic_counts_df$CONTIG == seg_contig) & (allelic_counts_df$POSITION >= seg_start) & (allelic_counts_df$POSITION <= seg_stop)] = seg_call
            modeled_segments_df$MEAN_LOG2_COPY_RATIO[(modeled_segments_df$CONTIG == seg_contig) & (modeled_segments_df$START >= seg_start) & (modeled_segments_df$END <= seg_stop)] = seg_log2_cr
        }
#source("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/loss_of_heterozygosity/HT065B1-S1H1A3N1/PlotAlternateAlleleFractionsWithModeledSegments.R")
        print("making p2")
	p2 = SetUpPlot(sample_name, "alternate-allele fraction", 0, 1.0, "contig", contig_names, contig_starts, contig_ends, contig_middles, TRUE)
        p2 = PlotAlternateAlleleFractionsWithModeledSegments(p2, allelic_counts_df, modeled_segments_df, contig_names, contig_starts, point_size_allele_fraction)
        p2 <- rasterize(p2, layers='Point', dpi=600)
        pdf(paste0(output_prefix,"_alleles.pdf"), width = 12, height = 3.5 * 1, useDingbats=F)
        print(p2)
        dev.off()
    }
    # p3 = SetUpPlot(sample_name, "denoised copy ratio", 0, maximum_denoised_copy_ratio, "contig", contig_names, contig_starts, contig_ends, contig_middles, TRUE)
    # p3 = p3 + 
    #     geom_rect(data = karyotype_df, aes(xmin = x.cytoband.start, xmax = x.cytoband.end, ymin = 0, ymax = maximum_denoised_copy_ratio, fill = arm)) +
    #     scale_fill_manual(values = c(
    #         "p" = "#7f7f7f",
    #         "centromere" = "#2efce7",
    #         "q" = "#cccccc"), ) +
    #     geom_rect(aes(xmin = contig_starts, xmax = contig_ends, ymin = 0, ymax = maximum_denoised_copy_ratio), fill = NA, color = "black") +
    #     #scale_x_continuous(breaks = contig_middles, labels = contig_names, guide = guide_axis(n.dodge=2)) +
    #     #scale_y_continuous(limits = c(y.min, y.max), breaks = seq(from = floor(y.min), to = ceiling(y.max), by = 1)) +theme_classic() +
    #     theme(text = element_text(color = "black", size = 8),
    #           #axis.text.x = element_text(color = "black", size = 8),
    #           #axis.text.y = element_text(color = "black", size = 8),
    #           legend.position="bottom",
    #           legend.text = element_text(color = "black", size = 8),
    #           axis.title=element_text(size=12),
    #           axis.ticks = element_line(color = "black"),
    #           panel.grid.major = element_blank(),
    #           panel.grid.minor = element_blank(),
    #           panel.background = element_blank(),
    #           panel.border = element_blank(),
    #           axis.ticks.x = element_blank(),
    #           axis.line.x = element_blank(),
    #           axis.title.x = element_blank(),
    #           #axis.title.y = element_blank(),
    #           #axis.ticks.y = element_blank(),
    #           #axis.line.y = element_blank(),
    #           title=element_blank(),
    #           axis.text.x = element_blank()) + # element_rect(colour = "black", fill=NA, linewidth=2)) +
    #     labs(title = paste0(sample_name, " copy ratio and modeled segments"),
    #          y = paste0("denoised copy ratio"),
    #          x = paste0(x.lab))
    #     
    # pdf(paste0(output_prefix,"_karyotype.pdf"), width = 12, height = 3.5 * 1, useDingbats=F)
    # p3
    # dev.off()
    total_plots = sum(c(exists("p1"),exists("p2"),exists("p3")))
    if (total_plots > 0) {
        if (is.na(opt[["retain_contigs"]])) {
            plot_width = 12
        } else {
            plot_width = (2 + (1*num_contigs))
        }
        if (all(file.exists(c(denoised_copy_ratios_file, allelic_counts_file, cytoband_file)))) {
            print("combining denoised_copy_ratios_file, allelic_counts_file, cytoband_file")
            print(c(rep(5,total_plots-1),1))
            pdf(paste0(output_prefix,"_complete.pdf"), width = plot_width, height = 3 * total_plots, useDingbats=F)
            print(ggarrange(p1,p2,p3, ncol = 1, nrow = total_plots, heights = c(rep(4,total_plots-1),1)))
            dev.off()
        } else if (all(file.exists(c(denoised_copy_ratios_file, cytoband_file)))) {
            print("combining denoised_copy_ratios_file, cytoband_file")
            print(c(rep(5,total_plots-1),1))
            pdf(paste0(output_prefix,"_complete.pdf"), width = plot_width, height = 3 * total_plots, useDingbats=F)
            print(ggarrange(p1,p3, ncol = 1, nrow = total_plots, heights = c(rep(4,total_plots-1),1)))
            dev.off()
        } else if (all(file.exists(c(allelic_counts_file, cytoband_file)))) {
            print("combining allelic_counts_file, cytoband_file")
            print(c(rep(5,total_plots-1),1))
            pdf(paste0(output_prefix,"_complete.pdf"), width = plot_width, height = 3 * total_plots, useDingbats=F)
            print(ggarrange(p2,p3, ncol = 1, nrow = total_plots, heights = c(rep(4,total_plots-1),1)))
            dev.off()
        } else if (all(file.exists(c(denoised_copy_ratios_file, allelic_counts_file)))) {
            print("combining denoised_copy_ratios_file, allelic_counts_file")
            print(c(1,1))
            pdf(paste0(output_prefix,"_complete.pdf"), width = plot_width, height = 3 * total_plots, useDingbats=F)
            print(ggarrange(p1,p2, ncol = 1, nrow = total_plots, heights = c(1,1)))
            dev.off()
        }
    }
    #check for created file and quit with error code if not found
    if (!file.exists(paste0(output_prefix,"_complete.pdf"))) {
        quit(save="no", status=1, runLast=FALSE)
    }
}

#WriteModeledSegmentsPlot = function(sample_name, allelic_counts_file, denoised_copy_ratios_file, modeled_segments_file, copy_ratio_call_file, cytoband_file, contig_names, genes_to_mark = c(), contig_lengths, output_dir, output_prefix)
WriteModeledSegmentsPlot(sample_name, allelic_counts_file, denoised_copy_ratios_file, modeled_segments_df, copy_ratio_call_file, cytoband_file, contig_names, contig_lengths, genes_to_mark, gene_location_file, output_dir, output_prefix)

