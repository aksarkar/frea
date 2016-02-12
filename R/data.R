#' Metadata for 127 epigenomes from Roadmap Epigenomics
#'
#' Unique identifiers, cell type descriptions, tissue level organization, and
#' coloring scheme.
#'
#' https://www.broadinstitute.org/~meuleman/reg2map/sample_info/sample_info.RData
'roadmap_sample_info'

#' Cluster density of enhancers (DNaseI regions selected with -log10(p) >= 2)
#'
#' Regions delineated using observed DNaseI data across 39 epigenomes,
#' annotated with the 5-mark 15-state model based on observed data across 111
#' Roadmap reference epigenomes.
#'
#' 2,328,936 putative enhancer regions (12.6385% of genome) clustered with
#' k=226, for a resolution of ~10k regions per cluster
'honeybadger_cluster_density'

#' Cluster density of enhancers (DNaseI peak -log10(p) >= 2, observed marks)
#'
#' Regions are delineated using observed DNaseI data across 53 epigenomes
#' annotated with the 5-mark 15-state model based on observed data across 127
#' epigenomes (Roadmap + ENCODE)
#'
#' 2,464,858 putative enhancer regions (14.29% of genome) are clustered with
#' k=246, for a resolution of ~10k regions per cluster
'honeybadger2_p2_cluster_density'

#' Cluster density of enhancers (DNaseI peak -log10(10) >= 10, observed marks)
#'
#' Regions are delineated using observed DNaseI data across 53 epigenomes
#' annotated with the 5-mark 15-state model based on observed data across 127
#' epigenomes (Roadmap + ENCODE)
#'
#' 618,429 putative enhancer regions (7.46% of genome) are clustered with k=62,
#' for a resolution of ~10k regions per cluster
'honeybadger2_p10_cluster_density'
