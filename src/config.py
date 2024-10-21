import secrets

from utils import *

RNASUITE_VERSION = "0.1.0"
RNASUITE_BUILD = "20231230"

RNASUITE_PORT = get_unused_port()
RNASUITE_TOKEN = secrets.token_hex(5)

#table type: 1 for input, 2 for output
class RNASuiteDataTables(ClassDict):
	read_count = RNASuiteTable(1, "Read Counts")
	gene_fpkm = RNASuiteTable(1, "FPKM")
	gene_tpm = RNASuiteTable(1, "TPM")
	sample_info = RNASuiteTable(1, "Sample Info")
	gene_annot = RNASuiteTable(1, "Gene Annotation")
	go_annot = RNASuiteTable(1, "GO Annotation")
	kegg_annot = RNASuiteTable(1, "KEGG Annotation")
	normal_count = RNASuiteTable(2, "Normalized Counts")
	degs_list = RNASuiteTable(2, "DEGs")
	go_enrich = RNASuiteTable(2, "GO Enrichment")
	kegg_enrich = RNASuiteTable(2, "KEGG Enrichment")

class RNASuiteErrors(ClassDict):
	import_error = RNASuiteError("Import file {} error", True)
	count_error = RNASuiteError("There are no read count matrix imported")
	sample_error = RNASuiteError("There are no sample information imported")
	degs_error = RNASuiteError("There are no identified DEGs")

class RNASuitePackages(ClassDict):
	_orders = ['Base', 'DEG', 'PCA', 'Enrichment', 'Plot', 'Annotation']
	Base = [
		RNASuitePackage('R', "R is a free software environment for statistical computing and graphics", 'CRAN', '4.3.2'),
		RNASuitePackage('httpgd', "A graphics device for R that is accessible via network protocols", 'CRAN', '2.0.2'),
		RNASuitePackage('reticulate', "Interface to Python modules, classes, and functions", 'CRAN', '1.34.0')
	]
	DEG = [
		RNASuitePackage('DESeq2', "Differential gene expression analysis based on the negative binomial distribution", 'Bioconductor', '1.42.0'),
		RNASuitePackage('edgeR', "Empirical Analysis of Digital Gene Expression Data in R", 'Bioconductor', '4.0.2'),
		RNASuitePackage('limma', "Linear Models for Microarray Data", 'Bioconductor', '3.58.1')
	]
	PCA = [
		RNASuitePackage('glmpca', "Dimension Reduction of Non-Normally Distributed Data", 'CRAN', '0.2.0'),
	]
	Enrichment = [
		RNASuitePackage('clusterProfiler', "A universal enrichment tool for interpreting omics data", 'Bioconductor', '4.10.0')
	]
	Plot = [
		RNASuitePackage('ggplot2', "Create Elegant Data Visualisations Using the Grammar of Graphics", 'CRAN', '3.4.4'),
		RNASuitePackage('ggpubr', "ggplot2 Based Publication Ready Plots", 'CRAN', '0.6.0'),
		RNASuitePackage('ggrepel', "ggrepel provides geoms for ggplot2 to repel overlapping text labels", 'CRAN', '0.9.4'),
		RNASuitePackage('ggvenn', "An easy-to-use way to draw pretty venn diagram by ggplot2", 'CRAN', '0.1.10'),
		RNASuitePackage('ComplexHeatmap', "Arrange multiple heatmaps and various annotation graphics", 'Bioconductor', '2.18.0'),
		RNASuitePackage('UpSetR', "Alternative to Venn and Euler Diagrams for Visualizing Intersecting Sets", 'CRAN', '1.4.0')
	]
	Annotation = [
		RNASuitePackage('org.Hs.eg.db', "Genome wide annotation for Human", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Mm.eg.db', "Genome wide annotation for Mouse", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Rn.eg.db', "Genome wide annotation for Rat", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Dm.eg.db', "Genome wide annotation for Fly", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Dr.eg.db', "Genome wide annotation for Zebrafish", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.At.tair.db', "Genome wide annotation for Arabidopsis", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Sc.sgd.db', "Genome wide annotation for Yeast", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Ce.eg.db', "Genome wide annotation for Worm", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Bt.eg.db', "Genome wide annotation for Bovine", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Gg.eg.db', "Genome wide annotation for Chicken", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Ss.eg.db', "Genome wide annotation for Pig", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Mmu.eg.db', "Genome wide annotation for Rhesus", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Cf.eg.db', "Genome wide annotation for Canine", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.EcK12.eg.db', "Genome wide annotation for E coli strain K12", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Xl.eg.db', "Genome wide annotation for Xenopus", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Pt.eg.db', "Genome wide annotation for Chimp", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Ag.eg.db', "Genome wide annotation for Anopheles", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.EcSakai.eg.db', "Genome wide annotation for E coli strain Sakai", 'Bioconductor', '3.18.0'),
		RNASuitePackage('org.Mxanthus.db', "Genome wide annotation for Myxococcus xanthus DK 1622", 'Bioconductor', '3.18.0'),
	]

RNASUITE_SETTINGS = {
	'R': {
		'binary': ("", str),
		'cran_mirror': ("http://cran.rstudio.com/", str),
		'bioc_mirror': ("", str)
	}
}
