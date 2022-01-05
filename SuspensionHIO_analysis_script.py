'''
BASIC SINGLE CELL ANALYSIS SCRIPT
by Josh Wu

Modified by Charlie Childs for data analysis of Capeling et al., 2022. 

Relies heavily on the Scanpy Python module developed by the Theis Lab
Read more about Scanpy at https://scanpy.readthedocs.io/en/latest/index.html

'''

from sca_run import *
#from tools.pipelines import *

figdir = './05202020_Figures/'
an_run = sca_run()
#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
an_run.storage_mount_point = 'Z:/'

## IDs of samples as represented in the metadata table
an_run.sample_list = ['234-1']

# List of interesting genes
an_run.add_gene_list(markers= ['MSLN', 'WT1', 'UPK3B', 'KRT19', 'PDPN', 'CDH1', 'EPCAM', 'SHH', 'IHH', 'DHH', 'PTCH1', 'PTCH2', 'SMO', 'GLI1', 'GLI2', 'GLI3', 'HHIP'],
					 label='basic_list')

an_run.add_gene_list(markers= ['MSLN', 'WT1', 'UPK3B', 'KRT19', 'PDPN', 'CDH1', 'EPCAM'],
					 label='dotplot_list')

an_run.add_gene_list(markers= ['GATA4', 'TCF21','TBX18'],
					 label='mesotheliumdevelopment_list')

an_run.add_gene_list(markers= ['CDX2', 'PDX1', 'MUC2', 'CHGA', 'LYZ', 'DPP4', 'MKI67', 'OLFM4', 'LGR5'],
 					 label='epithelial_list')

an_run.add_gene_list(markers= ['WNT1','WNT2','WNT2B','WNT3','WNT3A','WNT4','WNT5A','WNT5B','WNT6','WNT7A','WNT8A','WNT9A','WNT9B','WNT10A','WNT10B','WNT11','WNT16','FZD1','FZD2','FZD3','FZD4','FZD5','FZD6','FZD7','FZD8','FZD9','FZD10','LRP5','LRP6','CTNNB1','APC','GSK3A','GSK3B'],
				 label='WNT_list1')

an_run.add_gene_list(markers= ['CSNK1A1','CSNK1G1','BTRC','AXIN2','MYC'],
				 label='WNT_list2')

an_run.add_gene_list(markers= ['SHH', 'IHH', 'DHH', 'PTCH1', 'PTCH2', 'GAS1', 'CDON', 'BOC', 'SMO', 'GLI1', 'GLI2', 'GLI3', 'KIF7', 'SUFU', 'HHIP'],
					 label='HH_list')

## Parameters used to filter the data - Mainly used to get rid of bad cells
an_run.set_filter_params(min_cells = 0, # Filter out cells 
						 min_genes = 500, # Filter out cells with fewer genes to remove dead cells
						 max_genes = 8000, # Filter out cells with more genes to remove most doublets
						 max_counts = 50000, # Filter out cells with more UMIs to catch a few remaining doublets
						 max_mito = 0.1) # Filter out cells with high mitochondrial gene content

## Parameters used for initial clustering analysis
an_run.set_analysis_params(n_neighbors = 15, # Size of the local neighborhood used for manifold approximation
						   n_pcs = 12, # Number of principle components to use in construction of neighborhood graph
						   spread = 1, # In combination with min_dist determines how clumped embedded points are
						   min_dist = 0.4, # Minimum distance between points on the umap graph
						   resolution = 0.35) # High resolution attempts to increases # of clusters identified

## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
an_run.pipe_basic(figdir)

## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code

# New analysis parameters for the subset of parameters
analysis_params_ext = dict(n_neighbors = 15,
						n_pcs = 12,
						spread = 1,
						min_dist = 0.4,
						resolution = 0.3)

an_run.pipe_ext(analysis_params_ext, figdir=figdir, extracted=['2'], load_save='adata_save.p')


