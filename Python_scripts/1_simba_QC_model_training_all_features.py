#!/usr/bin/env python3
import os
import simba as si
si.settings.set_figure_params(dpi=80,
                              style='white',
                              fig_size=[5,5],
                              rc={'image.cmap': 'viridis'})

from IPython.display import set_matplotlib_formats
set_matplotlib_formats('retina')
import scanpy as sc
import sys
import argparse
import logging
import subprocess

def QC_RNA(args):
	logging.info("RNA QC")
	adata_CG = sc.read(args.input_rna)

	si.pp.filter_genes(adata_CG, min_n_cells=3)
	si.pp.cal_qc_rna(adata_CG)
	si.pp.normalize(adata_CG, method='lib_size')
	si.pp.log_transform(adata_CG)
	
	logging.info("Discretize RNA values in 5 bins")
	si.tl.discretize(adata_CG, n_bins=5)
	si.pl.discretize(adata_CG, kde=False)
	
	workdir = args.output
	
	# 🚫 Clean up reserved '_index' column if present in var
	if '_index' in adata_CG.var.columns:
		adata_CG.var.rename(columns={'_index': 'var_index'}, inplace=True)
	if '_index' in adata_CG.obs.columns:
		adata_CG.obs.rename(columns={'_index': 'obs_index'}, inplace=True)
	# 📌 Set raw after cleanup
	adata_CG.raw = adata_CG  # make sure .raw.var gets the clean version
	# 🔁 Just to be extra safe: check again in raw.var
	if '_index' in adata_CG.raw.var.columns:
		adata_CG.raw.var.rename(columns={'_index': 'raw_var_index'}, inplace=True)
	logging.info("Saving processed RNA")
	adata_CG.write(os.path.join(workdir, 'adata_CG.h5ad'))
	return adata_CG




def QC_ATAC(args):
	logging.info("ATAC QC")
	adata_CP = sc.read(args.input_atac)

	si.pp.filter_genes(adata_CP, min_n_cells=3)
	si.pp.cal_qc_rna(adata_CP)  # assuming ATAC is stored like RNA here
	si.pp.normalize(adata_CP, method='lib_size')
	si.pp.log_transform(adata_CP)

	logging.info("Discretize ATAC values in 5 bins")
	si.tl.discretize(adata_CP, n_bins=5)
	si.pl.discretize(adata_CP, kde=False)

	workdir = args.output

	# 🔒 Clean reserved column names before setting raw
	if '_index' in adata_CP.var.columns:
		adata_CP.var.rename(columns={'_index': 'var_index'}, inplace=True)

	if '_index' in adata_CP.obs.columns:
		adata_CP.obs.rename(columns={'_index': 'obs_index'}, inplace=True)

	adata_CP.raw = adata_CP  # ensure raw gets the cleaned var

	if '_index' in adata_CP.raw.var.columns:
		adata_CP.raw.var.rename(columns={'_index': 'raw_var_index'}, inplace=True)

	logging.info("Saving processed ATAC")
	adata_CP.write(os.path.join(workdir, 'adata_CP.h5ad'))

	return adata_CP

def read_kmers_and_motifs(args):
	logging.info("Read kmers and motifs from Rscript")
	adata_PK = si.read_hdf(args.input_kmers,'mat')
	adata_PM = si.read_hdf(args.input_motifs,'mat')
	logging.info("Convert byte string to string")
	adata_PK.obs.index = [x.decode('utf-8') for x in adata_PK.obs.index]
	adata_PK.var.index = [x.decode('utf-8') for x in adata_PK.var.index]
	adata_PM.obs.index = [x.decode('utf-8') for x in adata_PM.obs.index]
	adata_PM.var.index = [x.decode('utf-8') for x in adata_PM.var.index]
	si.pp.binarize(adata_PK)
	si.pp.binarize(adata_PM)
	workdir = args.output
	adata_PM.write(os.path.join(workdir,'adata_PM.h5ad'))
	adata_PK.write(os.path.join(workdir,'adata_PK.h5ad'))
	return adata_PK, adata_PM

def graph_generation(args, adata_CG, adata_CP, adata_PK, adata_PM):	
	logging.info("Generate Graph")	
	si.tl.gen_graph(list_CP=[adata_CP],
					list_CG=[adata_CG],
					list_PK=[adata_PK],
					list_PM=[adata_PM],
					copy=False,
					use_highly_variable=False,
					use_top_pcs=False,
					dirname='graph0')
	return

def train_model(args):
	logging.info("Training model")
	dict_config = si.settings.pbg_params.copy()
	dict_config['workers'] = args.threads
	si.tl.pbg_train(auto_wd=True, save_wd=True, output='model')
	si.pl.pbg_metrics(fig_ncol=1)
	return

def extract_embeddings(args):
	emb_path = args.path_emb
	logging.info("Set path to the directory for pbg")
	si.load_pbg_config(path=emb_path)
	workdir = args.output
	logging.info("Load pre-existing graph and model")	
	si.load_graph_stats(path=args.preexisting_graph)# load in graph ('graph0') info	
	si.load_pbg_config(path=args.preexisting_model)# load in model info for ('graph0')
	entity_alias_path = args.path_entity_alias
	entity_path = args.path_entity	
	logging.info("Load dictionary")
	#dict_adata = si.read_embedding(path_entity=entity_path, path_entity_alias = entity_alias_path)
	dict_adata = si.read_embedding()	
	keys_list = list(dict_adata.keys())
	print(keys_list)
	logging.info("extracting embeddings")	
	adata_C = dict_adata['C']  # embeddings for cells
	adata_G = dict_adata['G']  # embeddings for genes
	adata_P = dict_adata['P']  # embeddings for peaks
	adata_K = dict_adata['K']  # embeddings for kmers
	adata_M = dict_adata['M']  # embeddings for motifs
	#adata_M.obs.index = 'M_'+adata_M.obs.index #### to distinguish TF motif names from gene names in this case
	adata_M.obs.index = ['M_' + str(idx) for idx in adata_M.obs.index]
	adata_G.write(os.path.join(workdir,'adata_G.h5ad'))
	adata_P.write(os.path.join(workdir,'adata_P.h5ad'))
	adata_K.write(os.path.join(workdir,'adata_K.h5ad'))
	adata_M.write(os.path.join(workdir,'adata_M.h5ad'))
	return adata_C, adata_G, adata_P, adata_K, adata_M


def annotate_and_compare_cell_types(args, adata_CG, adata_C):
	workdir = args.output
	logging.info("Add annotation to cell embedding") ## Add annotation of celltypes (optional)	
	adata_C.obs['celltype'] = adata_CG[adata_C.obs_names,:].obs['refined_annotation_majority_vote'].copy()
	adata_C.obs['clone_line'] = adata_CG[adata_C.obs_names,:].obs['clone_line'].copy()
	adata_C.obs['Genotype'] = adata_CG[adata_C.obs_names,:].obs['Genotype'].copy()
	logging.info("UMAP of cell embedding") ## Add annotation of celltypes (optional)	
	si.tl.umap(adata_C,n_neighbors=15,n_components=2)
	adata_C.write(os.path.join(workdir,'adata_C.h5ad'))
	return adata_C

def gene_compare_entities(args, adata_C, adata_G):
	logging.info("compare gene and  cell entities") ## Add annotation of celltypes (optional)
	adata_cmp_CG = si.tl.compare_entities(adata_ref=adata_C, adata_query=adata_G)
	workdir = args.output
	adata_cmp_CG.write(os.path.join(workdir,'adata_cmp_CG.h5ad'))
	return adata_cmp_CG

def motif_compare_entities(args, adata_C, adata_M):
	logging.info("compare motif and  cell entities") ## Add annotation of celltypes (optional)
	adata_cmp_CM = si.tl.compare_entities(adata_ref=adata_C, adata_query=adata_M)
	workdir = args.output
	adata_cmp_CM.write(os.path.join(workdir,'adata_cmp_CM.h5ad'))
	return adata_cmp_CM

def peak_compare_entities(args, adata_C, adata_P):
	logging.info("compare peaks and  cell entities") ## Add annotation of celltypes (optional)
	adata_cmp_CP = si.tl.compare_entities(adata_ref=adata_C, adata_query=adata_P)
	workdir = args.output
	adata_cmp_CP.write(os.path.join(workdir,'adata_cmp_CP.h5ad'))
	return adata_cmp_CP

def global_embedding(args, adata_C, adata_G, adata_M, adata_K, adata_P):
	workdir = args.output
	logging.info("Put together all embeddings") ## Add annotation of celltypes (optional)
	adata_all = si.tl.embed(adata_ref=adata_C, list_adata_query=[adata_G, adata_M, adata_K, adata_P])
	adata_all.write(os.path.join(workdir,'adata_all.h5ad'))
	return adata_all

def main(args):
	logging.info('Start.')

	# Initialize critical variables

	adata_CG = adata_CP = adata_PK = adata_PM = None

	if not args.skip_RNA_QC:
		logging.info(f'Loading RNA data from {args.input_rna}')
		adata_CG = QC_RNA(args)
	if not args.skip_ATAC_QC:
		logging.info(f'Loading ATAC data from {args.input_atac}')
		adata_CP = QC_ATAC(args)
	if not args.skip_kmers_motifs:
		logging.info(f'Loading kmers data from {args.input_kmers}')
		logging.info(f'Loading motifs data from {args.input_motifs}')		
		logging.info('Read kmers and motifs file')	                
		adata_PK, adata_PM = read_kmers_and_motifs(args)
	if not args.skip_graph_generation and not args.skip_train_model:
		# Check if all required data are available for graph generation
		if adata_CG is None or adata_CP is None or adata_PK is None or adata_PM is None:
			logging.error("Missing required input data for graph generation.")
			sys.exit(1)
		logging.info('generating graph')                
		graph_generation(args, adata_CG, adata_CP, adata_PK, adata_PM)
		logging.info('training model')
		train_model(args)
	if not args.skip_compare_entities:
		logging.info(f'Loading pre existing graph from {args.preexisting_graph}')
		logging.info(f'Loading pre existing model from {args.preexisting_model}')
		# Check if pre-existing graph and model paths are provided
		if args.preexisting_graph is None or args.preexisting_model is None:
			logging.error("Pre-existing graph and model paths must be provided if skipping graph generation or model training.")
			sys.exit(1)		
		adata_C, adata_G, adata_P, adata_K, adata_M = extract_embeddings(args)
		# Check if RNA data (adata_CG) is available for annotation
		if adata_CG is None:
			logging.error("Annotation requires RNA data (adata_CG). Please provide RNA input file.")
			sys.exit(1)	
		logging.info('annotate cells with labels')
		adata_C = annotate_and_compare_cell_types(args, adata_CG, adata_C)
		adata_cmp_CG = gene_compare_entities(args, adata_C, adata_G)
		adata_cmp_CM = motif_compare_entities(args, adata_C, adata_M)
		adata_cmp_CP = peak_compare_entities(args, adata_C, adata_P)
	if not args.skip_global_embedding:
		workdir = args.output
		adata_C = si.read_h5ad(os.path.join(workdir,'adata_C.h5ad'))
		adata_G = si.read_h5ad(os.path.join(workdir,'adata_G.h5ad'))
		adata_M = si.read_h5ad(os.path.join(workdir,'adata_M.h5ad'))		
		adata_P = si.read_h5ad(os.path.join(workdir,'adata_P.h5ad'))
		adata_K = si.read_h5ad(os.path.join(workdir,'adata_K.h5ad'))
		adata_all = global_embedding(args, adata_C, adata_G, adata_M, adata_K, adata_P)
		
	logging.info('Finish.')
	return

def process_args():
	parser = argparse.ArgumentParser(description='Simba initial steps.')
	io_group = parser.add_argument_group('I/O arguments')
	io_group.add_argument('-rna', '--input-rna', required=True, type=str, metavar='FILE', help='Path to input RNA (.h5ad) file')
	io_group.add_argument('-atac', '--input-atac', required=True, type=str, metavar='FILE', help='Path to input ATAC (.h5ad) file')
	io_group.add_argument('-kmers', '--input-kmers', required=True, type=str, metavar='FILE', help='Path to kmer prediction file')
	io_group.add_argument('-motifs', '--input-motifs', required=True, type=str, metavar='FILE', help='Path to motif prediction file')	
	io_group.add_argument('-graph', '--preexisting-graph', required=False, type=str, help='graph preexisting dir')
	io_group.add_argument('-model', '--preexisting-model', required=False, type=str, help='model preexisting dir')
	io_group.add_argument('-entity_alias', '--path-entity_alias', required=False, type=str, help='path to entity_alias')
	io_group.add_argument('-entity', '--path-entity', required=False, type=str, help='path to entity')
	io_group.add_argument('-emb', '--path-emb', required=False, type=str, help='path to pbg model')
	io_group.add_argument('-o', '--output', required=True, type=str, help='Output directory to store processed files')
	



	model_group = parser.add_argument_group('Model arguments')
	model_group.add_argument('-t', '--threads', required=True, type=int, help='Number of threads to use for alignment')
		

	skip_group = parser.add_argument_group('Skip steps')
	skip_group.add_argument('--skip-RNA_QC', required=False, action='store_true', default=False, help='Skip RNA_QC step')
	skip_group.add_argument('--skip-ATAC_QC', required=False, action='store_true', default=False, help='Skip ATAC_QC step')
	skip_group.add_argument('--skip-kmers_motifs', required=False, action='store_true', default=False, help='Skip kmers and motifs step')
	skip_group.add_argument('--skip-graph-generation', required=False, action='store_true', default=False, help='Skip graph generation step')
	skip_group.add_argument('--skip-train-model', required=False, action='store_true', default=False, help='Skip training model step')
	skip_group.add_argument('--skip-compare_entities', required=False, action='store_true', default=False, help='Skip compare entitites step')
	skip_group.add_argument('--skip-global_embedding', required=False, action='store_true', default=False, help='Skip global embedding step')		
	return parser.parse_args()

if __name__ == '__main__':
	logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.INFO)
	args = process_args()
	main(args)

