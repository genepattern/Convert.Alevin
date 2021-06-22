import os, sys
import gc

import argparse
import shutil
import anndata
import scipy
import tarfile
import pandas as pd
import anndata as ad
from vpolo.alevin import parser as vparser

__author__ = "Anthony S. Castanza"
__email__ = "acastanza@ucsd.edu"
__version__="1.0.0"

def main():
	usage="%prog [options]" + "\n"
	ap = argparse.ArgumentParser()
	ap.add_argument("-i","--input-files",action="store",dest="input_file",help="List of .tar.gz Alevin output directories, with unique names.")
	ap.add_argument("-f","--features",action="store",nargs='?',dest="features",help="Two column .features.tsv from PreprocessVelocityTranscriptome")
	ap.add_argument("-m","--merge",action="store",dest="merge",help="Combine multiple Alevin runs into a single analysis object (True/False)")
	ap.add_argument("-o","--out",action="store",dest="output_ext",help="Output results as AnnData (h5ad) or Loom (loom)")
	options = ap.parse_args()

	infile=os.path.normpath(options.input_file)
	infile = open(infile, "r").read().splitlines()

	features = options.features
	outfile = ""
	merge = options.merge
	out_ext = options.output_ext

	if outfile == "":
		outfile = list(map(os.path.basename, infile))
		outfile = list(map(lambda st: str.replace(st, ".tgz", ""), outfile))
		outfile = list(map(lambda st: str.replace(st, ".gz", ""), outfile))
		outfile = list(map(lambda st: str.replace(st, ".tar", ""), outfile))
	#	outfile = outfile + ".loom"
	# else:
	# 	outfile = outfile.replace(".loom","")
	# 	outfile = outfile + ".loom"

	alevin_df = outfile.copy()
	adata_full = outfile.copy()

	for i in range(len(outfile)):
		tar = tarfile.open(os.path.normpath(infile[i]))
		tar.extractall(outfile[i])
		tar.close()

	if features != None:
		features = pd.read_table(features)

	for i in range(len(outfile)):
		if isinstance(features, pd.DataFrame):
			print("Reading in sample data:")
			adata_working=vparser.read_quants_bin(outfile[i])
			print("Getting spliced quantifications")
			alevin_spliced=adata_working[adata_working.columns.intersection(list(features.spliced))]
			print("Getting unspliced quantifications")
			alevin_unspliced=adata_working[adata_working.columns.intersection(list(features.intron))]
			del adata_working
			gc.collect()
			print("Processing spliced matrix")
			alevin_spliced_reindex=alevin_spliced.reindex(list(features.spliced), axis=1)
			alevin_spliced_reindex = alevin_spliced_reindex.fillna(0)
			alevin_spliced_mtx = scipy.sparse.csr_matrix(alevin_spliced_reindex.values)
			mtx_col_names = list(alevin_spliced_reindex.columns)
			mtx_row_names = list(alevin_spliced_reindex.index)
			del alevin_spliced_reindex
			gc.collect()
			print("Processing unspliced matrix")
			alevin_unspliced_reindex=alevin_unspliced.reindex(list(features.intron), axis=1)
			alevin_unspliced_reindex = alevin_unspliced_reindex.fillna(0)
			alevin_unspliced_reindex.columns = alevin_unspliced_reindex.columns.str.replace('-I', '')
			alevin_unspliced_mtx = scipy.sparse.csr_matrix(alevin_unspliced_reindex.values)
			del alevin_unspliced, alevin_unspliced_reindex
			gc.collect()
			ambiguous=pd.DataFrame(index=list(alevin_spliced.index), columns=list(features.spliced))
			ambiguous = ambiguous.fillna(0)
			ambiguous_mtx = scipy.sparse.csr_matrix(ambiguous.values)
			print("Creating anndata structure for sample: "+ str(outfile[i]))
			adata_full[i] = anndata.AnnData(X = alevin_spliced_mtx,
			                        layers = dict(spliced = alevin_spliced_mtx, unspliced = alevin_unspliced_mtx, ambiguous = ambiguous_mtx))
			del alevin_spliced_mtx, alevin_unspliced_mtx
			gc.collect()
			adata_full[i].var_names=mtx_col_names
			adata_full[i].var_names.set_names('Gene')
			adata_full[i].obs_names=mtx_row_names
			adata_full[i].obs_names.set_names('CellID')
		else: # IF no features database exists, parse identities from the alevin result
			print("Reading in sample data:")
			adata_working=vparser.read_quants_bin(outfile[i])
			spliced_members = list(set(adata_working.columns.str.replace('-I', '')))
			unspliced_members = [sub + "-I" for sub in spliced_members]
			print("Getting spliced quantifications")
			alevin_spliced=adata_working[adata_working.columns.intersection(spliced_members)]
			print("Getting unspliced quantifications:")
			alevin_unspliced=adata_working[adata_working.columns.intersection(unspliced_members)]
			del adata_working
			gc.collect()
			print("Processing spliced matrix")
			alevin_spliced_reindex=alevin_spliced.reindex(spliced_members, axis=1)
			alevin_spliced_reindex = alevin_spliced_reindex.fillna(0)
			alevin_spliced_mtx = scipy.sparse.csr_matrix(alevin_spliced_reindex.values)
			mtx_row_names = list(alevin_spliced_reindex.index)
			del alevin_spliced_reindex
			print("Processing unspliced matrix")
			gc.collect()
			alevin_unspliced_reindex=alevin_unspliced.reindex(unspliced_members, axis=1)
			alevin_unspliced_reindex = alevin_unspliced_reindex.fillna(0)
			alevin_unspliced_reindex.columns = alevin_unspliced_reindex.columns.str.replace('-I', '')
			alevin_unspliced_mtx = scipy.sparse.csr_matrix(alevin_unspliced_reindex.values)
			mtx_col_names = list(alevin_unspliced_reindex.columns)
			del alevin_unspliced, alevin_unspliced_reindex
			gc.collect()
			ambiguous=pd.DataFrame(index=list(alevin_spliced.index), columns=spliced_members)
			ambiguous = ambiguous.fillna(0)
			ambiguous_mtx = scipy.sparse.csr_matrix(ambiguous.values)
			print("Creating anndata structure for sample: "+ str(outfile[i]))
			adata_full[i] = anndata.AnnData(X = alevin_spliced_mtx,
			                        layers = dict(spliced = alevin_spliced_mtx, unspliced = alevin_unspliced_mtx, ambiguous = ambiguous_mtx))
		#	adata.Accession=list(alevin_spliced_reindex.columns)
		#	adata.Gene=list(alevin_unspliced_reindex.columns)
		#	adata.CellID=list(alevin_spliced_reindex.index)
			adata_full[i].var_names=mtx_col_names
			adata_full[i].var_names.set_names('Gene')
			adata_full[i].obs_names=mtx_row_names
			adata_full[i].obs_names.set_names('CellID')

	if(merge == "False"):
		if(out_ext == "loom"):
			if (len(adata_full) == 1):
				ad.AnnData.write_loom(adata_full[0],outfile[0]+"."+out_ext) #Write out Loom file
			elif (len(adata_full) > 1):
				for i in range(len(adata_full)):
					ad.AnnData.write_loom(adata_full[i],outfile[i]+"."+out_ext) #Write out Loom file
		if(out_ext == "h5ad"):
			if (len(adata_full) == 1):
				ad.AnnData.write(adata_full[0], compression="gzip", filename=outfile[0]+"."+out_ext) #Write out h5ad file
			elif (len(adata_full) > 1):
				for i in range(len(adata_full)):
					ad.AnnData.write(adata_full[i], compression="gzip", filename=outfile[i]+"."+out_ext) #Write out h5ad file
	if(merge == "True"):
		concat_samples = adata_full[0].concatenate(adata_full[1:len(adata_full)], join='outer', batch_categories=outfile, index_unique='-')
		if(out_ext == "loom"):
			ad.AnnData.write_loom(concat_samples, "Combined_Experiment"+"."+out_ext) #Write out Loom file
		if(out_ext == "h5ad"):
			ad.AnnData.write(concat_samples, compression="gzip", filename="Combined_Experiment"+"."+out_ext) #Write out h5ad file

	for f in outfile:
		shutil.rmtree(f)

if __name__ == '__main__':
	main()
