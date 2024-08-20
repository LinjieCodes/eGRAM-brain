import numpy as np
import pandas as pd
import re
import os
from scipy.stats import hypergeom
import sys, getopt


def read_gtf(gtfFile):
	'''
	Takes a Gene Annotation File (i.e., the gtf file from Genecode) and returns a gene dict and a gene set.
	
	The gene dict returned contains EnsemblID and symbol of each protein_coding gene.
	
	The gene set returned contains all protein_coding genes.
	'''
	gene_symbols = {}
	protein_genes = set()
	with open(gtfFile) as f:
		for line in f:
			if '\tgene\t' in line and 'gene_biotype "protein_coding"' in line:
				geneID = re.search('gene_id "(.+?)";', line).group(1)
				symbol = re.search('gene_name "(.+?)";', line).group(1)
				gene_symbols[geneID] = symbol
				protein_genes.add(symbol)
	return gene_symbols, protein_genes
	
	
def read_lncID(lncFile):
	'''
	Read a lncrna list file containing EnsemblID and symbol of each lncrna that is to be analyzed.

	This function is required only when the gene expression matrix is based on EnsemblID instead of gene symbol. 
	'''
	lnc_symbols = {}
	with open(lncFile) as f:
		f.readline()
		for line in f:
			cols = line.split('\t')
			symbol = cols[0]
			lncID = cols[1].split('.')[0]
			lnc_symbols[lncID] = symbol
	return lnc_symbols
	
		
def read_exp(expFile, targets, affinity):
	'''
	Read the gene expression matrix which is based on gene symbol.
	'''
	exp_dict = {}
	with open(expFile) as f:
		samples = f.readline().strip().split('\t')[1:]
		sampleIndexs = range(len(samples))
		for sample in samples:
			exp_dict[sample] = {}
		for line in f:
			cols = line.strip().split('\t')
			symbol = cols[0]
			if symbol in targets or symbol in affinity:
				values = [float(e) for e in cols[1:]]
				for sampleIndex in sampleIndexs:
					sample = samples[sampleIndex]
					exp_dict[sample][symbol] = values[sampleIndex]
	#gene row, sample column
	exp_df = pd.DataFrame(exp_dict)
	return exp_df
	
	
def cal_corr(exp_df):
	'''
	Take a dataframe of gene expression matrix and return a gene correlation dataframe.
	
	The correlation is calculated using Spearman method.
	'''
	#pandas corr's correlation calculating is performed between columns
	#so the gene expression dataframe should be "sample-row and gene-column"
	corr_df = exp_df.T.corr(method='spearman')
	return corr_df
	
	
def read_binding(bindFile, protein_genes):
	'''
	Read the lncrna-target binding matrix.

	Return a dict containg the binding affinity values of each lncrna-target pair, 
	and a gene set containg all protein_coding target genes.
	'''
	affinity = {}
	targets = set()
	with open(bindFile) as f:
		lncrnas = f.readline().strip().split('\t')[1:]
		lncIndexs = range(len(lncrnas))
		for line in f:
			cols = line.strip('\r\n').split('\t')
			symbol = cols[0]
			if symbol in protein_genes:
				targets.add(symbol)
				values = cols[1:]
				for lncIndex in lncIndexs:
					if values[lncIndex].isdigit():
						affi = float(values[lncIndex])
						lncrna = lncrnas[lncIndex]
						if lncrna not in affinity:
							affinity[lncrna] = {}
						affinity[lncrna][symbol] = affi
	return affinity, targets


def find_center(geneSet, corr_df, corrThres):
	'''
	Find the center gene of a given gene set upon the gene correlation relationships.
	
	Take a gene set, a correlation dataframe, and a correlation threshold, and return a center gene and a gene module.
	'''
	corrGeneNums = []
	module = set()
	for gene1 in geneSet:
		num = 0
		for gene2 in geneSet:
			if gene1 != gene2:
				corr_coef = corr_df[gene1][gene2]
				if corr_coef > corrThres:
					num += 1
		corrGeneNums.append((num, gene1))
	corrGeneNums.sort(reverse=True)
	centerNum, centerGene = corrGeneNums[0]
	module_pre = []
	for gene in geneSet:
		corr_coef = corr_df[centerGene][gene]
		if corr_coef > corrThres:
			module_pre.append((corr_coef, gene))
	module_pre.sort(reverse=True)
	for corr_coef, gene in module_pre:
		module.add(gene)
		if len(module)>300:	#small modules are preferred
			break			#when the number of correlated genes > 300, the most correlated one are adopted
	return centerGene, module
	

def exhaust_subset(geneSet):
	'''
	Exhaustively search all subset of a given gene set.
	
	Take a gene set, and return a list of all subsets of the given gene set.
	'''
	subsets = []
	N = len(geneSet)
	geneList = list(geneSet)
	for i in range(2**N):
		zj = []
		for j in range(N):
			if (i >> j) % 2 == 1:
				zj.append(geneList[j])
		if zj:
			subsets.append(tuple(sorted(zj)))
			
	subsets_nums = []
	for subset in subsets:
		subsets_nums.append((len(subset), subset))
	subsets_nums.sort(reverse=True)
	
	subsets = []
	for num, subset in subsets_nums:
		subsets.append(subset)
	return subsets
	

def determine_function(lncrna, module, corr_df):
	'''
	Determine the activator/inhibitior relationships by examining the lncrna is positively or 
	negativelycorrelated with the expression profiles of genes in the corresponding modules.
	
	Return 'A', 'I' or 'N'.
	
	'A' indicates activator, 'I' indicates inhibitor and 'N' indicates 'Not determined'.
	'''
	if lncrna not in corr_df.columns:
		return 'N'
	a_num = i_num = 0
	for target in module:
		if target in corr_df.columns and corr_df[lncrna][target] > 0:
			a_num += 1
		elif target in corr_df.columns and corr_df[lncrna][target] < 0:
			i_num += 1
	moduleSize = len(module)
	if a_num > 0.8*moduleSize:
		return 'A'
	elif i_num > 0.8*moduleSize:
		return 'I'
	else:
		return 'N'
		
		
def read_kegg(geneFile, pathFile, linkFile, categoryFile):
	'''
	Read kegg genes, pathways and pathway-gene links.
	'''
	gene_symbols = {}
	kegg_term = {}
	kegg_link = {}
	gene_id_set = set()
	
	pathway_categ = {}
	with open(categoryFile) as f:
		for line in f:
			category, pathwayName = line.strip().split('\t')
			pathway_categ[pathwayName] = category.split(':')[0]
		
	with open(geneFile) as f:
		for line in f:
			cols = line.strip('\r\n').split('\t')
			gene_id = cols[0]
			descri = cols[-1]
			if (',' in descri) or (';' in descri):
				symbol_descri = descri.split(';')[0]
				symbols = symbol_descri.split(", ")
				if gene_id not in gene_symbols:
					gene_symbols[gene_id] = set()
				for symbol in symbols:
					gene_symbols[gene_id].add(symbol)
				gene_id_set.add(gene_id)
					
	with open(pathFile) as f:
		for line in f:
			pathID, pathwayName = line.strip('\r\n').split('\t')
			pathID = pathID.replace("path:",'')
			pathwayName = pathwayName[:pathwayName.rfind(' -')]
			if pathwayName in pathway_categ:
				kegg_term[pathID] = pathwayName
			
	with open(linkFile) as f:
		for line in f:
			pathID, gene_id = line.strip('\r\n').split('\t')
			pathID = pathID.replace("path:",'')
			if pathID in kegg_term:
				if gene_id in gene_id_set:
					if pathID not in kegg_link:
						kegg_link[pathID] = set()
					kegg_link[pathID].add(gene_id)
				
	totalGeneNum = len(gene_id_set)
	
	kegg_gene = {}
	for pathID in kegg_link:
		for gene_id in kegg_link[pathID]:
			for gene_symbol in gene_symbols[gene_id]:
				kegg_gene[gene_symbol] = gene_id
	
	return kegg_gene, kegg_term, kegg_link, pathway_categ, totalGeneNum
	
	
def calcu_fdr(pValues):
	'''
	Benjamini-Hochberg method for p-value correction
	'''
	pValues = np.asfarray(pValues)
	by_descend = pValues.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(pValues)) / np.arange(len(pValues), 0, -1)
	fdr_values = np.minimum(1, np.minimum.accumulate(steps * pValues[by_descend]))
	fdrs = []
	for i in range(len(fdr_values[by_orig])):
		fdrs.append(fdr_values[by_orig][i])
	return fdrs
	
	
def download_kegg():
	'''
	Download KEGG data.
	'''
	if not os.path.exists('PathwayAnnotation/'):
		os.mkdir('PathwayAnnotation/')
		
	for species in ["hsa","mmu"]:
		gene_url = "http://rest.kegg.jp/list/%s" % species
		gene_file = "PathwayAnnotation/keggGene-%s" % species
		if not os.path.exists(gene_file):
			urlretrieve(gene_url,gene_file)
		
		term_url = "http://rest.kegg.jp/list/pathway/%s" % species
		term_file = "PathwayAnnotation/keggPathway-%s" % species
		if not os.path.exists(term_file):
			urlretrieve(term_url,term_file)
		
		link_url = "http://rest.kegg.jp/link/%s/pathway" % species
		link_file = "PathwayAnnotation/keggLink-%s" % species
		if not os.path.exists(link_file):
			urlretrieve(link_url,link_file)
	
	
def pathway_analysis(module, 
					 geneDict, 
					 pathDict, 
					 linkDict,
					 pathway_categ,
					 totalGeneNum):
	'''
	Perform pathway enrichment analysis.
	'''			
	results = []
	pvalues = []
	
	module_geneID = {}
	for gene in module:
		if gene in geneDict:
			module_geneID[geneDict[gene]] = gene
			
	for pathID in linkDict:
		interSet = set(module_geneID.keys()) & linkDict[pathID]
		hitNum = len(interSet)
		if hitNum:
			thisPathGeneNum = len(linkDict[pathID])
			moduleSize = len(module_geneID)
			
			#hypergeometric test
			pVal = hypergeom.sf(hitNum-1, totalGeneNum, thisPathGeneNum, moduleSize)
			pathNam = pathDict[pathID]
			hitGenes = ', '.join([module_geneID[geneID] for geneID in interSet])
			results.append([pathID, pathNam, hitGenes, pVal])
			pvalues.append(pVal)
		
	fdrs = calcu_fdr(pvalues)
	for i in range(len(fdrs)):
		results[i].append(fdrs[i])
		
	pathway_sorted = []
	for pathID, pathNam, hitGenes, pVal, fdr in results:
		if fdr < 0.01:
			pathway_sorted.append((fdr, pVal, pathNam, pathID, hitGenes))
	pathway_sorted.sort()
	
	if pathway_sorted:
		topFDR, topPVal, topPathway, topPathID, hitGenes = pathway_sorted[0]
		topPathCateg = pathway_categ[topPathway]
		return topFDR, topPVal, topPathway, topPathID, topPathCateg, hitGenes
	else:
		return 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'

	
def identify_module(exp_df, 
					affinity, 
					moduleFile,				
					kegg_gene, 
					kegg_term, 
					kegg_link,
					pathway_categ,					
					totalGeneNum,
					categoryFile,
					affiThres1, 
					affiThres2,
					moduleSize,
					corrThres):
	'''
	Exhaustively searche the space of all possible sets of lncrnas and target gene sets, 
	to identify all modules that are regulated by the given lncrnas.
	
	The results are written in a file named by the moduleFile param.
	'''
	#only expressed lncRNAs are adopted
	exp_lnc_affinity = {}
	for lncrna in affinity:
		if lncrna in exp_df.index:
			exp_lnc_affinity[lncrna] = affinity[lncrna]
	
	lncStrongTargets = {}
	for lncrna in exp_lnc_affinity:
		lncStrongTargets[lncrna] = set()
		for target in exp_lnc_affinity[lncrna]:
			if exp_lnc_affinity[lncrna][target] > affiThres1:
				lncStrongTargets[lncrna].add(target)
			
	lncWeakTargets = {}
	for lncrna in exp_lnc_affinity:
		lncWeakTargets[lncrna] = set()
		for target in exp_lnc_affinity[lncrna]:
			if affiThres2 < exp_lnc_affinity[lncrna][target] <= affiThres1:
				lncWeakTargets[lncrna].add(target)

	#the lncrnas and targets that have been explored will not be examined
	lnc_subSets_explored = set()
	lnc_target_explored = set()
	
	allExpTargets = set(exp_df.index) - set(exp_lnc_affinity.keys())
	allExpTargets_sorted = sorted(list(allExpTargets))
	
	corr_df = cal_corr(exp_df)
	
	modules = []
	for gene in allExpTargets_sorted:
		#lncrnas strongly binding to the gene
		lncSet = set()
		for lncrna in exp_lnc_affinity:
			if gene in exp_lnc_affinity[lncrna] and exp_lnc_affinity[lncrna][gene] > affiThres1:
				lncSet.add(lncrna)
				
		lnc_subSets = exhaust_subset(lncSet)
		
		for lnc_subSet in lnc_subSets:
			if lnc_subSet not in lnc_subSets_explored:
				#shared strong targets of lncrnas in the subset
				targetSet = allExpTargets
				for lncrna in lnc_subSet:
					targetSet = targetSet & lncStrongTargets[lncrna]
					
				#exclude lncrna-target that have been explored
				targetSet_filter = set()
				for target in targetSet:
					not_explore = 0
					for lncrna in lnc_subSet:
						if (lncrna, target) not in lnc_target_explored:
							not_explore = 1
							break
					if not_explore:
						targetSet_filter.add(target)		
				
				if len(targetSet_filter) >= moduleSize:
					centerGene, module = find_center(targetSet_filter, corr_df, corrThres)

					#weak targets with affinity of >affiThres2 with all the lncrnas in the subset
					weak_targetSet = allExpTargets - targetSet
					for lncrna in lnc_subSet:
						weak_targetSet = weak_targetSet & (lncWeakTargets[lncrna] | lncStrongTargets[lncrna])
						
					#exclude lncrna-target that have been explored
					weak_targetSet_filter = set()
					for target in weak_targetSet:
						not_explore = 0
						for lncrna in lnc_subSet:
							if (lncrna, target) not in lnc_target_explored:
								not_explore = 1
								break
						if not_explore:
							weak_targetSet_filter.add(target)	
						
					#search for weak targets with correlation >corrThres with the center gene
					corr_weak_targetSet = set()
					for target in weak_targetSet_filter:
						if corr_df[centerGene][target] > corrThres:
							corr_weak_targetSet.add(target)

					#the weak targets with mean affinity >affiThres1 are added to the module
					for target in corr_weak_targetSet:
						if np.mean([exp_lnc_affinity[lncrna][target] for lncrna in lnc_subSet]) > affiThres1:
							module.add(target)
					
					#modules with size >=moduleSize are accepted
					if len(module) >= moduleSize:
						lncrna_toWrite = []
						for lncrna in lnc_subSet:
							func = determine_function(lncrna, module, corr_df)
							lncrna_toWrite.append(lncrna+'(%s)' % func)
						lncrna_toWrite_str = ','.join(lncrna_toWrite) 
						target_toWrite_str = ','.join(module)
						
						#top enriched pathways
						topFDR, topPVal, topPathway, topPathID, topPathCateg, hitGenes = pathway_analysis(module, 
																										  kegg_gene, 
																										  kegg_term, 
																										  kegg_link, 
																										  pathway_categ, 
																										  totalGeneNum)
						
						modules.append((lncrna_toWrite_str, 
									    target_toWrite_str, 
										topPathway, 
										topPathID, 
										topFDR,
										topPVal,
										topPathCateg, 
										hitGenes))
						
						for target in module:
							for lncrna in lnc_subSet:
								lnc_target_explored.add((lncrna, target))
						
				lnc_subSets_explored.add(lnc_subSet)
		
	#write results
	with open(moduleFile, 'w') as f:
		f.write('\t'.join([
						   'Module', 
						   'LncRNA(A=activator, I=inhibitor, N=not determined)',
						   'Target gene',
						   'Top enriched pathway',
						   'PathwayID',
						   'P value(hypergeometric distribution test)',
						   'FDR',
						   'Functional category',
                           'Hit genes'])
						   +'\n')
		moduleNum = 0
		for lncrna_toWrite_str, target_toWrite_str, topPathway, topPathID, topFDR, topPVal, topPathCateg, hitGenes in modules:
			moduleNum += 1
			moduleName = 'Module_'+str(moduleNum)
			f.write('\t'.join([
							   moduleName,
							   lncrna_toWrite_str, 
							   target_toWrite_str, 
							   topPathway, 
							   topPathID,
							   str(topPVal),
							   str(topFDR),
							   topPathCateg,
                               hitGenes])
							   +'\n')
							   
				
def generate_cytosc_file(moduleFile, cytoscapeDir):
	'''
	Generate the input files for module visulization using cytoscape.

	The files will be stored in the 'cytoscapeFile/' directory.
	'''
	if '/' in moduleFile:
		fileName = moduleFile.split('/')[-1]
	else:
		fileName = moduleFile
	category_show = set(['Metabolism',
						 'Immune system',
						 'Nervous system',
						 'NA'])
	f_re1 = open(cytoscapeDir+fileName+'_egde', 'w')
	f_re1.write('\t'.join(['lncRNA', 'module', 'weight'])+'\n')
	f_re2 = open(cytoscapeDir+fileName+'_group', 'w')
	f_re2.write('\t'.join(['module', 'category'])+'\n')
	allLncs = set()
	with open(moduleFile) as f_mod:
		f_mod.readline()
		for line in f_mod:
			cols = line.strip('\n').split('\t')
			moduleID = cols[0]
			lncrnas = cols[1]
			category = cols[7]
			pathway = cols[3]
			for lncrna in lncrnas.split(','):
				lncName = lncrna[:lncrna.find('(')]
				if '(A)' in lncrna:
					weight = '1'
				elif '(I)' in lncrna:
					weight = '-1'
				else:
					weight = 'NA'
				f_re1.write('\t'.join([lncName, moduleID, weight])+'\n')
				allLncs.add(lncName)
			if category in category_show:
				f_re2.write('\t'.join([moduleID, pathway])+'\n')
			else:
				f_re2.write('\t'.join([moduleID, 'Other'])+'\n')
		for lncName in allLncs:
			f_re2.write('\t'.join([lncName, 'LncRNA_regulator'])+'\n')
	f_re1.close()
	f_re2.close()    
					
					
if __name__ == '__main__':
	try:
		opts, args = getopt.getopt(sys.argv[1:],'',['h','t1=','t2=','m=','c=','f1=','f2=','f3=', 's=','o='])
	except getopt.GetoptError:
		print ('Usage: python eGRAM.py --t1 100 --t2 60 --m 5 --c 0.5 --f1 bindingMatrix --f2 expressionMatrix --f3 KeggCategory --s human --o output')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '--h':
			print ('Usage: python eGRAM.py --t1 100 --t2 60 --m 5 --c 0.5 --f1 bindingMatrix --f2 expressionMatrix --f3 KeggCategory --s human --o output')
			sys.exit(2)
		elif opt == '--t1':
			affiThres1 = int(arg)
		elif opt == '--t2':
			affiThres2 = int(arg)
		elif opt == '--m':
			moduleSize = int(arg)
		elif opt == '--c':
			corrThres = float(arg)
		elif opt == '--f1':
			bindFile = arg
		elif opt == '--f2':
			expFile = arg
		elif opt == '--f3':
			categoryFile = arg
		elif opt == '--s':
			species = arg.lower()
		elif opt == '--o':
			outFile = arg

	resultDir = 'moduleResults/'
	cytoscapeDir = 'cytoscapeFile/'
	
	if not os.path.exists(resultDir):
		os.mkdir(resultDir)
	if not os.path.exists(cytoscapeDir):
		os.mkdir(cytoscapeDir)
		
	moduleFile = resultDir + outFile
	
	if species == 'human':
		keggGeneFile = 'PathwayAnnotation/keggGene-hsa'
		keggPathFile = 'PathwayAnnotation/keggPathway-hsa'
		keggLinkFile = 'PathwayAnnotation/keggLink-hsa'
		gtfFile = 'Homo_sapiens.GRCh38.101.gtf'
	elif species == 'mouse':
		keggGeneFile = 'PathwayAnnotation/keggGene-mmu'
		keggPathFile = 'PathwayAnnotation/keggPathway-mmu'
		keggLinkFile = 'PathwayAnnotation/keggLink-mmu'
		gtfFile = 'Mus_musculus.GRCm38.101.gtf'

	kegg_gene, kegg_term, kegg_link, pathway_categ, totalGeneNum = read_kegg(keggGeneFile, keggPathFile, keggLinkFile, categoryFile)
	
	gene_symbols, protein_genes = read_gtf(gtfFile)
	
	affinity, targets = read_binding(bindFile, protein_genes)
	
	exp_df = read_exp(expFile, targets, affinity)
	print('Gene expression matrix:')
	print(exp_df)
	
	identify_module(exp_df, 
					affinity, 
					moduleFile,
					kegg_gene, 
					kegg_term, 
					kegg_link,
					pathway_categ,						
					totalGeneNum,
					categoryFile,
					affiThres1,
					affiThres2,
					moduleSize,
					corrThres)
					
	generate_cytosc_file(moduleFile, cytoscapeDir)