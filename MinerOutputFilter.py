#!/usr/bin/env python2

from __future__ import division
import argparse
import os
import shutil
from Bio import AlignIO
from itertools import cycle, islice
import itertools
import re

# regex ecpressions to identify and exclude refernce sequences in prep_fastas function
REFERENCES = {'Alyrata': 'AL[0-9]G[0-9]{5}',
             'Athaliana': 'AT[0-9]G[0-9]{5}',
             'Bdistachyon': 'BD[0-9]G[0-9]{5}',
             'Cpapaya': 'CP[0-9]{5}G[0-9]{5}',
             'Fvesca': 'FV[0-9]G[0-9]{5}', 
             'Gmax': 'GM[0-9]{2}G[0-9]{5}',
             'Mdomestica': 'MD[0-9]{2}G[0-9]{6}',
             'Mesculenta': 'ME[0-9]{5}G[0-9]{5}',
             'Mtrunculata': 'MT[0-9]G[0-9]{6}', 
             'Osativa': 'OS[0-9]{2}G[0-9]{5}',
             'Ptrichocarpa': 'PT[0-9]{2}G[0-9]{5}',
             'Rcommunis': 'RC[0-9]{5}G[0-9]{5}',
             'Sbicolor': 'SB[0-9]{2}G[0-9]{6}',
             'Tcacao': 'TC[0-9]{2}G[0-9]{6}',
             'Vvinifera': 'VV[0-9]{2}G[0-9]{5}',
             'Zmays': 'ZM[0-9]{2}G[0-9]{5}'}

AVAILLABLE = [k for k in REFERENCES]
#############################
#functions

def roundrobin(*iterables):
    #roundrobin('ABC', 'D', 'EF') --> A D E B F C
    # Recipe credited to George Sakkis
    pending = len(iterables)
    nexts = cycle(iter(it).next for it in iterables)
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1
            nexts = cycle(islice(nexts, pending))
			
def reduce_coordinates(coordinates):
    # takes a list of coordinates (form exon1_start, intron_start, ex...) and ...
	# ... excludes all regions with "gaps" shorter than 37 bp, as these are not true introns
    # the following line will split the list of coordinates into a list of tuples of two numbers
    zipped = zip(coordinates[0::2], coordinates[1::2])
    to_pop = []
    for pos in range(len(zipped)-1):
        # this part loops through the list of tuples, converts each tuple and the immediately following one
        #... back into a list and compares the last element of list 1 with the first element of list 2 
        #... these are the introns: if the result is < 37 both numbers will be deleted from the original list
        a = list(zipped[pos])
        b = list(zipped[pos +1]) 
        if b[0] - a[1] < 37:
            to_pop.append(a[1]) 
            to_pop.append(b[0])
    for item in range(len(to_pop)):
        if to_pop[item] in coordinates:
            coordinates.remove(to_pop[item])
    # if the seq ends on a long intron, the list will have an uneven number of entries,.... 
    #.....with the last number representing...the end of the an intron that cannot be removed;
    #... the next part will remove it
    if len(coordinates)%2 != 0:
        del coordinates[-1]
    return coordinates
	

def define_lentgh(positions, minlen, maxlen):
    # this function excudes all exons exceeding a maximum and minimum length
    zipped = zip(positions[0::2], positions[1::2])
    to_pop = []
    for pos in range(len(zipped)):
        lst = list(zipped[pos])
        if lst[1] - lst[0] < minlen + 1 or lst[1] - lst[0] > maxlen + 1:
            to_pop.append(lst[0])
            to_pop.append(lst[1])
    for item in range(len(to_pop)):
        if to_pop[item] in positions:
            positions.remove(to_pop[item])
    return positions


def get_differences(seqs):
    positions = []
    for subset in itertools.combinations(seqs, 2):
        for base in range(len(subset[0])-1):
            if subset[0][base] == subset[1][base]:
                 pass
            else:
                positions.append(base)     
    return len(set(positions))*100/len(seqs[0])   


def create_geneList(MinerDir, n):
    # creates List of gene_names supported by at least n input transcriptomes
    # the number of transcriptomes is in the single_copy_gene_txt file as last digit before the last ';'
    # initialise empty list
    SCG_list = []
    with open(os.path.join(MinerDir, 'single_copy_genes.txt')) as f:
        first_line = f.readline()
        for line in f:
        # split line at ; and assort 1st part (gene, function etc) to key and 2nd (Datasets in which hit occurs) to value
            try:
                key, value = line.split(';')
            except:
                new_line = line.replace(';', ',', line.count(';') - 1)
                key, value = new_line.split(';')
            # check number of transcriptomes, and append orthologous gene name to SCG_list
            if int(value.split()[0]) >= n :
                SCG_list.append(key.split('\t')[0])
                
    # sort the list numerically
    sorted_SCG_list = sorted(SCG_list)
    print str(len(sorted_SCG_list)) + ' reference loci supported by at least ' + str(n) + ' datasets'

    with open(os.path.join(MinerDir, 'GeneList.txt'), 'w') as genes:
        for item in sorted_SCG_list:
            genes.write("%s\n" % item)
    return sorted_SCG_list
	
def prep_fastas(MakerMinerDir, selected_dir, genelist, regex):
# this code will read the MarkerMiner output files and write a new non-interleaved file  without the Atha sequence...
# it is a prerequisite to the counts function
	fastafolder = os.path.join(MakerMinerDir,'MAFFT_ADD_REF_ALIGN_FASTA')
	for element in genelist:
		for f in os.listdir(fastafolder):
			if element == f.split('.')[0]:
				dest = os.path.join(selected_dir, element + '.ingroup.fasta')
				with open(fastafolder + '\\' + f) as infi:
					with open(dest, 'w+') as newfasta:
						for record in AlignIO.read(infi, "fasta"):
							if re.search(regex, record.id):
								pass
							else:
								sequence = str(record.seq)
								newfasta.write('>' + record.id + '\n' + sequence + '\n')	
								
def find_exons(fastas, outdir, minlen, maxlen):
	# this code will read every non-interleaved fasta file and identify the exonic regions
	# the coordinates (0-based indices) will be written to a csv file for every gene in the new Counts directory
	# the output numbers are the indices of the bases in the string seqeunce; 
	# the actual position of each base is the index + 1!
	for f in os.listdir(fastas):
		with open(os.path.join(fastas, f), 'r') as fasta:
			for line in fasta:
				start_exon = []
				end_exon = []
				if line.startswith('>'):
					header = line[1:].rstrip()
				else :
					if line[1] != '-':
						start_exon.append(1)
					else:
						end_exon.append(1)
					for base in range(len(line) -1):
						if line[base -2: base] == '--' and line[base] != '-' :
								start_exon.append(base)
						elif line[base + 1 : base + 3] == '--' and line[base] != '-' :
								end_exon.append(base +1)
					if line[(len(line) - 1)] != '-':
							start_exon.append(len(line) - 1)
					else:
							end_exon.append(len(line) - 1)
					if start_exon[0] == 1 :
						result = list(roundrobin(start_exon, end_exon))
					else: 
						result = list(roundrobin(start_exon, end_exon[1:]))
					min_intron_len = reduce_coordinates(result)
					min_exon_len = define_lentgh(min_intron_len, minlen, maxlen)
					final = zip(min_exon_len[0::2], min_exon_len[1::2])
					if len(final) > 0:
						with open(os.path.join(outdir, f.split('.')[0] + '_count.csv'), 'a') as csv:
							for element in final:
								csv.write('{},{},{}\n'.format(f.split('.')[0], element[0], element[1]))
					else:
						continue
	print str(len(os.listdir(outdir))) + ' exons fit length criteria'

def pop_doublettes(csv):					
	# the next three loops read in the csv file and exclude all exons, that are still present as 
	# double entries due to one record being shorter than the other in the original alignment
			#initialise empty dicts
	end_doubles = {}
	final = {}
	start_doubles = {}
	for line in open(csv, 'r'):
		start = int(line.split(',')[1])
		end = int(line.split(',')[2].rstrip())
		if start not in start_doubles:
			start_doubles[start] = end
		else:
			for key, v in start_doubles.items():
				if end - v < 0:
					start_doubles[start] = end

	for start, end in start_doubles.items():
		if end not in end_doubles.keys():
			end_doubles[end] = start
		else:
			if start > end_doubles[end]:
				end_doubles[end] = start
				
	del start_doubles
	
	for k, v in end_doubles.items():
		final[v] = k
		
	del end_doubles
	
	return final

def get_variability(Listsdir, fastas_dir, maxvar, minvar, Probesdir):
	for csv in os.listdir(Listsdir):		
		for fasta in os.listdir(fastas_dir):
			if csv.split('_')[0] == fasta.split('.')[0]:
				gene = fasta.split('.')[0]
				alignment = AlignIO.read(os.path.join(fastas_dir, fasta), 'fasta')
				exons = pop_doublettes(os.path.join(Listsdir, csv))
				for start, end in exons.items():
					subaln = [str(record.seq[start: end]) for record in alignment]
					diff = get_differences(subaln)
					if diff <= maxvar and diff >= minvar:
						probe = os.path.join(Probesdir, str(diff)[:4] + '_' + gene + '_' + str(start) + '.fasta')
						with open(probe, 'a') as p:
							with open(os.path.join(fastas_dir, fasta)) as fas:
								for l in fas:
									if l[0] == '>':
										p.write('{}_{}_{}\n'.format(l[:5], gene, str(start)))
									else:
										p.write(l[start:end] + '\n')
	print str(len(os.listdir(Probesdir))) + ' exons in final Probe set.'									
								
								
								
#################################################################


def main():
#parse arguments
##############################
	parser = argparse.ArgumentParser(description='filter MarkerMiner (Chamala et al., 2015) output for exon length and variability to create fasta files with probes')
	parser.add_argument('MarkerMinerOutputDir', type=str, help='full path to extracted MarkerMiner output directory')
	parser.add_argument('--reference', type=str, help='reference used in MarkerMiner analysis. ', required=True, choices=AVAILLABLE)
	parser.add_argument('--minlength', type=int, help='minimum required exon length. Default: 120', default=120)
	parser.add_argument('--maxlength', type=int, help='maximum allowed exon length. Default: 1000', default=1000)
	parser.add_argument('--maxvar', type=float, help='maximum allowed exon variability (discarding reference). Default: 10', default=10)
	parser.add_argument('--minvar', type=float, help='minimum required exon variability (discarding reference). Default: 0', default=0)
	parser.add_argument('--cleanup', help='if selected, all intermediates files will be deleted', default=False, action='store_true')
	
	args = parser.parse_args()
#############################	
	
	SCGs = os.path.join(args.MarkerMinerOutputDir, 'single_copy_genes.txt')
	selected_fastas = os.path.join(args.MarkerMinerOutputDir, 'SelectedFastas')
	ExonsDir = os.path.join(args.MarkerMinerOutputDir, 'ExonLists')
	ProbesDir = os.path.join(args.MarkerMinerOutputDir, 'Probes_' + str(args.minvar) + '-' + str(args.maxvar))

	folders = [selected_fastas, ExonsDir, ProbesDir]
	for new_dir in folders:
		if not os.path.exists(new_dir):
			os.mkdir(new_dir)

	SCGlist = create_geneList(args.MarkerMinerOutputDir, 2)

	prep_fastas(args.MarkerMinerOutputDir, selected_fastas, SCGlist, REFERENCES[args.reference])

	find_exons(selected_fastas, ExonsDir, args.minlength, args.maxlength)

	get_variability(ExonsDir, selected_fastas, args.maxvar, args.minvar, ProbesDir)
	
	if args.cleanup == True:
		os.remove(os.join.path(MarkerMinerOutputDir, 'GeneList.txt'))
		shutil.rmtree(ExonsDir, ignore_errors=True)
		shutil.rmtree(selected_fastas, ignore_errors=True)

if __name__ == '__main__':
	main()
