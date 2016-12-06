# -*- coding: utf-8 -*-
from os import listdir
import os,sys,tempfile
from os.path import isfile, join
from PDBUtils import PDBchain
import pickle
import numpy as np

#remove the insertion codes from Chothia codes
def stripChothia(res):
	icode = "ZXCVBNMLKJHGFDSAWQERTYUIOP"
	if res[len(res)-1] in icode:
		res = res[0:(len(res)-1)]
	return res

#Check if a residue is a CDR residue according to the given definition (deff)
def is_CDR(res,deff):
	#return "H3"
	_res = stripChothia(res)
	for CDR in definitions[deff]:
		if _res in definitions[deff][CDR]:
			return CDR
	return False

# Check if a residue is at the variable-constant interface#
def is_not_VC_interface(res,deff):
	if deff != 'chothia': raise NotImplementedError("Check for non-VC-interface currently only implemented using chothia numbering. Requested numbering: '%s'"%deff)
	_res = stripChothia(res)
	return _res not in definition_chothia_VCinterface

# This was defined by taking all residues, from a set of 10 antibody Fab models with varying sequence and affinity, which were within 8 angstroms of any constant chain atom.
definition_chothia_VCinterface = ['H10', 'H107', 'H108', 'H109', 'H11', 'H110', 'H111', 'H112', 'H113', 'H12', 'H13', 'H14', 'H15', 'H41', 'H42', 'H7', 'H8', 'H83', 'H84', 'H85', 'H87', 'H88', 'H89', 'H9', 'L10', 'L103', 'L104', 'L105', 'L106', 'L107', 'L108', 'L109', 'L11', 'L110', 'L12', 'L13', 'L14', 'L15', 'L39', 'L40', 'L41', 'L80', 'L81', 'L82', 'L83', 'L84']

##MAIN##

#AA mapping
Naive_HH = {
'ALA': 'B',
'ARG': 'C',
'ASN': 'U',
'ASP': 'C',
'CYS': 'S',
'GLU': 'C',
'GLN': 'U',
'GLY': 'S',
'HIS': 'C',
'ILE': 'H',
'LEU': 'H',
'LYS': 'C',
'MET': 'B',
'PHE':	'H',
'PRO': 'S',
'SER': 'U',
'THR': 'U',
'TRP': 'H',
'TYR': 'H',
'VAL': 'H'}

coloring = {'H':100,'C':25,'U':50,'S':75}


#AA mapping
AAtable = {
'ALA': 'A',
'ARG': 'R',
'ASN': 'N',
'ASP': 'D',
'CYS': 'C',
'GLU': 'E',
'GLN': 'Q',
'GLY': 'G',
'HIS': 'H',
'ILE': 'I',
'LEU': 'L',
'LYS': 'K',
'MET': 'M',
'PHE': 'F',
'PRO': 'P',
'SER': 'S',
'THR': 'T',
'TRP': 'W',
'TYR': 'Y',
'VAL': 'V'}

#CDR Definitions - Chothia+2 on either side.
definitions = {		 
		
		"chothia" : {"L1" : ["L22","L23","L24", "L25", "L26", "L27", "L28", "L29", "L30", "L31", "L32", "L33", "L34","L35","L36"], 
					"L2" : ["L48","L49","L50", "L51", "L52", "L53", "L54", "L55", "L56","L57","L58"],
					"L3" : ["L87","L88","L89", "L90", "L91", "L92", "L93", "L94", "L95", "L96", "L97","L98","L99"], 
					"H1" : ["H24","H25","H26", "H27", "H28", "H29", "H30", "H31", "H32","H33","H34"], 
					"H2" : ["H50","H51","H52", "H53", "H54", "H55", "H56","H57","H58"] ,
					"H3" : ["H93","H94","H95", "H96", "H97", "H98", "H99", "H100", "H101", "H102","H103","H104"]} 
		
		}

#Get the total area of the side-chain
def parsePSAArea(psafile):
	dicto = dict()
	curr_chain = 'H'
	dicto['H'] = dict()
	dicto['L'] = dict()
	last = -1
	for line in open(psafile):
			line =	line.strip()
			if line[0:6] =='ACCESS':
				sid_int = int(line[6:11])
				if sid_int<last:
					curr_chain = 'L'
				last = sid_int
				sid = line[6:12].replace(" ","")
				dicto[curr_chain][sid] = float(line[55:61])
				
	return dicto

#Parses the PSA output - it assumes that the chains come in heavy light order...
def parsePSA(psafile):
	dicto = dict()
	curr_chain = 'H'
	dicto['H'] = dict()
	dicto['L'] = dict()
	last = -1
	for line in open(psafile):
			line =	line.strip()
			if line[0:6] =='ACCESS':
				sid_int = int(line[6:11])
				if sid_int<last:
					curr_chain = 'L'
				last = sid_int
				sid = line[6:12].replace(" ","")
				dicto[curr_chain][sid] = float(line[61:67])
				#float(line[61:67])
				
	return dicto

def runPSA(pdb_file,where_to):
	dobavka = ""
	if (sys.platform=='darwin'):
		dobavka = "_mac"
	os.system(local_path+"/aux/./psa"+dobavka+" -t "+pdb_file+" > "+where_to)

#Create a temporary folder in the user directory - return it so that we can do whatever with it
def create_temp_folder():
	dirpath = tempfile.mkdtemp()
	return dirpath

#Normalize the hydrophobicity to a z-score
def normalize(matrix,index):
	#find max-min
	vals = []
	for res in matrix:
		vals.append(matrix[res][index])
	annotation = dict()
	mu = np.mean(vals)
	sigma = np.std(vals)
	
	for res in matrix:
		#annotation[res] = (matrix[res][index]-mu)/sigma
		#Normalizing to between 1 and 2
		annotation[res]= ((matrix[res][index]+(0-min(vals)))/((max(vals)-min(vals))))+1.0
	
	return annotation


#For simplicity, we are not doing charge at this point
def read_charge(f_in):
	nm = f_in.replace('.pdb','')
	ch_map = dict()
	pka_map = dict()
	pka_max = -10000
	pka_min = 1000000
	ch_max = -10000
	ch_min = 1000000
	
	#Get the individual residue maps
	#DIFFERENT PH	
	for line in open(join('surfacech','results',nm,nm+'-pH05.00.tsv')):
		if 'NAME' not in line:
			line = line.split('\t')
			r_id = line[2]+line[1]
			ch_map[r_id] = line[4]
			pka_map[r_id] = line[3]
			pka = float(line[3])
			ch = float(line[4])
			if pka > pka_max:
				pka_max = pka
			if pka < pka_min:
				pka_min = pka
			if ch > ch_max:
				ch_max = ch
			if ch < ch_min:
				ch_min = ch
					
	summary_charge = dict()
	#Get the summary
	for line in open(join('surfacech','summaries',nm)):
		if 'pH' in line:
			line = line.strip()
			line = line.replace("# pH:","")
			line = line.replace("Net charge:","")
			line = line.split(';')
			ph = float(line[0])
			ch=  float(line[1])
			summary_charge[ph] = ch
	return ch_map,pka_map,summary_charge,pka_max,pka_min,ch_max,ch_min

#Which of the four classes do we belong to: positive, negative, hydrophobic or none
residue_classes = ['hydrophobic','positive','negative']

def CreateAnnotation(annotation_index,which_ph,input_file,output_fold):

	results_folder = join(local_path,output_fold,str(annotation_index))
	
	#Initialize the hydrophobic mappings
	hydrophobics = dict()
	for line in open('Hydrophobics.txt'):
		line = line.strip().split('\t')
		hydrophobics[line[0].upper()] = (float(line[1]),float(line[2]),float(line[3]),float(line[4]),float(line[5]))
	#Normalize the hydrophobics
	hydrophobicity_annotation = normalize(hydrophobics,annotation_index)
	#Where the input files are. TODO parametrize.
	mypath = join(local_path,input_fold)
	onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) and f.endswith(".pdb") ]
	done = 0
	stats = dict()
	
	if not os.path.exists(results_folder):
		os.mkdir(results_folder)
	if not os.path.exists(local_path+'/contact_maps'):
		os.mkdir(local_path+'/contact_maps')
	if not os.path.exists(local_path+'/surfaces'):
		os.mkdir(local_path+'/surfaces')
	for f_in in onlyfiles:

		
		done+=1
		
		print done,'/',len(onlyfiles)
		print f_in
		stats[f_in] = dict()
		#Get surface exposed residues by running PSA
		runPSA( join(local_path,input_fold,f_in),'surfaces/'+f_in)
		asa = parsePSA('surfaces/'+f_in)
		asa_area = parsePSAArea('surfaces/'+f_in)
		pdbchain = PDBchain( join(local_path,input_fold,f_in),'H')
		#READ CHARGE DATA
		try:
			ch_map,pka_map,summary_charge,pka_max,pka_min,ch_max,ch_min = read_charge(f_in)
		except IOError:
			#TODO background issue
			print "IOError",f_in
			
			continue		
		#Get the residue objects
		ab_structure = {'H':PDBchain( join(local_path,input_fold,f_in),'H'),'L':PDBchain( join(local_path,input_fold,f_in),'L')}

		residues = dict()
		for chain in ['H','L']:
			for elem in sorted(ab_structure[chain].residues):

				#Ignore entries which are on the VC interface.
				if is_not_VC_interface(chain+str(elem[0]),'chothia')==False:
					#print chain+str(elem[0]),'on the vc interface, ignoring...'
					continue

				#Get the residue type here.
				typ = ab_structure[chain].residues[elem].res_type
				r_elem = {'res':ab_structure[chain].residues[elem],'typ':typ,'charged':None,'hydrophobic':None,'asa':None,'asa_sum':None}
				bfact = 0
				res_asa = 0
				res_asa_area = 0
				try:
					res_asa = asa[chain][str(elem[0])+elem[1]]
					res_asa_area = asa_area[chain][str(elem[0])+elem[1]]
				except KeyError:
				#	#TODO Sometimes ASA skips residues -- especially the terminal ones.
					print "ASA ERROR",f_in
					continue
				if res_asa< 7.5:#Do not bother with buried residues.
					continue
				#print elem,res_asa, "Residue is exposed, continuing"
				r_elem['asa'] = res_asa
				r_elem['asa_sum'] = res_asa_area
				#Do not bother if the entries are not CDRs.
				#TODO: parametrize on whether we want to look at CDR data.			
				#if is_CDR(chain+str(elem[0]),'chothia')== False:
					#for atom in ab_structure[chain].residues[elem].atoms:
						#f_out.write(pdbchain.color_line(atom.pdb_line, str(bfact)+'.00')+'\n')
				#	continue	
				#Hydrophobicity value
				r_elem['hydrophobic'] = float(hydrophobicity_annotation[typ])
		
				#Charge
				#If a residue can be associated with charge-- make it charged. Everything else is just different shades of hydrophbic
				chentry = chain+str(elem[0])+str(elem[1])
				if chentry in ch_map:
					#print  'Charge', ch_map[chentry]
					charge = float(ch_map[chentry])
					if charge>0:
						r_elem['positive'] = charge
					else:
						r_elem['negative'] = charge
					
				
				residues[(chain,elem)] = r_elem
		
		#Define CDR neighborhood- so as not to calculate CDRs only...
		cdr_vicinity = []

		#Calculate cdrs and their neighbors.
		number_cdr_residues = 0
		for r1 in residues:
			chid = r1[0]+str(r1[1][0])
			if is_CDR(chid,'chothia')==False:
				continue
			else:
				number_cdr_residues+=1
				#Add the cdrs in naturally
				if r1 not in cdr_vicinity:
					cdr_vicinity.append(r1)

			for r2 in residues:
				#Skip if we are the same residue or if the other is a CDR as well
				if r1==r2 or r2 in cdr_vicinity:
					continue
				#If we are here, we know r1 is CDR and r2 is not. See if r2 is in the neighborhood then?
				res1 = residues[r1]['res']
				res2 = residues[r2]['res']
				d = res2.distance(res1,res2)
				if d<4.0:#Let's say 4.0 is close enough to be in the CDR neighborhood.
					if r2 not in cdr_vicinity:
						cdr_vicinity.append(r2)
						break
		print "We have",len(cdr_vicinity),"in CDR area of which",number_cdr_residues,'are cdr'
		print "There are total of ",len(residues),'residues that we take into account'
		
		
		
		#We will collect the following data:
		#Number of CDR residues
		#Number of residues in CDR vicinity
		#Number of residues on surface
		#Surface Area- entire mol [v]
		#Surface Area, hydrophobicity-weighted- entire mol [v]
		#Summary Charge [v]
		#Positive charge sum- entire mol
		#Negative charge sum- entire mol	
		#Adjecency matrix- entire molecule
		print "Summary charge",summary_charge[7.5]
		#Calculate the total surface area
		total_area = 0.0
		#hydrophobicity-scaled area
		h_asa_total = 0.0
		#positive charge sum
		pos_charge_total = 0.0
		#negative charge sum
		neg_charge_total = 0.0
		for r1 in residues:
			total_area+=residues[r1]['asa_sum']
			h_asa_total+=residues[r1]['asa_sum']*residues[r1]['hydrophobic']
			if 'negative' in residues[r1]:
				neg_charge_total+=residues[r1]['negative']
			if 'positive' in residues[r1]:
				pos_charge_total+=residues[r1]['positive']
		#Calculate CDR only area
		cdr_area = 0.0
		#Calcualte the CDR only hasa
		h_asa_cdr = 0.0
		#positive charge sum
		pos_charge_cdr = 0.0
		#negative charge sum
		neg_charge_cdr = 0.0
		#Only look at the residues in the immediate vicinity of the CDRs
		for r1 in cdr_vicinity:
			
			cdr_area+=residues[r1]['asa_sum']
			h_asa_cdr+=residues[r1]['asa_sum']*residues[r1]['hydrophobic']
			if 'negative' in residues[r1]:
				neg_charge_cdr+=residues[r1]['negative']
			if 'positive' in residues[r1]:
				pos_charge_cdr+=residues[r1]['positive']
		print "Total Area",total_area
		print "CDR area",cdr_area
		print "Total hASA",h_asa_total
		print "CDR hASA",h_asa_cdr
		print "Negative charge, total",neg_charge_total
		print "Negative charge, CDR",neg_charge_cdr
		print "Positive charge, total",pos_charge_total
		print "Positive charge, CDR",pos_charge_cdr
		
		#PREPARE OUTPUT DATA- ENTIRE MOLECULE
		#ADJECENCY MATRIX FOR THE ENTIRE MOLECULE
		adj_matrix = dict()
		for i in residue_classes:
			adj_matrix[i] = dict()
			for j in residue_classes:
				adj_matrix[i][j] = 0

		for r1 in residues:
			for r2 in residues:
				if r1 != r2:
					#Get the distance
					res1 = residues[r1]['res']
					res2 = residues[r2]['res']
					d = res2.distance(res1,res2)
					if d>7.5:
						continue
					for c1 in residue_classes:	
						for c2 in residue_classes:
							if c1 in residues[r1] and c2 in residues[r2]:
								#The Coulombic measure.
								
								coulombic_score = (residues[r1][c1]*residues[r2][c2])/(d*d)
								adj_matrix[c1][c2] +=coulombic_score
		#ADJECENCY MATRIX FOR THECDRs only
		adj_matrix_cdr = dict()
		for i in residue_classes:
			adj_matrix_cdr[i] = dict()
			for j in residue_classes:
				adj_matrix_cdr[i][j] = 0

		for r1 in cdr_vicinity:
			for r2 in cdr_vicinity:
				if r1 != r2:
					#Get the distance
					res1 = residues[r1]['res']
					res2 = residues[r2]['res']
					d = res2.distance(res1,res2)
					if d>7.5:
						continue
					for c1 in residue_classes:	
						for c2 in residue_classes:
							if c1 in residues[r1] and c2 in residues[r2]:
								#The Coulombic measure.
								coulombic_score = (residues[r1][c1]*residues[r2][c2])/(d*d)
								adj_matrix_cdr[c1][c2] +=coulombic_score
								

								
		print "Adjecency matrix total:",adj_matrix
		print "Adjecency matrix CDR:",adj_matrix_cdr
		#Full statistics used for calculations.		
		stats = dict()
		stats['adj_cdr'] = adj_matrix_cdr
		stats['adj_tot'] = adj_matrix
		stats['area_tot'] = total_area
		stats['cdr_area'] = cdr_area
		stats['hasa_tot'] = h_asa_total
		stats['hasa_cdr'] = h_asa_cdr		
		stats['neg_tot'] = neg_charge_total
		stats['neg_cdr'] = neg_charge_cdr
		stats['pos_tot'] = pos_charge_total
		stats['pos_cdr'] = pos_charge_cdr
		stats['total_exposed'] = len(residues)
		stats['num_cdrs'] = number_cdr_residues
		stats['tot_vicinity'] = len(cdr_vicinity)
		stats['total_charge'] = summary_charge[7.5]
		
 		pickle.dump(stats,open(join(results_folder,f_in),'w'))
		#f_out.close()
			
		#quit()
		
#MAIN#
if __name__ == '__main__':

	#0: Gravy - Kyte & Doolttle (1982)
	#1: Wimley & White (1996)
	#2: Hessa et al. (2005)
	#3: Eisenberg & McLachlan (1986)
	#4: Black & Mould (1990)

	#Get the local path.
	local_path = os.path.dirname(os.path.realpath(__file__))
	#Which ph value to use.
	#We used to run the charge calculations on 	
	#pH = sys.argv[1]
	input_fold = sys.argv[2]
	output_fold = sys.argv[3]
	annotation =0
	CreateAnnotation(annotation,float(pH),input_fold,output_fold)
