import pickle
from os.path import join
from os import listdir
import numpy as np
from matplotlib import pyplot
import math

#Plot the data obtained -- the linear score
#which_parameter: according to which annotation score
#input_set: the pickle you obtained from Annotator.py
def plot_linear(which_parameter,input_set):
	hscores = { 0: 'Gravy - Kyte & Doolttle (1982)',
		    1: 'Wimley & White (1996)',
	 	    2: 'Hessa et al. (2005)',
	            3: 'Eisenberg & McLachlan (1986)',
		    4: 'Black & Mould (1990)'}
	
		
	model = pickle.load(open(input_set))
	val = model[which_parameter]['hasa_cdr']
	
	#Plot the reference values
	#Plot the NGS reference data
	models = load_data('NGS','linear',which_parameter)
	mx_ngs = max(models)
	plot_models(models,'NGS','b')
	#Plot the therapeutic reference data.
	models = load_data('Therapeutic','linear',which_parameter)
	mx_ter = max(models)
	
	#Reference max value
	mx = max([mx_ter,mx_ngs])
	plot_models(models,'Therapeutics','g')

	pyplot.plot([val,val],[0,0.01],'k',linewidth=3,label='Your antibody')
	pyplot.xlabel('Linear H-ASA')
	pyplot.legend()
	pyplot.show()	
	#pyplot.savefig('/homes/krawczyk/Downloads/Cristian.png')

#Plot the data obtained -- the Patch score
#which_parameter: according to which annotation score
#input_set: the pickle you obtained from Annotator.py
def plot_patch(which_parameter,input_set):
	hscores = { 0: 'Gravy - Kyte & Doolttle (1982)',
		    1: 'Wimley & White (1996)',
	 	    2: 'Hessa et al. (2005)',
	            3: 'Eisenberg & McLachlan (1986)',
		    4: 'Black & Mould (1990)'}
	
		
	model = pickle.load(open(input_set))
	val = model[which_parameter]['adj_cdr']['hydrophobic']['hydrophobic']
	
	#Plot the reference values
	#Plot the NGS reference data
	models = load_data('NGS','patch',which_parameter)
	mx_ngs = max(models)
	plot_models(models,'patch','b')
	#Plot the therapeutic reference data.
	models = load_data('Therapeutic','patch',which_parameter)
	mx_ter = max(models)
	
	#Reference max value
	mx = max([mx_ter,mx_ngs])
	plot_models(models,'patch','g')

	pyplot.plot([val,val],[0,0.03],'k',linewidth=3,label='Your antibody')
	pyplot.legend()
	pyplot.show()

#Plot a single dataset given:
#models: the histogram data
#lab: the label for the legend
#fc: facecolor for the data.
def plot_models(models,lab,fc):
	pyplot.hist(models,bins=min(len(models),20), alpha=0.4, label=lab,facecolor=fc,normed=True)
	
def load_data(which_data,which_score,hydrophobicity):
	return pickle.load(open(join('plotting',which_data+'_'+which_score+'_'+str(hydrophobicity)+'.txt')))

#Get the reference data, without plotting.	
def get_patch_data_without_plotting(which_data,h_score,parameter,e1,e2):
	model_dir = join(which_data,str(h_score))#
	models = []
	for pdb in listdir(model_dir):
		
		try:
			model = pickle.load(open(join(model_dir,pdb)))
		except IOError:
			continue
		
		#What is the entry we are after.
		model_data = model[parameter][e1][e2]

		models.append(model_data)
	return models

#Get the reference data, without plotting.	
def get_linear_data_without_plotting(which_data,h_score,parameter):
	model_dir = join(which_data,str(h_score))#
	models = []
	for pdb in listdir(model_dir):
		
		try:
			model = pickle.load(open(join(model_dir,pdb)))
		except IOError:
			continue
		
		#What is the entry we are after.
		model_data = model[parameter]
	
		models.append(model_data)
	return models

if __name__ == '__main__':

	import sys
	#USAGE
	#python PlotResults.py [results_file] [linear/patch] [hydrophobicity parameter]
	#e.g. python PlotResults.py results.txt linear 0
	hydrophobicity_score = int(sys.argv[3])
	print "Hydrophobicity score",hydrophobicity_score
	print "Performing",sys.argv[1]
	if sys.argv[2] == 'linear':
		print "Plotting linear"
		plot_linear(hydrophobicity_score,sys.argv[1])
	if sys.argv[2] == 'patch':
		plot_patch(hydrophobicity_score,sys.argv[1])
