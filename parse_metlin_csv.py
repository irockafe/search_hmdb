from bs4 import BeautifulSoup
import csv
from sys import argv
import pandas
#BeautifulSoup also requires libxml-dev and libxslt-dev packages!

'''
This code takes the csv file downloaded from at 
METLIN batch search. It will search HMDB for METLIN hits
and add information, like if the metabolite is from urine or 
feces, to a new csv file

Metlin does not output masses that where nothing was found.
So, users should add their input masses and I can infer that
nothing was found if the METLIN csv has nothing in it for that mass
'''
#You have to point python to your metlin csv in the
#command line
#script, metlin_csv = argv

#Path to csv file
#bigger toy metlin_csv = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/metlin_toy_data.csv'
metlin_csv = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/smaller_toy_data.csv'

hmdb_all_metabolites = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/hmdb_metabolites/hmdb_metabolites.xml'
toy_hmdb_all_metabolites = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/toy_hmdb/toy_hmdb.xml'

hmdb_directory = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/hmdb_metabolites/'
toy_hmdb_directory = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/toy_hmdb/'
output_filename = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/toy_hmdb/HOLYSHIT_IT_WORKED.csv'
#import csv to pandas, with first row as the header
metlin_data = pandas.read_csv(metlin_csv, header=0)
masses = metlin_data['mass']
#Count number of repeat input mass values - the number of isomers
isomers = metlin_data['inputmass'].value_counts()

#These are the masses you should search HMDB for
#We're not searching for individual compounds, instead
#we want to find anything that can match our mass
#masses_to_search = metlin_data['mass'].value_counts()
masses = metlin_data['mass']
#mol_ids = metlin_data['molid']
#adduct_masses_df = metlin_data[['inputmass', 'mass']]
#Convert dataframe into dictionary object
#molid_masses_dict = molid_masses_df.set_index('inputmass').T.to_dict('list')

toy_masses = [102.031694058, 74.08439833]
metabolites = BeautifulSoup(open(hmdb_all_metabolites), "lxml").find_all('metabolite')
#save memory by not opening two giant files
#hmdb_soup = BeautifulSoup(open(hmdb_all_metabolites), "lxml")
#metabolites = hmdb_soup.find_all('metabolite')

def search_HMDB(masses, metabolites):
	'''
	INPUT - a dictionary of {mass: molid} of the masses METLIN identified as 
	compounds and their molids, and an xml tree of metabolites

	TODO: Modify code to work on a single metabolite

	FUNCTION - parse a downloaded HMDB directory 
	OUTPUT - a dictionary of masses and HMDB ID numbers that match
	this mass ex: {mass: [HMDB_ID_1, HMDB_ID_2]}
	 These will be used to access individual HMDB_ID 
	files to make searching faster (??)
	'''
	print 'Beginning to search HMDB'
	hmdb_matches = {}
	for metabolite in metabolites:
		MW = float(metabolite.monisotopic_moleculate_weight.contents[0])
		#TODO - may need to not do exact matching, but define a 0.1ppm 
		#threshold to determine if MW is one of our target masses.
		if MW in masses:
			print 'Found a match for mass %s: %s!' % (MW, metabolite)
			#unicode string
			hmdb_id = metabolite.accession.contents[0]
			#Add {MW:id} to output dictionary
			if MW in hmdb_matches:
				hmdb_matches[MW].append(hmdb_id)
			else:
				hmdb_matches[MW] = [hmdb_id]

	return hmdb_matches

def get_hmdb_info(hmdb_identified):
	'''
	INPUT - dictionary of {mass1:[HMDB_ID_1, HMDB_ID_2, ...], 
						   mass2:...}
	FUNCTION - Open the correct HMDB metabolite's file and get
		information that you want from it. i.e. biofluid, 
		disease type, etc.
	OUTPUT - {Mass_1: {HMDB_ID_1: {Urine: 0, Feces: 0, Others}, 
			  		   HMDB_ID_2: ...},
			  Mass_2: ...
			  }
	'''
	print 'Finished Identifying Potential Compounds.\n Moving on to gathering compound information.'
	output = {}
	for mass in hmdb_identified.keys():
		output[mass] = {}
		#print mass
		for entry in hmdb_identified[mass]:
			#print entry
			with open(hmdb_directory+entry+'.xml') as metabolite:
				#Get biofluid status
				metabolite_soup = BeautifulSoup(metabolite)
				biofluid_info = which_biofluids(metabolite_soup)
			output[mass][entry] = biofluid_info
	return output

def which_biofluids(hmdb_metabolite):
	'''
	INPUT -	A BeautifulSoup object from an individual HMDB accession number's xml file
	FUNCTION - gets the various biofluid locations from the xml
	OUTPUT - {Name: compound_name, Urine: 0, Feces:1, Other_fluids: 'string'}
	'''
	compound_name = hmdb_metabolite.find('name').contents[0]
	#print compound_name
	biofluids_present = {'Name': compound_name, 'Urine': 0, 'Feces': 0, 'Other Biofluids': 'None Listed'}
	Other_fluids = []
	#Get all the biofluid locations listed and run
	#through them, looking especially for Urine and Feces
	biofluid_locations = hmdb_metabolite.biofluid_locations.find_all('biofluid')
	for biofluid in biofluid_locations:
		#print biofluid.contents
		if biofluid.contents[0] == 'Urine':
			biofluids_present['Urine'] = 1
			break
		elif biofluid.contents[0] == 'Feces':
			biofluids_present['Feces'] = 1
			break
		else:
			Other_fluids.append(biofluid.contents[0])
	if len(Other_fluids) > 0:
		#Turn list into string
		fluids = ', '.join(Other_fluids)
		biofluids_present['Other Biofluids'] = fluids

	return biofluids_present

def combine_info(biofluids, metlin_dataframe, isomers):
	'''
	INPUT - Dataframe of HMDB, Urine, feces, and mass info
	FUNCTION - Merge the metlin mass, adduct, and dppm with
		# isomers, Urine, Feces, and Other biofluid presence
	OUTPUT - Pandas dataframe with [mass, adduct, ppm, Urine,
		Feces, Other] as headers.
	'''
	#make the mass values indices for the dataframe
	metlin_dataframe = metlin_dataframe.set_index('mass')
	curated_df = pandas.concat([metlin_dataframe['adduct'],
								metlin_dataframe['dppm']],
								axis=1)
	curated_df['Isomers'] = pandas.Series(index=curated_df.index)
	#If isomer mass is equal to row mass, isomer category
	# in dataframe = mass from isomers

	#Add the number of isomers to dataframe
	#print curated_df.index
	curated_df = add_isomers_to_data(curated_df, isomers)
	curated_df = curated_df.drop_duplicates()

	#Testing to see how program reacts to two compounds with same mass
	'''
	fuck_shit_up = {'Other Biofluids': u'SHIT', 'Feces': 'crap', 'Urine': 'poo', 'Name': u'Trying to fuck it up'}
	#print '\n fuck it up \n', fuck_shit_up
	fuck_up_df = pandas.DataFrame(fuck_shit_up, index=[74.08439833])
	#print fuck_up_df

	shitty_biofluids = pandas.concat([biofluids, fuck_up_df], axis=0)
	#print shitty_biofluids
	'''
	#Iterate through biofluids rows, extract mass. Iterate through 
	#metlin info. if they're equal, add metlin stuff to biofluids
	combined_df = pandas.DataFrame()

	#row_id is the mass, row is the urine, feces, etc info.
	for row_id, row in biofluids.iterrows():
		#mass is the mass, met_row is the row info - isomers, adduct, etc
		for mass, met_row in curated_df.iterrows():
			#if the hmdb mass matches a mass from metlin
			#(Which it should have to) Then add the metlin isomer
			#and adduct information to the biofluid information
			if row_id == mass:
				concat_data = pandas.concat([row, met_row])
				combined_df = pandas.concat([combined_df, concat_data], axis=1)
	
	nearly_final_output = combined_df.T.set_index('HMDB_ID')
	#Get the list of columns
	column_order_original = nearly_final_output.columns.tolist()
	column_order_final = ['mass', 'dppm', 'adduct', 'Isomers', 'Name',
				    'Urine', 'Feces', 'Other Biofluids']

	final_output = nearly_final_output[column_order_final]
	return final_output


def add_isomers_to_data(metlin_data_curated, isomers_data):
	'''Takes a main dataframe with an empty isomers, a dataframe 
	containing hte number of isomers at each mass, and adds the 
	isomer information to the main dataframe'''
	for mass_value in metlin_data_curated.index:
		#Add number of isomers at a given mass to each row
		#representing that mass
		for metlin_mass in isomers.index:
			if mass_value == metlin_mass:
				metlin_data_curated.loc[mass_value, 'Isomers'] = isomers[metlin_mass]
	return metlin_data_curated

def nested_dict_to_dataframe(nested_dict):
	#print nested_dict
	#iterate over masses:
	output = pandas.DataFrame()
	for mass in nested_dict:
		#print mass
		#iterate over multiple(?) hmdb ids per mass
		for hmdb_id in nested_dict[mass]:
			#print nested_dict[mass][hmdb_id]
			#Create dataframe and add columns to it
			data = pandas.DataFrame(nested_dict[mass][hmdb_id], index=[mass])
			data['HMDB_ID'] = hmdb_id
			data['mass'] = mass
			#print data
			output = pandas.concat([output, data])
			#print output
	return output




'''
hmdb_identified = search_HMDB(toy_masses, metabolites) 
print 'Found HMDB info: ', hmdb_identified

biofluids_nested_dict = get_hmdb_info(hmdb_identified)
print '\n Biofluid_info as nested dict',biofluids_nested_dict
biofluid_dataframe = nested_dict_to_dataframe(biofluids_nested_dict)
print '\nBiofluid dataframe\n', biofluid_dataframe
print biofluid_dataframe

final_output = combine_info(biofluid_dataframe, metlin_data, isomers)
print '\nfinal data to output to file\n', final_output

#output_file = open(output_file)
final_output.to_csv(output_filename)
'''


