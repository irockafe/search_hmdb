import csv
from sys import argv
import pandas
from lxml import etree
from collections import defaultdict
#You might also requires libxml-dev and libxslt-dev packages!

'''
INPUT    - A csv file downloaded from a METLIN batch search. 
		   concatenate multiple METLIN files into one if you 
		   are searching more than 500 compounds (METLIN batch 
		   search limits you to 500 compounds per query)

FUNCTION - This script will search a downloaded HMDB database
		   (xml-format) for the masses identified by METLIN and add 
		   HMDB information, such as if the metabolite is from urine or 
		   feces, to a new csv file

OUTPUT - A csv file containing columns for:
		['inputmass', 'mass', 'dppm', 'adduct', 'Isomers', 
		'Name','Urine', 'Feces', 'Other Biofluids'] 

NOTE - If METLIN did not find any hits for a mass you searched,
it will not save that mass to the metlin-csv, and this code will not 
search HMDB for that mass.
'''
#Path to csv file
metlin_csv = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/sept28_data/negativePeaks_sept28_METLIN.csv'
#Path to HMDB xml database
hmdb_all_metabolites = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/hmdb_metabolites/removed_xml_hmdb_metabolites.xml'
#Path to output file
output_filename = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/sept28_data/negative_ion_HMDB_hits.csv'
mass_spec_csv = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/sept28_data/prelim_negativePeaks_Sept28_instrument.csv'
merged_csv = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/sept28_data/hmdb_mass_spec_merged.csv'

'''
#For debugging/testing purposes. Tests not included yet in git
metlin_csv = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/metlin_toy_data.csv'
#metlin_csv = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/smaller_toy_data.csv'
hmdb_all_metabolites = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/toy_hmdb/removed_xml_toy_hmdb.xml'
toy_hmdb_directory = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/toy_hmdb/'
output_filename = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/sept28_data/negative_ion_HMDB_hits_formula.csv'
mass_spec_csv = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/sept28_data/prelim_negativePeaks_Sept28_instrument.csv'
merged_csv = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/sept28_data/hmdb_mass_spec_merged.csv'
'''



def get_data_from_hmdb(xml_file, masses, mass_deviation_tolerance):
	'''
	INPUT - path to an xml database file, a list of masses you
			want to search, and the maximum deviation tolerance 
			for considering	floating-point integers to be equal
	FUNCTION - Find entries in hmdb that match our masses,
				get relevant info - name, biofluid location

	OUTPUT - A nested dictionary:
			{Mass1: {HMDB_ID-1:{Urine: x, Feces: x, 
								Other Biofluids: 'string'},
					 HMDB_ID-2:..},
			 Mass2: {..},...}
	'''

	log_file = open('log_parse_hmdb.txt', 'w')
	output_nested_dict = defaultdict(dict)
	#Go through each metabolite, one at a time
	xml_tree = etree.iterparse(xml_file, tag='metabolite')
	count = 0
	for event, element in xml_tree:
		count += 1
		if count % 5000 == 0:
			print 'Working on compound #%s' % count
		#Get the MW
		monisotopic_tag_text = element.findtext('monisotopic_moleculate_weight')
		#If no MW listed, skip over that metabolite and write to log 
		#file that you did so.
		if monisotopic_tag_text == '':
			compound_name = element.findtext('name')
			hmdb_id = element.findtext('accession') 
			log_file.write('Skipped over %s, %s, because no MW was listed\n' % (hmdb_id, compound_name))
			#Go to next iteration of for loop
			continue

		MW = float(element.findtext('monisotopic_moleculate_weight'))
		hmdb_id = element.findtext('accession')

		for mass in masses:
			#floating-point math is annoying
			if abs(MW - mass) < mass_deviation_tolerance:
				hmdb_id = element.findtext('accession')
				compound_name = element.findtext('name')

				metabolite_info = {'Chemical Formula': element.findtext('chemical_formula'),
						'Chemspider ID' : element.findtext('chemspider_id'),
						'Pubchem ID': element.findtext('pubchem_compound_id')
						}
				#print metabolite_info
				#print 'Found Mass %s, hmdb: %s, name: %s'% (MW, hmdb_id, compound_name)

				#Get a list of biofluid locations
				biofluid_locations =  element.find('biofluid_locations').findall('biofluid')
				biofluid_dict = get_biofluid_locations(biofluid_locations, compound_name)
				#Add the biofluid info to the chemical formula, etc. info
				hmdb_metabolite_info = metabolite_info.copy()
				hmdb_metabolite_info.update(biofluid_dict)


				#Add new HMDB entries of the same MW to the same 
				#dictionary key
				output_nested_dict[MW].update({hmdb_id: hmdb_metabolite_info})
		#clean up memory once you're done with each metabolite 		
		element.clear()
	log_file.close()
	#Raise error if no hmdb matches were found
	if not output_nested_dict:
		raise ValueError('No hits were found in HMDB for the masses\
						 you searched: %s' % masses)
	return output_nested_dict
	

def get_biofluid_locations(biofluid_xml, compound_name):
	'''
	INPUT - xml tree of the <biofluid_locations> tag from hmdb 
	FUNCTION - parse that xml and get the text within each child tag 
			   to see which biofluids it is present in (This should 
			   	only be <biofluid> tags. but could change in 
			   	future database releases)
	OUTPUT - {Urine: (1/0), Feces: (1/0), Other Biofluids: ['string']}
	'''
	biofluid_info = {'Name': compound_name, 'Urine': 0, 'Feces': 0, 'Other Biofluids': []}
	other_fluids = []
	for fluid in biofluid_xml:
		if fluid.text == 'Urine':
			biofluid_info['Urine'] = 1
		elif fluid.text == 'Feces':
			biofluid_info['Feces'] = 1
		else:
			other_fluids.append(fluid.text)
	#Add the non-urine/feces biofluids
	biofluid_info['Other Biofluids'] = ', '.join(other_fluids)

	return biofluid_info

def nested_dict_to_dataframe(nested_dict):
	'''
	INPUT - A nested dictionary containing masses, HMDB_IDs, and 
			biofluid info:
			{Mass1: {HMDB_ID-1:{Urine: x, Feces: x, 
								Other Biofluids: 'string'},
					 HMDB_ID-2:..},
			 Mass2: ..}
	FUNCTION - Convert a nested dict into a dataframe
	OUTPUT - A Pandas dataframe with rows labeled by mass,
			 and columns labeled by HMDB_ID, Urine, Feces, Other Biofluids
	'''
	#Convert a nested dictionary into dataframe
	output = pandas.DataFrame()
	for mass in nested_dict:
		#print mass
		#iterate over multiple(?) hmdb ids per mass
		for hmdb_id in nested_dict[mass]:
			#print nested_dict[mass][hmdb_id]
			#Create dataframe of the Urine, Feces, other biofluids info,
			#Then add the HMDB_ID and mass to that.
			data = pandas.DataFrame(nested_dict[mass][hmdb_id], index=[mass])
			data['HMDB_ID'] = hmdb_id
			data['mass'] = mass
			output = pandas.concat([output, data])
	return output

def combine_hmdb_metlin(hmdb_df, metlin_df, MW_deviation_tolerance):
	'''
	INPUT - A pandas dataframe of hmdb information, metlin information, and 
			the maximum difference that two floating-point integers can 
			have and still be considered the same number.
	FUNCTION - Merge the information found in HMDB with the inputmass,
			   adduct, and dppm values from METLIN
	OUTPUT - A pandas dataframe containing rows labeled by HMDB_ID
			 and columns labeled by ['inputmass', 'mass', 'dppm', 
			 'adduct', 'Isomers', 'Name', 'Urine', 'Feces', 
			 'Other Biofluids']
	'''
	#TODO - pre-allocate memory for this dataframe.
	combined_df = pandas.DataFrame()
	#Combine the metlin and hmdb data into one dataframe
	for hmdb_mass, hmdb_row in hmdb_df.iterrows():
		for metlin_mass, metlin_row in metlin_df.iterrows():
			#if the masses in hmdb match a mass from metlin, which 
			#they have to, add the metlin info to hmdb_info
			if abs(hmdb_mass - metlin_mass) < MW_deviation_tolerance:
				concat_row = pandas.concat([metlin_row, hmdb_row])
				combined_df = pandas.concat([combined_df, concat_row], axis=1)

	#transpose dataframe and set HMDB_ID as index
	combined_df = combined_df.T
	combined_df = combined_df.set_index('HMDB_ID')
	#Re-arrange the columns
	column_order_final = ['inputmass', 'mass', 'dppm', 'adduct', 'Isomers', 'Name',
				    'Chemical Formula', 'Chemspider ID', 'Pubchem ID','Urine', 'Feces', 'Other Biofluids']
	final_df = combined_df[column_order_final]
	return final_df

def merge_hmdb_instrument_info(mass_spec_csv, hmdb_hits_csv, merged_csv):
	'''
	INPUT - list of file paths
	FUNCTION - Merges information from the mass spec and info from
				hmdb into one document on their m/z
	OUTPUT - csv file
	'''
	mass_spec_df = pandas.read_csv(mass_spec_csv)
	hmdb_df = pandas.read_csv(hmdb_hits_csv)
	#rename a header from inputmass to m/z
	hmdb_df = hmdb_df.rename(columns={'inputmass':'m/z'})
	hmdb_df = hmdb_df.rename(columns={'mass':'Metlin Mass'})
	#merge dataframes based on m/z value
	merged_df = pandas.merge(hmdb_df, mass_spec_df, how='outer', on=['m/z'])
	#reorder the output
	column_order = [u'RT', u'm/z', u'Metlin Mass', u'adduct', u'HMDB_ID', u'MS2? (1= YES)', 
					u'dppm', u'Isomers', u'Name', u'Chemical Formula', 'Urine', 
					'Feces', u'Chemspider ID', u'Pubchem ID',
					'Other Biofluids', u'std UW', u'Proc blank PPL', u'Proc blank', 
					u'1a', u'4p', u'1p', u'Inst blank', u'1p PPL', u'4p PPL', u'1a PPL']
	merged_df = merged_df[column_order]
	merged_df.to_csv(merged_csv, encoding='utf-8')

#This is the maximum difference between the MW we are 
#searching for and the MW we find in HMDB that we will consider a match 
MW_deviation_tolerance = 1e-9

#import METLIN-csv to pandas, with first row as the header
metlin_data = pandas.read_csv(metlin_csv, header=0)

#Count number of repeat input mass values - the number of isomers
isomers = metlin_data['mass'].value_counts()
masses = isomers.index
print '\nSearching HMDB for %s distinct masses' % masses.size

#Choose which METLIN data you will add to the HMDB data
selected_metlin_data = pandas.concat([metlin_data['inputmass'],metlin_data['mass'],metlin_data['adduct'],
									  metlin_data['dppm']],
									  axis=1)
#print '\nselected data:\n %s' % selected_metlin_data

#Remove duplicate masses and set row indices to mass values
selected_metlin_data = selected_metlin_data.drop_duplicates().set_index('mass')
selected_metlin_data['Isomers'] = isomers
#print '\nSelected Data + Isomers:\n%s'% selected_metlin_data

#Find hmdb matches
hmdb_matches = get_data_from_hmdb(hmdb_all_metabolites, 
								 masses, MW_deviation_tolerance)
#print '\n I found the following compounds in HMDB: \n %s' % nested_dict

#Convert hmdb_matches to a pandas dataframe
hmdb_df = nested_dict_to_dataframe(hmdb_matches)
#print hmdb_df
#print '\nHMDB_datafame: \n %s' % hmdb_df 

#Combine METLIN data with HMDB data
combined_df = combine_hmdb_metlin(hmdb_df, selected_metlin_data, MW_deviation_tolerance)
print '\nI found %s compounds in HMDB\n' % combined_df.shape[0]
print 'I might have skipped some HMDB entries. See the log file: log_parse_hmdb.txt\n'
print 'Saved the output at %s' % output_filename

#Save METLIN/HMDB data to file
combined_df.to_csv(output_filename, encoding='utf-8')

merge_hmdb_instrument_info(mass_spec_csv, output_filename, merged_csv)
print '\nMerged Mass-spec instrument data with hmdb data at %s' % merged_csv
