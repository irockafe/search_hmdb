import csv
from sys import argv
import pandas
from lxml import etree
from collections import defaultdict
#You might also requires libxml-dev and libxslt-dev packages!

'''
This code takes the csv file downloaded from at 
METLIN batch search. It will search a downloaded HMDB database
(xml-format) for the masses identified by METLIN and add 
HMDB information, like if the metabolite is from urine or 
feces, to a new csv file

Metlin does not output masses where nothing was found.
'''
#You have to point python to your metlin csv in the
#command line
#script, metlin_csv = argv

#Path to csv file
metlin_csv = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/negativePeaks_sept28_METLIN.csv'
#metlin_csv = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/metlin_toy_data.csv'
#metlin_csv = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/smaller_toy_data.csv'

hmdb_all_metabolites = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/hmdb_metabolites/removed_xml_hmdb_metabolites.xml'
toy_hmdb_all_metabolites = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/toy_hmdb/removed_xml_toy_hmdb.xml'

hmdb_directory = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/hmdb_metabolites/'
toy_hmdb_directory = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/toy_hmdb/'

output_filename = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/toy_hmdb/all_mass_tests.csv'


#metabolites = BeautifulSoup(open(hmdb_all_metabolites), "lxml").find_all('metabolite')

#save memory by not opening two giant files
#hmdb_soup = BeautifulSoup(open(hmdb_all_metabolites), "lxml")
#metabolites = hmdb_soup.find_all('metabolite')

def get_data_from_hmdb(xml_file, masses, mass_deviation_tolerance):
	'''
	INPUT - path to an xml database file and a list of masses you
			want to search
	FUNCTION - Find entries in hmdb that match our masses,
				get relevant info - name, biofluid location

	OUTPUT - {Mass1: {HMDB_ID-1:{Urine: x, Feces: 0, Others: 'string'}}}
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
		#print '%s -- %s\n\n' % (event, element.findtext('monisotopic_moleculate_weight'), )
		#Get the MW
		monisotopic_tag_text = element.findtext('monisotopic_moleculate_weight')
		'''
		print monisotopic_tag_text
		print monisotopic_tag_text == ''
		'''
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
			'''
			if hmdb_id == 'HMDB11635':
				print '\n',hmdb_id
				print mass
				print type(masses[0])
				print MW
				print type(MW)
				print abs(MW-mass) < mass_deviation_tolerance, '\n'
			'''
			if abs(MW - mass) < mass_deviation_tolerance:

				hmdb_id = element.findtext('accession')
				
				compound_name = element.findtext('name')
				#print 'Found Mass %s, hmdb: %s, name: %s'% (MW, hmdb_id, compound_name)

				#Get a list of biofluid locations
				biofluid_locations =  element.find('biofluid_locations').findall('biofluid')
				
				biofluid_dict = get_biofluid_locations(biofluid_locations, compound_name)

				#Add new entries of the same MW to that MW's dictionary
				#entry
				output_nested_dict[MW].update({hmdb_id: biofluid_dict})

				#hmdb_properties = {hmdb_id: biofluid_info}

			#output_nested_dict[MW] = {hmdb_id:biofluid_info}

		element.clear()
	log_file.close()
	#Raise error if no hmdb matches were found
	if not output_nested_dict:
		raise ValueError('No hits were found in HMDB for the masses\
						 you searched: %s' % masses)
	return output_nested_dict
			#output_nested_dict[MW] = 
			#Get 			

def get_biofluid_locations(biofluid_xml, compound_name):
	'''
	INPUT - xml of the <biofluid_locations> tag from hmdb 
	FUNCTION - parse that xml and return the biofluids biofluids present 
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
	#append the non-urine/feces biofluids
	biofluid_info['Other Biofluids'] = ', '.join(other_fluids)

	return biofluid_info

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

def combine_hmdb_metlin(hmdb_df, metlin_df, MW_deviation_tolerance):
	'''
	duplicate_mass =pandas.DataFrame({'Feces': 'shit', 
		'Urine':'piss', 'Other Biofluids': "you don't wanna know", 
		'HMDB_ID': 'asdf', 'mass': 74.084398}, index=[74.08439833])
	
	duplicate_mass_test = pandas.concat([hmdb_df,  duplicate_mass], axis=0)
	print duplicate_mass_test
	'''
	
	combined_df = pandas.DataFrame()
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
	column_order_final = ['inputmass', 'mass', 'dppm', 'adduct', 'Isomers', 'Name',
				    'Urine', 'Feces', 'Other Biofluids']
	final_df = combined_df[column_order_final]


	#output_df = pandas.concat([metlin_df, duplicate_mass_test],axis=1)
	return final_df

#This is the maximum difference between a the MW we are 
#searching for and the MW we find in HMDB. 
#We have to do this because we cannot compare 
#floating-point integers via equality - 29.0 could 
#be represented as 28.999999999 and 29.000000000
MW_deviation_tolerance = 1e-9
#import csv to pandas, with first row as the header
metlin_data = pandas.read_csv(metlin_csv, header=0)
#metlin_data = metlin_data.set_index('mass')
#print '\n Metlin Data: \n %s' % metlin_data
#Count number of repeat input mass values - the number of isomers
isomers = metlin_data['mass'].value_counts()
#print '\n Isomers \n%s' % isomers
masses = isomers.index
print '\nSearching HMDB for %s distinct masses' % masses.size
#print '\n Masses \n%s' % masses

selected_metlin_data = pandas.concat([metlin_data['inputmass'],metlin_data['mass'],metlin_data['adduct'],
									  metlin_data['dppm']],
									  axis=1)
#print '\nselected data:\n %s' % selected_metlin_data

#Remove duplicate masses and set row indices to mass values
selected_metlin_data = selected_metlin_data.drop_duplicates().set_index('mass')
selected_metlin_data['Isomers'] = isomers
#print '\nSelected Data + Isomers:\n%s'% selected_metlin_data

#pandas.set_option('precision', 10)
nested_dict = get_data_from_hmdb(hmdb_all_metabolites, 
								 masses, MW_deviation_tolerance)
#print '\n I found the following compounds in HMDB: \n %s' % nested_dict

hmdb_df = nested_dict_to_dataframe(nested_dict)
#print '\nHMDB_datafame: \n %s' % hmdb_df 


combined_df = combine_hmdb_metlin(hmdb_df, selected_metlin_data, MW_deviation_tolerance)
print '\nI found %s compounds in HMDB\n' % combined_df.shape[0]
print 'Some metabolites may have been skipped. See the log file: log_parse_hmdb.txt\n'
print 'Saved the output at %s' % output_filename


combined_df.to_csv(output_filename, encoding='utf-8')
