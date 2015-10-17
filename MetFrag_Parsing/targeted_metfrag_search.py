import pandas as pd 
import scipy.io as sio


#Real data
#mass_spec_hmdb = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/sept28_data/hmdb_mass_spec_merged.csv'
ms2_data_matlab = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/ms2_data/neg_10p_MS2_aligned.2015-09-23.mat'
metfrag = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/MetFrag_Parsing/MetFrag2.2-CL.jar'
parameter_file = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/MetFrag_Parsing/example_parameter_file.txt'
output_mz_intensity_path = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/MetFrag_Parsing/mz_intensities/'
#Toy Data
mass_spec_hmdb = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/MetFrag_Parsing/toy_mass_spec_hmdb.csv'


#0.) Go through mass_spec/hmdb info file with ID, chemical 
#     database ID's. Pull out ID'd compounds with MS2 values 
def get_mz_intensities(mass_spec_ID, ms2_data):
	'''INPUT - the sample number for mass mass_spec
	OUTPUT - Two numpy data arrays, mz values and intensities
	'''
	#TODO Do we want to take values from average mz and intensity 
	#columns or from the ms2_mz and MS2_intensities columsn???
	mz_values = ms2_data['MS2_data']['averageMS2'][0][0][mass_spec_ID]['MS2_mz'][0]
	intensities = ms2_data['MS2_data']['averageMS2'][0][0][mass_spec_ID]['MS2_intensity'][0]
	return mz_values, intensities

def write_mz_intensities_to_file(mz, intensities, output_path):
	with open(output_path,'w') as f:
		for i in range(0,len(intensities)):
			f.write('%s\t%s\n' % (mz[i][0], intensities[i][0]))

mass_spec_df = pd.read_csv(mass_spec_hmdb, header=0)
#print mass_spec_df

#load matlab with ms2 data
ms2_data = sio.loadmat(ms2_data_matlab)

#Go through each identified hmdb compound. If it has an MS2, 
#Run that through MetFrag and ge the results
for i in range(0, len(mass_spec_df)):
	if mass_spec_df.ix[i]['MS2'] == 1:
		#minus 1 because python indexing starts at 0
		mass_spec_ID = mass_spec_df.ix[i]['ID'] - 1
		print 'ID: %s' % mass_spec_ID
		#Get mz values and intensities
		mz_values, intensities = get_mz_intensities(mass_spec_ID, ms2_data) #TODO
		mz_intensity_path = output_mz_intensity_path + 'mz_intensity_%s.txt' % mass_spec_ID
		#TODO write mz_values and intensities to a txt file
		write_mz_intensities_to_file(mz_values, intensities, mz_intensity_path)


#1.) Get the m/z and intensities from each ID
#1.5) Write m/z values to file
#1.75) Sum the intensities to compare how well MetFrag matches
#2.) Modifty the metfrag parameter file to include your search criteria
#     and new data.txt file (m/z and intensities file)

#3.) Evaluate the output of Metfrag, scoring % peaks explained
#     and % of Intensities explained