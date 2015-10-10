Tested in Python v2.7.6
- Required Python libraries are listed in 'requirements.txt'
- You may also need to install libxml2-dev, libxslt-dev,
    and python-dev on your machine.

Before searching HMDB, you must modify the HMDB database XML file. 
You have to remove the repetitve <?xml> headers and add a new root tag. 

Do this by running "python remove_excess_xml_declarations.py"
  - Make sure to update the "hmdb_all_metabolites" and 
    "output_file" variables to the paths on your computer.

To search HMDB:
  0 - Make sure you HMDB xml file has been reformatted (see above)
  1 - Download the results of METLIN batch search, in csv format.
  2 - Update the paths in parse_metlin_csv.py to your HMDB database, 
      METLIN csv file, and output file.
  3 - Run "python parse_metlin_csv.py"
