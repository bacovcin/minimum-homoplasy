# Minimum Homoplasy
A python package for implementing a minimum homoplasy method for finding
optimal phylogenetic trees. 

# Dependencies
This code requires Python3 and has been verified to work with Python 3.6.2.
It also requires:
	[BioPython](http://biopython.org/) for pretty tree printing: pip install biopython
	[PAUP*](http://paup.phylosolutions.com/): Used version 4a157


# Usage
python fit_data.py CSV_FILE_WITH_DATA

# Input Format
Columns for characters, rows for terminals (terminal name in first row)

Two optional rows:
InheritedValue - Gives the state for the proto-language (when independently recoverable)
Weighting - Gives the weight for states in this character

Special State Value:
NA - A state named NA is ignored by the software (useful for dealing with polymorphic characters)

# Using Indo-European Data
One of our examples uses Indo-European data from the project (Computational Phylogenetics In Historical Linguistics)[https://www.cs.utexas.edu/users/tandy/histling.html]

The data can be downloaded from this link: (http://www.cs.rice.edu/~nakhleh/CPHL/IEDATA_112603)[http://www.cs.rice.edu/~nakhleh/CPHL/IEDATA_112603]

The data must then be saved in the home directory of the repository as IEDATA.

The python script convert\_IE\_data.py can be used to convert the data to the format used by the Minimum Homoplasy algorithm.

## Conversion Flags (takes the form flag=value)
P=Weight for phonological characters
M=Weight for morphological characters
L=Weight for lexical characters
exinh=Exceptional ancestral values (Example ex=M1=1 would set the ancestral value for M1 to 1), note that phonological characters are assumed to have an ancestral value of 1
exc=Exclusion of character (Example exc=M11 would exclude M11 from the output dataset)
specw=Special weight for specific character (Example specw=M10=1000 would set the weight for M10 to 1000 while all other morphological characters would get the M weight)

# PAUP-Based Weighted Maximum Parsimony and Weighted Maxmium Compatibility
We have implemented a Python script (run\_paup.py) that takes as a command line argument a csv file in the Maximum Synapomorphy format and
performs the following operations:
1) Converts the data into a Nexus file for use in PAUP (state values are converted to A-Z to avoid encoding issues)
2) Updates the PAUP templates in the PAUP template folder to use the data from the new Nexus file
3) Runs the PAUP templates (and uses some intermediate output files to calculate weighted compatibility)
4) Outputs the majority rule consensus tree for maximum parsimony and weighted maximum compatibility (midpoint rooting is used in the output tree, but the methods actually only optimize an unrooted tree)

# Makefile
The exact commands used to generate data for the paper "A New Method for Computational Cladistics: An Afro-Asiatic Case Study" can be found in the Makefile. To replicate our analysis: 
1) Download the Indo-European data as described above
2) Run: make casestudy
