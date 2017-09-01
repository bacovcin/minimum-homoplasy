# Maximum Synapomorphy
A python package for implementing a maximum synapomorphy method for finding
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
