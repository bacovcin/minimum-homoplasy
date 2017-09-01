# Maximum Synapomorphy
A python package for implementing a maximum synapomorphy method for finding
optimal phylogenetic trees. 

# Dependencies
This code requires Python3 and has been verified to work with Python 3.6.2.

# Usage
python fit_data.py CSV_FILE_WITH_DATA

# Input Format
Rows for characters, columns for languages

Two optional columns:
InheritedValue - Gives the state for the proto-language (when independently recoverable)
Weighting - Gives the weight for states in this character

Special State Value:
NA - A state named NA is ignored by the software (useful for dealing with polymorphic characters)
