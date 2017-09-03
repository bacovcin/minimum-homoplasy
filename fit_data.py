from itertools import *
from time import perf_counter as perfc
import sys
import copy
from random import randint, shuffle, random, sample
from Bio.Phylo import BaseTree
from Bio import Phylo

def prettyprint_tree(tree, file=None):
    # Convert the "tree" object (list of clades) to a BioPython tree 
    # to take advantage of their output methods
    def create_ntree(tree):
        ntree = BaseTree.Clade()
        for key in tree:
            el = tree[key]
            if len(el.values()) > 0:
                ntree.clades.append(create_ntree(el))
            else:
                ntree.clades.append(BaseTree.Clade(name=list(key)[0]))
        return ntree

    # Recursively add a node to a dictionary representation of a tree 
    # (each node is a key containing a dictionary with subnodes)
    def create_tree_dict(curdict,curnode):
        for key in curdict:
            if curnode.issubset(key):
                curdict[key] = create_tree_dict(curdict[key],curnode)
                break
        else:
            curdict[curnode] = {}
        return curdict

    # Sort the clades from largest to smallest
    new_tree = sorted(tree, key=lambda x:-len(x))
    # Build a dictionary representation of the tree
    tree_dict = {}
    for clade in new_tree:
        tree_dict = create_tree_dict(tree_dict,clade)
    # Convert the dictionary representation to a BioPython Tree object
    ntree = BaseTree.Tree(create_ntree(tree_dict))
    # Use the BioPython print method
    Phylo.draw_ascii(ntree,file=file)
    try:
        Phylo.draw(ntree)
    except:
        pass
    return

def readin_chars(arg):
    # Creates data dictionary object (with four keys chars, inhchar, weights,
    # and expected_failures)
    infile = open(arg)
    # Get a list of character names
    names = infile.readline().rstrip().split(',')
    # Initialize the data dictionary and language list
    data = {}
    lgg = []
    # Initialize characters, known inherited features, and weights
    data['chars'] = {} # Characters
    data['inhchar'] = {} # Known inherited characters
    data['weights'] = {} # Character weights
    for name in names:
        if name != '':
            data['chars'][name] = {}
            data['inhchar'][name] = ''
            data['weights'][name] = 1
    data['expected_failures'] = 0 # Sum of character weights
    for line in infile:
        # Read in the row
        s = line.rstrip().split(',')
        if s[0] == 'Weighting':
            for i in range(len(s[1:])):
                data['weights'][names[i+1]] = int(s[i+1])
        elif s[0] == 'InheritedValue':
            for i in range(len(s[1:])):
                data['inhchar'][names[i+1]] = s[i+1]
        elif s[0] not in ['Weighting','InheritedValue']:
            # Add the language name to the language list
            lgg.append(s[0])
            for i in range(len(s[1:])):
                try:
                    data['chars'][names[i+1]][s[i+1]] += [s[0]]
                except KeyError:
                    data['chars'][names[i+1]][s[i+1]] = [s[0]]
                if s[i+1] != data['inhchar'][names[i+1]]:
                    data['expected_failures'] += data['weights'][names[i+1]]
    data['lgg'] = lgg
    return data

def find_possible_clades(data):
    # Create a dictionary of possible clades with their supporting data/score
    posclades = {}
    # Initialise clades for all terminal nodes
    for lg in data['lgg']:
        curclade = frozenset([lg])
        posclades[curclade] = {}
        posclades[curclade]['evid'] = ['Terminal']
        posclades[curclade]['score'] = 1
    # Examine each character
    for char in data['chars']:
        # Examine each value
        for val in data['chars'][char]:
            # Do not include if it is an inherited character
            if val == data['inhchar'][char] or val in ('NA','N/A'):
                continue
            curclade = frozenset(data['chars'][char][val])
            try:
                posclades[curclade]['evid'].append(char)
                posclades[curclade]['score'] += data['weights'][char]
            except KeyError:
                posclades[curclade] = {}
                posclades[curclade]['evid'] = [char]
                posclades[curclade]['score'] = data['weights'][char]
    return posclades 

def find_conflicts(clades):
    # Return a dictionary that maps all possible clades to other possible
    # clades that they are incompatible with
    def testTreeSuit(c1,c2):
        # Test is two clades are compatible
        if c1.intersection(c2) in [c1,c2,set([])]:
            return True
        else:
            return False
    # Initialise dictionary
    conflicts = {}
    # Go through all clade pairs testing for conflicts
    for clade1 in clades:
        conflicts[clade1] = []
        for clade2 in clades:
            if not testTreeSuit(clade1, clade2):
                conflicts[clade1].append(clade2)
    return conflicts 

def build_tree_from_list(clist,conflicts,posclades):
    # Extract a tree (list of clades) from
    # an ordered list of (possibly conflicting) clades
    clist = copy.copy(clist)

    # Initialize list and score
    tree = []
    score = 0
    # Go through list
    while len(clist) > 0:
        # Remove the next clade on the list
        clade = clist.pop(0)
        # Add the clade to the tree
        tree.append(clade)
        # Remove from the list of clades any clades incompatible with the clade
        # just added to the tree and increment the score with the weights
        # relevant to the excluded clades
        for incomp_clade in conflicts[clade]:
            try:
                clist.remove(incomp_clade)
                score += posclades[incomp_clade]['score']
            except ValueError:
                continue
    return (tree,score)

def sort_clades(clades, conflicts, posclades):
    # Rank each clade by (1) strength of incompatible clades and (2) number of
    # memebers
    scores = [(x,sum([posclades[y]['score'] for y in conflicts[x]])/len(x))
              for x in clades]
    return [x[0] for x in sorted(scores,key=lambda x:x[1])]

def branch_and_bind(sorted_clades, posclades, conflicts, 
                    optima, curscore, curtree):
    # Use the branch_and_bind method to find the optimal trees for a given
    # dataset.
    # Recursively move through the list of clades
    # For each clade:
    #  Try keeping the clade (removing all remaining incompatible with penalty) 
    #  Try not keeping the clade (with penalty for removed clade)

    # Extract the current best from current optimum list
    curbest = optima[0][0]

    # Get the current worked on clade
    clade = sorted_clades.pop(0)

    # Print progress
    print('Current Best = ' + str(curbest) + ': ' + str(len(sorted_clades)))

    # Create a copy of the sorted clades for removal of clades incompatible
    # with current clade
    keeplist = copy.copy(sorted_clades)

    # Initialise the current score
    keep_score = curscore

    # Add the clade to the keep tree
    keep_tree = curtree + [clade]

    # Remove all incompatible clades
    for incomp_clade in conflicts[clade]:
        try:
            keeplist.remove(incomp_clade)
            keep_score += posclades[incomp_clade]['score']
        except:
            continue

    # Check and see if the search needs to be bounded
    if keep_score <= curbest:
        # If there are more nodes to search, continue the search
        if len(keeplist) > 0:
            new_optima = branch_and_bind(keeplist, posclades, conflicts,
                                         optima, keep_score, keep_tree)
            if new_optima[0][0] == curbest:
                for x in new_optima:
                    if x not in optima:
                        optima.append(x)
            elif new_optima[0][0] < curbest:
                optima = new_optima
        # Since this is the end of the search, update the optima
        else:
            if (keep_score == curbest and 
                True not in [keep_tree == x[1] for x in optima]):
                optima.append((keep_score, keep_tree))
            else:
                optima = [(keep_score, keep_tree)]

    # Update the current best score for possible new optima
    curbest = optima[0][0]

    # Calculate the score if the current clade has been removed
    drop_score = curscore + posclades[clade]['score']

    # Check and see if the search needs to be bounded
    if drop_score <= curbest and conflicts[clade] != []:
        # If there are more nodes to search, continue the search
        if len(sorted_clades) > 0:
            new_optima = branch_and_bind(sorted_clades, posclades, conflicts,
                                         optima, drop_score, curtree)
            if new_optima[0][0] == curbest:
                for x in new_optima:
                    if x not in optima:
                        optima.append(x)
            elif new_optima[0][0] < curbest:
                optima = new_optima
        # Since this is the end of the search, update the optima
        else:
            if (drop_score == curbest and 
                True not in [curtree == x[1] for x in optima]):
                optima.append((drop_score,curtree))
            else:
                optima = [(drop_score,curtree)]
    return optima

def write_output(optima, posclades, arg, maxscore, elapsed_time):
    # Get the output filename/directory from input filename/directory
    s = '/'.join(arg.split('.')[0].split('/')[1:])
    outfilename = 'outputs/' + s + '.txt'
    outfile = open(outfilename,'w')
    # Write out headers (time and # of optima)
    outfile.write('Elapsed Time: ' + str(elapsed_time)+' seconds\n')
    outfile.write('# of trees: ' + str(len(optima)) + '\n')
    # Print out all the trees, their associated scores, what percentage of the
    # weighted score was excluded and a list included/excluded clades and the
    # evidence for each clade
    i = 1
    for optimum in optima:
        print('Tree #'+str(i)+':')
        outfile.write('Tree #' + str(i) +':\n')
        i += 1
        prettyprint_tree(optimum[1])
        prettyprint_tree(optimum[1], outfile)
        outfile.write('\tScore: ' + str(optimum[0]))
        outfile.write('\tExclusion rate: ' +
                      str(round(float(optimum[0])/float(maxscore),5))+'\n')
        outfile.write('\tIncluded Groups:\n')
        for clade in sorted(optimum[1], key =lambda x:len(x)):
            evid = posclades[clade]['evid']
            outfile.write('\t\t'+', '.join(clade)+':\n')
            outfile.write('\t\t\tSupported by ' + str(len(evid))+
                          ' Characters: '+ ', '.join(evid)+'\n')
        outfile.write('\tExcluded Groups:\n')
        for clade in sorted(posclades.keys(), key=lambda x:len(x)):
            if clade not in optimum[1]:
                evid = posclades[clade]['evid']
                outfile.write('\t\t'+', '.join(clade)+':\n')
                outfile.write('\t\t\tSupported by ' + str(len(evid))+
                              ' Characters: ' + ', '.join(evid)+'\n')
    outfile.close()
    return

if __name__ == "__main__":
    debug = False
    # Go through commandline arguments
    for arg in sys.argv:
        # Set the debug flag
        if arg == 'debug':
            debug = True
        # Only csv files are valid datasets
        if arg[-3:] != 'csv':
            continue
        # Save the start time to get accurate run times
        start_time = perfc()
        # Parse the dataset
        data = readin_chars(arg)
        # Extract the list of possible clades
        posclades = find_possible_clades(data)
        # Print out a list of clades and their evidence (for debugging)
        if debug:
            for clade in sorted(posclades.keys(),key = lambda y:len(y)):
                print('Clade: ' + str(clade))
                print('\tEvidence: ' + str(posclades[clade]['evid']))
        # Calculate which clades conflict with each other
        conflicts = find_conflicts(posclades.keys())
        # Print out the list of conflicting clades (for debugging)
        if debug:
            for clade in conflicts:
                print('Clade = ' + str(clade))
                print('\tConflicts: ' + str(conflicts[clade]))
        # Sort the clades for optimal algorithm completion
        sorted_clades = sort_clades(posclades.keys(), conflicts, posclades)
        # Create the initial guess for an optimal tree to set limit on bound
        guess = build_tree_from_list(sorted_clades, conflicts, posclades)
        # Perform branch and bind optimization
        optima = branch_and_bind(sorted_clades, posclades, 
                                 conflicts, [(guess[1], guess[0])], 0, [])
        # Get the amount of time from start to end of optimization
        elapsed = perfc() - start_time
        # Write the results to file
        write_output(optima, posclades, arg, data['expected_failures'],elapsed)
