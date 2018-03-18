import string
from Bio import Phylo
import os
import sys
import timeit

def weighted_compatibility_filter(treefile,scorefile,weights):
    infile = open(scorefile)
    header = infile.readline()

    scores = []
    curtree = []
    for line in infile:
        if line[0] != '\t' and curtree != []:
            scores.append((len(scores),sum(curtree)))
            curtree = []
        else:
            s = line.rstrip().split('\t')
            compscore = int(s[-1])
            charid = int(s[-2])-1
            curtree.append(compscore*weights[charid])

    infile = open(treefile)
    treeheader = ''
    trees = []
    treestart = False
    for line in infile:
        if line[:4] != 'tree' and not treestart:
            treeheader += line
        else:
            treestart = True
            trees.append(line)

    maxcomp = 0
    for score in scores:
        if score[1] > maxcomp:
            maxcomp = score[1]

    bestscores = [x for x in scores if x[1] == maxcomp]
    besttrees = [trees[x[0]] for x in bestscores]

    outfile = open('tmp.trees','w')
    outfile.write(treeheader)
    for tree in besttrees:
        outfile.write(tree)
    outfile.write('End;')
    print('DONE!!!')

def addDataToPAUPFile(infilename,outfilename,newline):
    with open(infilename) as infile:
        with open(outfilename,'w') as outfile:
            outfile.write(infile.readline())
            outfile.write(newline)
            for line in infile:
                outfile.write(line)

def convertWeightedCharactersToNexus(infile,outfile,
                                     outfilename):
    # Convert the data from our format into a Nexus file
    # For weighted characters, repeat the character a number of times
    # equal to its weight

    # Save the alphabet for naming states
    alphabet = string.ascii_uppercase

    # Extract the names of the characters
    names = infile.readline().rstrip().split(',')

    # Create a dictionary with the character data
    lgg = {}

    # Fill in the dictionary row by row
    for line in infile:
        # Split the line by commas
        s = line.rstrip().split(',')
        lgg[s[0].replace(' ','_').replace("'",'')] = s[1:]

    # Get a list of the taxa
    taxa = [x for x in lgg.keys() if x not in ['Weighting','InheritedValue']]

    # Get the character names and possible state labels
    cstates = []
    astates = []
    for taxum in [x for x in lgg.keys() if x != 'Weighting']:
        for i in range(len(lgg[taxum])):
            newstate = lgg[taxum][i]
            try:
                try:
                    lgg[taxum][i] = cstates[i][newstate]
                except KeyError:
                    cstates[i][newstate] = alphabet[astates[i]]
                    lgg[taxum][i] = cstates[i][newstate]
                    astates[i] += 1
            except IndexError:
                cstates.append({'NA':'-','N/A':'-'})
                astates.append(0)
                try:
                    lgg[taxum][i] = cstates[i][newstate]
                except KeyError:
                    cstates[i][newstate] = alphabet[astates[i]]
                    lgg[taxum][i] = cstates[i][newstate]
                    astates[i] += 1

    # Write the data to a Nexus file
    outfile.write('#NEXUS\n')
    outfile.write('BEGIN TAXA;\n')
    outfile.write('  DIMENSIONS NTAX='+str(len(taxa))+';\n')
    nextstring = '  TAXLABELS\n'
    j = 1
    for taxum in taxa:
        nextstring += "[" +str(j)+"] '" + taxum + "'\n"
        j+=1
    outfile.write(nextstring.rstrip()+'\n;\n')
    outfile.write('END;\n')
    outfile.write('BEGIN CHARACTERS;\n')
    outfile.write('  DIMENSIONS NCHAR='+str(len(names)-1)+';\n')
    outfile.write('  FORMAT DATATYPE=Standard GAP=- SYMBOLS="'+
                  ''.join(alphabet)+'";\n')

    outfile.write('  MATRIX')
    for taxum in taxa:
        outfile.write("\n    '"+taxum+"' "+''.join(lgg[taxum]))
    else:
        outfile.write(';\nEND;\n')
    if 'Weighting' in lgg.keys() or 'InheritedValue' in lgg.keys():
        outfile.write('BEGIN ASSUMPTIONS;\n')
        try:
            outfile.write('  WTSET mywts (VECTOR)='+
                          ' '.join(lgg['Weighting'])+';\n')
        except KeyError:
            pass
        try:
            outfile.write('  ANCSTATES ancestor (VECTOR)='+
                          ''.join(lgg['InheritedValue'])+';\n')
        except KeyError:
            pass
        outfile.write('END;\n')
    addDataToPAUPFile('pauptemplates/template1.nex',
                      'temp1.nex',
                      'EXECUTE '+outfilename+';\n')
    addDataToPAUPFile('pauptemplates/template2.nex',
                      'temp2.nex',
                      'EXECUTE '+outfilename+';\n')
    addDataToPAUPFile('pauptemplates/template3.nex',
                      'temp3.nex',
                      'EXECUTE '+outfilename+';\n')
    addDataToPAUPFile('pauptemplates/template4.nex',
                      'temp4.nex',
                      'EXECUTE '+outfilename+';\n')
    addDataToPAUPFile('pauptemplates/template5.nex',
                      'temp5.nex',
                      'EXECUTE '+outfilename+';\n')
    return [int(x) for x in lgg['Weighting']]

if __name__ == "__main__":
    # Load in the datasets
    for arg in sys.argv:
        # Only csv files are valid datasets
        if arg[-3:] != 'csv':
            continue
        print('Creating nexus file...')
        outdir = '/'.join(arg.split('.')[0].split('/')[1:-1])
        filename = arg.split('.')[0].split('/')[-1]
        outfilename = 'nexfiles/' + outdir+'/'+filename + '.nex'
        with open(outfilename,'w') as outfile:
            weights = convertWeightedCharactersToNexus(open(arg),
                                                       outfile,
                                                       outfilename)
        os.system('paup temp1.nex')
        weighted_compatibility_filter('tmp.trees','tmp.scores',weights) 
        os.system('paup temp2.nex &&'+
                  'paup temp3.nex &&'+
                  'paup temp4.nex &&'+
                  'paup temp5.nex')
        os.remove('tmp_comp_con.trees')
        os.remove('tmp_wp_con.trees')
        os.remove('tmp.trees')
        os.remove('temp1.nex')
        os.remove('temp2.nex')
        os.remove('temp3.nex')
        os.remove('temp4.nex')
        os.remove('temp5.nex')
        trees = list(Phylo.parse('out_comp.trees','newick'))
        with open('outputs/'+outdir+'/comp/'+filename+'.txt','w') as outfile:
            i = 0
            for tree in trees:
                i += 1
                Phylo.draw_ascii(tree,file=outfile)
        os.remove('out_comp.trees')
        trees = list(Phylo.parse('out_wp.trees','newick'))
        with open('outputs/'+outdir+'/wmp/'+filename+'.txt','w') as outfile:
            i = 0
            for tree in trees:
                i += 1
                outfile.write('Tree '+str(i)+':\n')
                Phylo.draw_ascii(tree,file=outfile)
        os.remove('out_wp.trees')
        os.rename('wp_scores.txt','outputs/'+outdir+'/wmp/'+filename+'_scores.txt')
        trees = list(Phylo.parse('out_wp_con.trees','newick'))
        with open('outputs/'+outdir+'/wmp/'+filename+'_con.txt','w') as outfile:
            i = 0
            for tree in trees:
                i += 1
                outfile.write('Tree '+str(i)+':\n')
                Phylo.draw_ascii(tree,file=outfile)
        os.remove('out_wp_con.trees')
        os.remove('tmp.scores')
