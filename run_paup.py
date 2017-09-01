import string
from Bio import Phylo
import os
import sys
import timeit

def convertWeightedCharactersToNexus(infile,outfile,keyfile,
                                     outfilename, outputname):
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
        for i in range(len(names[1:])):
            keyfile.write(names[i+1]+':\n')
            for key in cstates[i]:
                keyfile.write('\t'+cstates[i][key]+':'+key+'\n')
        tmpfile = open('temp1.nex','w')
        tmpfile.write('BEGIN PAUP;\n')
        tmpfile.write('  execute '+outfilename+';\n')
        tmpfile.write('  set criterion=parsimony maxtrees=100 increase=no;\n')
        tmpfile.write('  hsearch start=stepwise addseq=random ')
        tmpfile.write('nreps=25 swap=tbr;\n')
        tmpfile.write('  filter best=yes;\n')
        tmpfile.write('  set maxtrees=100 increase=no;\n')
        tmpfile.write('  hsearch start=current swap=tbr hold=1 nbest=1000;\n')
        tmpfile.write('  pscores all/ ci ri rc hi compatibility=yes ')
        tmpfile.write('scorefile=tmp.scores replace=yes;\n')
        tmpfile.write('  savetrees file=tmp.trees replace=yes ')
        tmpfile.write('format=nexus;\n')
        tmpfile.write('  QUIT;\n')
        tmpfile.write('END;\n')
        tmpfile = open('temp2.nex','w')
        tmpfile.write('BEGIN PAUP;\n')
        tmpfile.write('  execute '+outfilename+';\n')
        tmpfile.write('  gettrees file=tmp.trees;\n')
        tmpfile.write('  treeinfo;\n')
        tmpfile.write('  pscores all/ ci ri rc hi compatibility=yes ')
        tmpfile.write('scorefile='+outputname+
                      '_comp.scores replace=yes;\n')
        tmpfile.write('  savetrees file='+outputname+
                      '_comp.trees replace=yes ')
        tmpfile.write('format=newick;\n')
        tmpfile.write('  contree /treeFile='+outputname+
                      '_comp_con.trees replace=yes semistrict=yes;\n')
        tmpfile.write('  QUIT;\n')
        tmpfile.write('END;')
        tmpfile = open('temp3.nex','w')
        tmpfile.write('BEGIN PAUP;\n')
        tmpfile.write('  execute '+outfilename+';\n')
        tmpfile.write('  ASSUME wtset=mywts ancstates=ancestor;\n')
        tmpfile.write('  set criterion=parsimony maxtrees=100 increase=no;\n')
        tmpfile.write('  hsearch start=stepwise addseq=random ')
        tmpfile.write('nreps=25 swap=tbr;\n')
        tmpfile.write('  filter best=yes;\n')
        tmpfile.write('  set maxtrees=100 increase=no;\n')
        tmpfile.write('  hsearch start=current swap=tbr ')
        tmpfile.write('hold=1 nbest=100;\n')
        tmpfile.write('  filter best=yes;\n')
        tmpfile.write('  pscores all/ ci ri rc hi compatibility=yes ')
        tmpfile.write('scorefile='+outputname+
                      '_wp.scores replace=yes;\n')
        tmpfile.write('  savetrees file='+outputname+
                      '_wp.trees replace=yes ')
        tmpfile.write('format=newick;\n')
        tmpfile.write('  contree /treeFile='+outputname+
                      '_wp_con.trees replace=yes semistrict=yes;\n')
        tmpfile.write('  QUIT;\n')
        tmpfile.write('END;\n')
        tmpfile = open('temp4.nex','w')
        tmpfile.write('BEGIN PAUP;\n')
        tmpfile.write('  gettrees file='+outputname+'_comp_con.trees;\n')
        tmpfile.write('  savetrees file='+outputname+
                      '_comp_con.trees replace=yes ')
        tmpfile.write('format=newick;\n')
        tmpfile.write('  QUIT;\n')
        tmpfile.write('END;\n')
        tmpfile = open('temp5.nex','w')
        tmpfile.write('BEGIN PAUP;\n')
        tmpfile.write('  gettrees file='+outputname+'_wp_con.trees;\n')
        tmpfile.write('  savetrees file='+outputname+
                      '_wp_con.trees replace=yes ')
        tmpfile.write('format=newick;\n')
        tmpfile.write('  QUIT;\n')
        tmpfile.write('END;\n')
    return

if __name__ == "__main__":
    # Load in the datasets
    for arg in sys.argv:
        # Only csv files are valid datasets
        if arg[-3:] != 'csv':
            continue
        print('Creating nexus file...')
        s = '/'.join(arg.split('.')[0].split('/')[1:])
        outfilename = 'nexfiles/' + s + '.nex'
        outputname = 'outputs/' + s 
        keyfilename = 'nexfiles/' + s + '.key'
        with open(outfilename,'w') as outfile:
            with open(keyfilename,'w') as keyfile:
                convertWeightedCharactersToNexus(open(arg),
                                                 outfile,
                                                 keyfile,
                                                 outfilename,
                                                 outputname)
        os.system('paup temp1.nex &&'+
                  'python paup_maximum_compatibility.py &&'+
                  'paup temp2.nex &&'+
                  'paup temp3.nex &&'+
                  'paup temp4.nex &&'+
                  'paup temp5.nex')
        os.remove('tmp.trees')
        os.remove('temp1.nex')
        os.remove('temp2.nex')
        os.remove('temp3.nex')
        os.remove('temp4.nex')
        os.remove('temp5.nex')
        trees = list(Phylo.parse(outputname+'_comp.trees','newick'))
        with open(outputname+'_comp.trees','w') as outfile:
            i = 0
            for tree in trees:
                i += 1
                outfile.write('Tree #'+str(i)+':\n')
                Phylo.draw_ascii(tree,file=outfile)
        trees = list(Phylo.parse(outputname+'_comp_con.trees','newick'))
        with open(outputname+'_comp_con.trees','w') as outfile:
            i = 0
            for tree in trees:
                i += 1
                outfile.write('Tree #'+str(i)+':\n')
                Phylo.draw_ascii(tree,file=outfile)
        trees = list(Phylo.parse(outputname+'_wp.trees','newick'))
        with open(outputname+'_wp.trees','w') as outfile:
            i = 0
            for tree in trees:
                i += 1
                outfile.write('Tree #'+str(i)+':\n')
                Phylo.draw_ascii(tree,file=outfile)
        trees = list(Phylo.parse(outputname+'_wp_con.trees','newick'))
        with open(outputname+'_wp_con.trees','w') as outfile:
            i = 0
            for tree in trees:
                i += 1
                outfile.write('Tree #'+str(i)+':\n')
                Phylo.draw_ascii(tree,file=outfile)
