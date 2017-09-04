import sys

if __name__ == '__main__':
    for arg in sys.argv:
        if arg[-3:] in ('nex','nxs'):
            infile = open(arg)
            header = infile.readline().rstrip()
            if header != '#NEXUS':
                continue
            blocks = {}
            curblock = []
            curline = ''
            for line in infile:
                print(line)
                if line[:5] == 'BEGIN':
                    curblock = [line.split(' ')[1].rstrip(';\n'),
                                []]
                    curline = ''
                elif line[:3] == 'END':
                    blocks[curblock[0]] = curblock[1]
                elif line[-2] == ';':
                    curline += line[:-2].lstrip()
                    curblock[-1].append(curline.upper())
                    curline = ''
                else:
                    curline += line[:-1].rstrip().lstrip()+' '
            try:
                blocks['DATA'] = blocks['TAXA'] + blocks['CHARACTERS']
            except KeyError:
                pass

            dim = {}
            taxa = []
            myformat = {}
            matrix = {}
            for line in blocks['DATA']:
                s=line.rstrip().lstrip().split(' ')
                if s[0] == 'DIMENSIONS':
                    for x in s:
                        if x.split('=')[0] == 'NTAX':
                            dim['taxa'] = int(x.split('=')[1])
                        elif x.split('=')[0] == 'NCHAR':
                            dim['char'] = int(x.split('=')[1])
                elif s[0] == 'TAXLABELS':
                    for x in s[1:]:
                        if x[0] == '[':
                            continue
                        else:
                            taxa.append(x)
                elif s[0] == 'FORMAT':
                    for x in s[1:]:
                        y = x.split('=')
                        myformat[y[0]] = y[1]
                elif s[0] == 'MATRIX':
                    i = 1
                    while i < len(s):
                        matrix[s[i]] = list(s[i+1])
                        try:
                            matrix[s[i]] = [x.replace(myformat['GAP'],'NA')
                                            for x in matrix[s[i]]]
                        except KeyError:
                            pass
                        i += 2
            weightings = ['1'] * dim['char']
            inhchar = ['-1'] * dim['char']
            try:
                for line in blocks['ASSUMPTIONS']:
                    s = line.rstrip().lstrip().split(' ')
                    if s[0] == 'WTSET':
                        if 'VECTOR' in line:
                            weightings = line.rstrip().split('=')[1].split(' ')
                        else:
                            s = line.rstrip().split('=')[1].split(',')
                            for x in s:
                                y = s.split(':')
                                for z in y[1].split(' '):
                                    if '-' in z:
                                        a = z.split('-')
                                        for myint in range(int(a[0]),
                                                           int(a[1])+1):
                                            weightings[myint-1] = y[0]
                                    else:
                                        weightings[int(z)-1] = y[0]
                    elif s[0] == 'ANCSTATES':
                        if 'VECTOR' in line:
                            inhchar = list(line.rstrip().split('=')[1])
                        else:
                            s = line.rstrip().split('=')[1].split(',')
                            for x in s:
                                y = s.split(':')
                                for z in y[1].split(' '):
                                    if '-' in z:
                                        a = z.split('-')
                                        for myint in range(int(a[0]),
                                                           int(a[1])+1):
                                            inhchar[myint-1] = y[0]
                                    else:
                                        inhchar[int(z)-1] = y[0]
            except KeyError:
                pass
            outfilename = 'rawdata/converts/'+arg.split('/')[1].split('.')[0]+'.csv'
            outfile = open(outfilename,'w')
            outfile.write(','.join(['char'+str(x) for x in
                                    range(dim['char']+1)])+'\n')
            outfile.write('Weighting,'+','.join(weightings)+'\n')
            outfile.write('InheritedValue,'+','.join(inhchar)+'\n')
            for taxum in matrix:
                outfile.write(taxum+','+','.join(matrix[taxum])+'\n')
