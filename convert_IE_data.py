import sys
def isInt(x):
    try:
        int(x)
        return True
    except:
        return False

if __name__ == '__main__':
    weights = {'P':'1','M':'1','L':'1'}
    exceptioninh = {}
    exclusions = []
    for arg in sys.argv:
        s = arg.split('=')
        if s[0] == 'P' and isInt(s[1]):
            weights['P'] = str(s[1])
        elif s[0] == 'M' and isInt(s[1]):
            weights['M'] = str(s[1])
        elif s[0] == 'L' and isInt(s[1]):
            weights['L'] = str(s[1])
        elif s[0] == 'exinh':
            exceptioninh[s[1]] = s[2]
        elif s[0] == 'exc':
            exclusions.append(s[1])
        elif s[0] == 'specw':
            weights[s[1]] = s[2]

    infile = open('IEDATA')
    names = infile.readline().rstrip().split(' ')
    chars = []
    weightings = []
    inhval = []
    for line in infile:
        try:
            s = line.rstrip().split(' ')
            if s[1] in exclusions:
                continue
            if s[1][0] == 'P':
                try:
                    weightings.append(weights[s[1]])
                except KeyError:
                    weightings.append(weights['P'])
                try:
                    inhval.append(exceptioninh[s[1]])
                except KeyError:
                    inhval.append('1')
            elif s[1][0] == 'M':
                try:
                    weightings.append(weights[s[1]])
                except KeyError:
                    weightings.append(weights['M'])
                try:
                    inhval.append(exceptioninh[s[1]])
                except KeyError:
                    inhval.append('-1')
            else:
                try:
                    weightings.append(weights[s[1]])
                except KeyError:
                    weightings.append(weights['L'])
                try:
                    inhval.append(exceptioninh[s[1]])
                except KeyError:
                    inhval.append('-1')
            charname = s[0]+s[1]
            if '.' in s[0] and s[0].split('.')[0] == prevnum:
                news = s[2:]
                for i in range(len(prev)):
                    if news[i] == prev[i]:
                        news[i] = 'NA'
                chars.append((charname,news))
            else:
                prev = s[2:]
                prevnum = s[0].split('.')[0]
                chars.append((charname,s[2:]))
        except IndexError:
            pass

    if len(list(weights.keys())) == 3:
        outfile = open('rawdata/IE/IE-convert-'+
                       weights['L']+
                       weights['P']+
                       weights['M']+
                       '.csv','w')
    else:
        outfile = open('rawdata/IE/IE-convert-special.csv','w')
    header = ','
    for char in chars:
        header += char[0]+','
    outfile.write(header.rstrip(',')+'\n')
    outfile.write('Weighting,'+','.join(weightings)+'\n')
    outfile.write('InheritedValue,'+','.join(inhval)+'\n')

    for i in range(len(names)):
        nextline = names[i]+','
        for char in chars:
            nextline += char[1][i] + ','
        outfile.write(nextline.rstrip(',')+'\n')
