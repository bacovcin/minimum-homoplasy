import os

if __name__ == '__main__':

    infile = open('tmp.scores')
    header = infile.readline()

    scores = []
    for line in infile:
        scores.append(line.rstrip().split('\t'))

    infile = open('tmp.trees')
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
        if int(score[2]) > maxcomp:
            maxcomp = int(score[2])

    bestscores = [x for x in scores if int(x[2]) == maxcomp]
    besttrees = [trees[int(x[0])-1] for x in bestscores]

    outfile = open('tmp.trees','w')
    outfile.write(treeheader)
    for tree in besttrees:
        outfile.write(tree)
    outfile.write('End;')
    print('DONE!!!')

    os.remove('tmp.scores')
