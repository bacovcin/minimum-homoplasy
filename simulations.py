import multiprocessing as mp
from fit_data import *
from copy import copy
from random import shuffle, random, randint, sample
import string, os

def build_dicttree_from_list(clist):
    listtree = []
    for clade in clist:
        if clade == frozenset():
            continue
        if False not in [testTreeSuit(clade,x) for x in listtree]:
            listtree.append(clade)
    new_tree = sorted(listtree, key=lambda x:-len(x))
    tree_dict = {}
    for clade in new_tree:
        tree_dict = create_tree_dict(tree_dict,clade)
    return tree_dict

def getLChar(tree,curval,availables):
    output = {}
    for key in tree:
        q = random()
        if q < 0.2:
            newval = availables.pop(0)
        else:
            newval = curval
        if tree[key] == {}:
            output[key] = newval
        else:
            newdict = getLChar(tree[key],newval,availables)
            for key2 in newdict:
                output[key2] = newdict[key2]
    return output

def generate_tree(taxa):
    posclades = [frozenset([x]) for x in taxa]+[frozenset(taxa)]
    for i in range(1000):
        posclades.append(frozenset(sample(taxa,randint(2,len(taxa)-1))))
    return build_dicttree_from_list(posclades)

def flattenDictTree(tree):
    output = []
    for key in tree:
        output.append(key)
        if tree[key] != {}:
            output += flattenDictTree(tree[key])
    return output

def generate_data(ntaxa, nchar, inhrate):
    data = {}
    data['lgg'] = string.ascii_uppercase[:ntaxa]
    data['chars'] = {}
    data['inhchar'] = {}
    data['weights'] = {}
    data['expected_failures'] = 0
    tree = generate_tree(data['lgg'])
    data['tree'] = flattenDictTree(tree)
    for c in range(nchar):
        cname = 'char' + str(c)
        data['chars'][cname] = {}
        lchar = getLChar(tree,0,list(range(1,100)))
        for l in lchar:
            try:
                data['chars'][cname][lchar[l]] += [list(l)[0]]
            except KeyError:
                data['chars'][cname][lchar[l]] = [list(l)[0]]
        q = random()
        if q < inhrate:
            data['inhchar'][cname] = 0
            data['expected_failures'] -= 1
        else:
            data['inhchar'][cname] = -1
        data['weights'][cname] = 1
    for char in data['chars']:
        for state in data['chars'][char]:
            if state != data['inhchar'][char]:
                data['expected_failures'] += data['weights'][char]
    return data

def getConfusionMatrix(tree1,tree2):
    set1 = set(tree1)
    set2 = set(tree2)
    confmat = {}
    confmat['actual_clade_num'] = len(set1)
    confmat['tp'] = len(set1 & set2)
    confmat['fn'] = len(set1 - set2)
    confmat['fp'] = len(set2 - set1)
    return confmat

def run_simulation(mysim):
    data = generate_data(mysim['ntaxa'],mysim['nchar'],mysim['inhrate'])
    starttime = perfc()
    posclades = find_possible_clades(data)
    mysim['time-findclades'] = perfc() - starttime
    starttime = perfc()
    mysim['possible_clade_num']=len(posclades.keys())
    conflicts = find_conflicts(posclades.keys())
    mysim['time-findconflicts'] = perfc() - starttime
    starttime = perfc()
    sorted_clades = sort_clades(posclades.keys(),conflicts,posclades)
    mysim['time-sorting'] = perfc() - starttime
    starttime = perfc()
    optima = branch_and_bind(sorted_clades, posclades, 
                             conflicts, 
                             [(data['expected_failures'],
                               sorted_clades)],
                             0, [], False)
    mysim['time-bandb'] = perfc() - starttime
    mysims = [copy(mysim)]
    treenum = 1
    mysims[-1]['treenum'] = treenum

    confmat = getConfusionMatrix(data['tree'],optima[0][1])
    for key in confmat:
        mysims[-1][key] = confmat[key]
    consensustree = set(optima[0][1])
    for optimum in optima[1:]:
        consensustree &= set(optimum[1])
        mysims.append(copy(mysim))
        treenum += 1
        if treenum > 1:
            mysims[-1]['time-findclades'] = 'NA'
            mysims[-1]['time-findconflicts'] = 'NA'
            mysims[-1]['time-sorting'] = 'NA'
            mysims[-1]['time-bandb'] = 'NA'
        mysims[-1]['treenum'] = treenum
        confmat = getConfusionMatrix(data['tree'],optimum[1])
        for key in confmat:
            mysims[-1][key] = confmat[key]
    else:
        if treenum > 1:
            mysims.append(copy(mysim))
            mysims[-1]['time-findclades'] = 'NA'
            mysims[-1]['time-findconflicts'] = 'NA'
            mysims[-1]['time-sorting'] = 'NA'
            mysims[-1]['time-bandb'] = 'NA'
            mysims[-1]['treenum'] = 'consensus'
            confmat = getConfusionMatrix(data['tree'],list(consensustree))
            for key in confmat:
                mysims[-1][key] = confmat[key]

    return mysims

def simworker(queue):
    while True:
        mysim = queue.get(True)
        print('Worker '+ str(os.getpid())+
              ' running simulation #' + str(mysim['simnum']) + '; '+
              'NTaxa: '+str(mysim['ntaxa'])+'; '+
              'NChar: '+str(mysim['nchar'])+'...')
        mysim['pid'] = os.getpid()
        newsims = run_simulation(mysim)
        for sim in newsims:
            with open('simulations/simulations.csv','a') as outfile:
                outfile.write(','.join([str(sim[x])
                                        for x in names])+'\n')
        queue.task_done()

if __name__ == '__main__': 
    names = ['pid',
             'simnum',
             'treenum',
             'ntaxa',
             'nchar',
             'inhrate',
             'time-findclades',
             'time-findconflicts',
             'time-sorting',
             'time-bandb',
             'tp',
             'fn',
             'fp',
             'actual_clade_num',
             'possible_clade_num']
    with open('simulations/simulations.csv','w') as outfile:
        outfile.write(','.join(names)+'\n')
    simulation_queue = mp.JoinableQueue()
    workers = mp.Pool(6,simworker,(simulation_queue,))
    counter = 0
    for l in range(100):
        for j in [600,50,300,100,500,200,400]:
            for i in [25,6,8,21,9,10,17,12,11,14,13,7,5]:
                counter += 1
                simulation_queue.put({'simnum':counter,
                                      'ntaxa':i,
                                      'nchar':j,
                                      'inhrate':0.1})
    with open('counter','w') as counterfile:
        counterfile.write(str(counter))
    simulation_queue.join()
