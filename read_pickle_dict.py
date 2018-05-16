import cPickle as pickle
import sys 

fin = pickle.load(open(sys.argv[1], 'rb'))

print fin

#for element in fin.keys():
#	print element, fin[element]