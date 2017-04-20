"""
Compare analysed models using the stats.json files (which contain the evidence,
but also marginal summaries)
"""
import json
import sys
import numpy
from math import log, log10

prefixes = sys.argv[1:]

models = dict([(f, json.load(open(f if f.endswith('stats.json') else f + "stats.json"))['global evidence']) for f in prefixes])

best = max(models, key=models.__getitem__)
Zbest = models[best]


for m in models: 
	models[m] -= Zbest
Ztotal = log(sum(numpy.exp([Z for Z in models.values()])))
limit = 30 # for example, Jeffreys scale for the Bayes factor

print
print 'Model comparison'
print '****************'
print
for m in sorted(models, key=models.__getitem__):
	Zrel = models[m]
	print 'model %-25s: log10(Z) = %6.1f %s' % (m.replace('stats.json',''), Zrel / log(10),
		' XXX ruled out' if Zrel < Ztotal - log(limit) else '   <-- GOOD' )

print
print 'The last, most likely model was used as normalization.'
print 'Uniform model priors are assumed, with a cut of log10(%s) to rule out models.' % limit
print

try:
	AICs = dict([(f, 2*len(json.load(open(f.replace('stats.json', 'params.json')))) + numpy.loadtxt(f.replace('stats.json', '') + '.txt')[:,1].min()) for f in prefixes])

	fout = open('modelcompare.rst', 'w')
	fout.write('+-%-40s-+-%-15s-+-%-15s-+\n' % ('-'*40, '-'*15, '-'*15))
	fout.write('| %-40s | %+15s | %+15s |\n' % ('model', 'ln Z', 'AIC'))
	#print '+-%-40s-+-%-15s-+-%-15s-+ ' % ('-'*40, '-'*15, '-'*15)
	fout.write('+=%-40s=+=%-15s=+=%-15s=+\n' % ('='*40, '='*15, '='*15))
	for m in sorted(models, key=models.__getitem__):
		mname = m.replace('stats.json','').rstrip('_')
		Zrel = models[m]
		fout.write('| %-40s | %+15s | %+15s |\n' % (mname, 
			'%.1f' % Zrel, '%.1f' % AICs[m]))
		fout.write('+-%-40s-+-%-15s-+-%-15s-+\n' % ('-'*40, '-'*15, '-'*15))
except IOError as e:
	print 'Warning:', e
        pass
	

