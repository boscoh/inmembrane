import sys
import os
import glob


names = """
AE004092 Streptococcus pyogenes 
AE017198 Lactobacillus johnsonii
P000033 Lactobacillus acidophilus 
CP000413 Lactobacillus gasseri 
CR954253 Lactobacillus delbrueckii 
"""
tokens = [l.split() for l in names.splitlines()]
gbk = { ts[0]: ' '.join(ts[1:]) for ts in tokens if ts}


for genome in gbk:
    pse = genome + '.pse'
    if not os.path.isfile(pse):
      os.system("python ../../inmembrane.py %s.fasta > %s" % (genome, pse))

    if base in gbk:
        print "%s: %s" % (base, gbk[base])
    print "--"
    
    category_count = {}
    for line in open(pse):
        if line.startswith("#"):
            continue
        tokens = line.split()
        category = tokens[1]
        if category not in category_count:
            category_count[category] = 0
        category_count[category] += 1

    n = sum(category_count.values())

    category_count["TOTAL"] = n
    for category in ['CYTOPLASM', 'MEMBRANE', 'PSE', 'SECRETED', 'TOTAL']:
        count = category_count[category]
        percent = 100.0*count/float(n)
        print "%-9s %4d %3.f%%" % (category, count, percent)

    print
    print
