#!/usr/bin/env python

import sys
import gzip
import re
import subprocess 

if len(sys.argv) != 3:
    sys.stderr.write("Usage:  gatkLowQualityFilter ref.fasta calls.vcf\n")
    sys.exit(-1)

ref = sys.argv[1]
vcffile = sys.argv[2]

debug = False

lastchr = ""
lastseq = ""

faillines = []
if vcffile.endswith(".gz"):
    fp = gzip.open(vcffile, "rt")
else:
    fp = open(vcffile)
for line in fp:
    if line.startswith("#"):
        continue

    f = line.rstrip("\n").split("\t")
    if f[6] != "PASS":
        faillines.append(f)
        continue


    match = re.search(";MQ=([^;]+);", f[7])
    if match is not None and float(match.group(1)) < 30.0:
        faillines.append(f)
        continue

    match = re.search(";QD=([^;]+);", f[7])
    if match is not None and float(match.group(1)) < 2.0:
        faillines.append(f)
        continue

    f2 = f[8].split(":")
    adidx = f2.index("AD")
    gtidx = f2.index("GT")
    gqidx = f2.index("GQ")
    plidx = f2.index("PL")

    goodcnt = 0
    majoritycnt = 0
    for i in range(9, len(f)):
        f2 = f[i].split(":")
        if adidx >= len(f2) or (f2[adidx] != "." and sum([ int(x) for x in f2[adidx].split(",") ]) == 0):
            continue

        ad = [ int(x) for x in f2[adidx].split(",") ]
        refcnt = ad[0]
        altcnt = sum(ad[1:])

        if len(f[3]) == 1 and len(f[4]) == 1:
            if altcnt >= 3 and altcnt * 1.0 / (altcnt + refcnt) >= 0.2:
                goodcnt += 1
        else:
            if altcnt >= 4 and altcnt * 1.0 / (altcnt + refcnt) >= 0.25:
                goodcnt += 1

        if altcnt > refcnt:
            majoritycnt += 1

    if goodcnt == 0:
        faillines.append(f)
        continue
        
    badplcnt = 0
    homcnt = 0
    for i in range(9, len(f)):
        f2 = f[i].split(":")
        if adidx >= len(f2) or (f2[adidx] != "." and sum([ int(x) for x in f2[adidx].split(",") ]) == 0):
            continue

        if plidx < len(f2) and sum([ 1 for x in f2[plidx].split(",") if x == "0" ]) > 1:
            badplcnt += 1

        gt = f2[gtidx]
        if gt[0] == "." or gt[2] == ".":
            continue

        if gt[0] == gt[2] and gt[0] != '0':
            homcnt += 1

    if badplcnt > 2 and badplcnt * 4 > len(f) - 9:
        faillines.append(f)
        continue

    chr = f[0]
    pos = int(f[1])

    if chr != lastchr:
        sys.stderr.write("Reading %s..." % chr)
        sys.stderr.flush()

        cmd = "samtools faidx %s %s" % (ref, chr)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=-1, universal_newlines=True)
        chrfp = p.stdout
    
        l = []
        for x in chrfp:
            if x.startswith(">"): continue
            l.append(x.strip().upper())
        chrfp.close()
        p.wait()

        lastseq = "".join(l)
        lastchr = chr

        sys.stderr.write("done\n")
        sys.stderr.flush()

    cnts = [ 0, 0, 0, 0 ]
    for ch in lastseq[pos:pos+15]:
        if ch == 'A':
            cnts[0] += 1
        elif ch == 'C':
            cnts[1] += 1
        elif ch == 'G':
            cnts[2] += 1
        elif ch == 'T':
            cnts[3] += 1

    if homcnt * 5 < len(f) - 9 and max(cnts) >= 13 and abs(len(f[3]) - len(f[4])) <= 1 and majoritycnt == 0:
        faillines.append(f)
        continue
fp.close()

failset = {}
for f in faillines:
    qual = float(f[5])
    match = re.search(";QD=([^;]+);", f[7])
    if match is None:
        #sys.stderr.write("\t".join(f) + "\n")
        continue
    qd = float(match.group(1))

    if qd >= 20.0:
        continue

    gtidx = f[8].split(":").index("GT")
    plidx = f[8].split(":").index("PL")
    pllist = []
    for t in f[9:]:
        f2 = t.split(":")
        if len(f2) > plidx and f2[gtidx][0] != '.' and (f2[gtidx][0] != '0' or f2[gtidx][2] != '0'):
            pllist.append(int(f2[plidx].split(",")[0]))

    if len(pllist) > 0 and sum(pllist) / len(pllist) >= 800:
        continue

    key = ":::".join([ f[0], f[1], f[3], f[4] ])
    failset[key] = f[6]

if vcffile.endswith(".gz"):
    fp = gzip.open(vcffile, "rt")
else:
    fp = open(vcffile)
for line in fp:
    if line.startswith("#"):
        sys.stdout.write(line)
        continue

    f = line.rstrip("\n").split("\t")
    
    key = ":::".join([ f[0], f[1], f[3], f[4] ])
    if key in failset:
        if failset[key] == "PASS":
            f[6] = "LowQual"
        else:
            f[6] = failset[key]
    else:
        f[6] = "PASS"

    sys.stdout.write("\t".join(f) + "\n")
