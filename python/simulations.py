from random import randint as rand
import sedef

#global vars
minLen = 1000 #min sd length
maxLen = 100000 #max sd length
maxSED = 15 #maximum SNP/single indel error (perecnt)
maxLED = 15 #Maximum gap error (percent)

letter = 'ATCGATCG'

def loadSeq(filename):
    """get sequence from filename"""
    with open(filename, 'r') as f:
        seq = ''.join(l.strip() for l in f) # seq is your seq now
    return seq

def randSeq(length):
    """ make a random seq of specifies length """
    arrout = []
    for i in range(length):
        arrout.append(letter[rand(0,3)])

    return ''.join(arrout)

def makeSmall(sequence, error):
    """ modifies a sequence with single inserts, deletes, and SNPs """
    def SNP(a):
        # returns a random base different from input
        return letter[letter.find(a) + rand(1,3)]

    smallED = error
    loc = 0
    arrout = []
    while loc < len(sequence):
        action = rand(1, 100)
        if action <= smallED//3: # delete
            pass
        elif action <= 2*smallED//3: # insert
            arrout.append(letter[rand(0,3)])
            arrout.append(sequence[loc])
        elif action <= smallED: # SNP
            arrout.append(SNP(sequence[loc]))
        else:
            arrout.append(sequence[loc])
        loc += 1

    return (''.join(arrout), smallED)

def makeLarge(sequence, error):
    """ modifies a sequence with large indels """
    def noIntersection(start, end):
        # returns true if start and end position do not intersect with a deleted or outside region
        if end > length:
            return False
        for i in inserts:
            if start <= i[0] <= end:
                return False
        for d in deletes:
            if start <= d[0] <= end or d[0] <= start <= d[0] + d[1]:
                return False
        return True

    length = len(sequence)
    maxLargeED = error*length/100
    largeED = 0
    arrout = []
    inserts = [] # contains (location, length) of insert
    deletes = [] # contains (location, length) of delete
    data = (inserts, deletes)

    counter = 0 # randomly generate indels
    while maxLargeED > 50 and counter < 10:
        counter += 1
        gapLen = rand(50, maxLargeED)
        action = rand(0,1)
        location = rand(0, length)
        if noIntersection(location, location + action*gapLen):
            maxLargeED -= gapLen
            largeED += gapLen
            data[action].append((location, gapLen))

    inserts.sort()
    deletes.sort()

    iloc = 0
    ilen = len(inserts)

    dloc = 0
    dlen = len(deletes)

    loc = 0
    tempIns = inserts + [(float('inf'), float('inf'))]
    tempDel = deletes + [(float('inf'), float('inf'))]

    while iloc != ilen and dloc != dlen: # create the sequence with the random indels
        if tempIns[iloc][0] < tempDel[dloc][0]:
            arrout.append(sequence[loc:tempIns[iloc][0]])
            arrout.append(randSeq(tempIns[iloc][1]))
            loc = tempIns[iloc][0]
            iloc += 1
        elif tempIns[iloc][0] > tempDel[dloc][0]:
            arrout.append(sequence[loc:tempDel[dloc][0]])
            loc = tempDel[dloc][0] + tempDel[dloc][1]
            dloc += 1
        else:
            print('error')
            print(data)
    arrout.append(sequence[loc:])
    return (''.join(arrout), (largeED, data))

def generate_random_sd(error, seq = None):
    """ generates random sd with error% error rate
    If seq is specified, random sd  is generated from a substring of it."""
    if seq == None:
        seq1 = randSeq(rand(minLen, maxLen))
    else:
        length = rand(minLen, maxLen)
        start = rand(0, len(seq) - length - 1)
        seq1 = seq[start, start + length]
    sED = rand(max(0, error - maxLED),min(maxSED, error))
    seq2 = makeSmall(seq1, sED)[0]
    seq2 = makeLarge(seq2, error-sED)[0]
    return seq1, seq2, sED

#####above, generating sequences. below, running sedef on them and outputting in tables ####

def combine(arr):
    """ makes overlapping sequences 1 sequence """
    def first(item):
        return item[0]
    def second(item):
        return item[1]

    if len(arr) == 0 or len(arr) == 1:
        return arr

    sarr = []
    for c, val in enumerate(arr):
        sarr.append((val[0], val[1], c))

    sarr = sorted(sarr, key = second)
    sarr = sorted(sarr, key = first)

    chains = [[sarr[0][0], sarr[0][1], [sarr[0][2]]]]

    for s, e, c in sarr[1:]: #start, end, counter
        if s <= chains[-1][1] +1:
            chains[-1][1] = max(e, chains[-1][1])
            chains[-1][2].append(c)
        else:
            chains.append([s, e, [c]])

    return chains

def calculateSum(arr):
    """input array of starts and ends, output length. """

    if len(arr) == 0:
        return 0
    elif len(arr) == 1:
        return arr[0][1] - arr[0][0]

    out = 0
    for s, e, c in combine(arr):
        out += e - s

    return out

def resultsTable(runs, maxerror = 30, output = "output.txt", seq = None, freeroom = 20):
    """ outputs results in tables. WIll simulate # on runs for each values between 0 and maxerror. inclusive.
    output is output file.  Can specify that a sequence should be a subsequence of seq.
    Can specify what percentage of free rooms to call hits vs. partials  """
    filename = "detailed-"+ output
    array = []
    badNum = 0
    misses  = 0
    with open(filename, 'w') as f:
        f.write('sep=;\n')
        f.write('length1; length2; aveLen; TED; SED; LED; jSum1; jSum2; jMax1; jMax2; '+
            'cSum1; cSum2; cMax1; cMax2; cAlig; jcSum1; jcSum2; jcMax1; jcMax2; jcAlig; ' +
            'jSumVal; jMaxVal; cSumVal; cMaxVal; jcSumVal; jcMaxVal; jcSum012\n')

        for error in range(maxerror + 1):
            out = [0,0,0] # hit, miss, partial
            for i in range(runs):
                print('***',i,'***')

                seq1, seq2, sED = generate_random_sd(error, seq)
                len1 = len(seq1)
                len2 = len(seq2)
                alen = float(len1 + len2)/2

                jSum = [0,0]
                jMax = [0,0]
                cSum = [0,0]
                cMax = [0,0]
                aMax = [0,0]
                jcSum = [0,0]
                jcMax = [0,0]


                temp1 = []
                temp2 = []


                aln = sedef.PyAligner()

                hits1 = aln.chain_align(seq1.upper(), seq2.upper())

                print('chains')
                for h in hits1:
                    print('  {}..{} -> {}..{}  (sz={} g={} m={})'.format(h.query_start(), h.query_end(), h.ref_start(), h.ref_end(),
                        h.alignment_size(), h.gaps(), h.mismatches()))
                    temp1.append((h.query_start(), h.query_end()))
                    temp2.append((h.ref_start(), h.ref_end()))
                    cMax[0] = max(cMax[0],  h.query_end() -  h.query_start())
                    cMax[1] = max(cMax[1],  h.ref_end() -  h.ref_start())
                    aMax[0] = max(aMax[0], h.alignment_size())
                cSum = [calculateSum(temp1), calculateSum(temp2)]


                temp1 = []
                temp2 = []


                hits2 = aln.jaccard_align(seq1.upper(), seq2.upper())

                print('jaccard')
                for h in hits2:
                    print('  {}..{} -> {}..{} === {} (sz={} g={} m={})'.format(h.query_start(), h.query_end(), h.ref_start(), h.ref_end(),
                        h.cigar(), h.alignment_size(), h.gaps(), h.mismatches()))

                    temp1.append((h.query_start(), h.query_end()))
                    temp2.append((h.ref_start(), h.ref_end()))
                    jMax[0] = max(jMax[0],  h.query_end() -  h.query_start())
                    jMax[1] = max(jMax[1],  h.ref_end() -  h.ref_start())
                jSum = [calculateSum(temp1), calculateSum(temp2)]


                print('chain on jaccard')

                def extend(query_start, query_end, ref_start, ref_end):
                    """extend each jaccard hit to length 15000 or 10 times itself"""
                    w = max(query_end - query_start, ref_end - ref_start);
                    w = min(15000, int(5 * w));
                    return (max(0, query_start - w), min(query_end + w, len1), max(0, ref_start - w), min(ref_end + w, len2))

                for i in range(len(temp1)):
                    extended = extend(temp1[i][0], temp1[i][1], temp2[i][0], temp2[i][1])
                    temp1[i] = (extended[0], extended[1])
                    temp2[i] = (extended[2], extended[3])

                combined = (combine(temp1), combine(temp2))

                if len(combined[0]) == 1 and len(combined[1]) == 1:
                    temp1 = []
                    temp2 = []

                    hits3 = aln.chain_align(seq1[combined[0][0][0]:combined[0][0][1]], seq2[combined[1][0][0]:combined[1][0][1]])
                    for h in hits3:
                        print('  {}..{} -> {}..{}  (sz={} g={} m={})'.format(h.query_start(), h.query_end(), h.ref_start(), h.ref_end(),
                            h.alignment_size(), h.gaps(), h.mismatches()))
                        temp1.append((h.query_start(), h.query_end()))
                        temp2.append((h.ref_start(), h.ref_end()))
                        jcMax[0] = max(jcMax[0],  h.query_end() -  h.query_start())
                        jcMax[1] = max(jcMax[1],  h.ref_end() -  h.ref_start())
                        aMax[1] = max(aMax[1], h.alignment_size())
                    jcSum = [calculateSum(temp1), calculateSum(temp2)]

                elif len(combined[0]) == 0 or len(combined[1]) == 0:
                    misses += 1

                else: #non overlappping thingsd
                    groups = [[],[]]

                    print('*(87E98008343247@&^@(&*#&^#^NOT IMPLEMENTED@#&(#^&#@#$@#(&^@#*)&^)#(&*^')
                    print(len(combined[0]))
                    print(len(combined[1]))
                    badNum += 1
                    # should probably figure out how to make this case work

                #'jSumVal; jMaxVal; cSumVal; cMaxVal; jcSumVal; jcMaxVal; jcSum012
                if jcSum[0]+ jcSum[1] == 0:
                    val = 0
                    out[1] += 1
                elif 50*float(jcSum[0]+ jcSum[1])/alen < 100 - freeroom - error:
                    vel = 1
                    out[2] += 1
                else:
                    val = 2
                    out[0] += 1

                with open(filename, 'a') as f:
                    f.write('{}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; \n'.format(len1,
                    len2, alen, error, sED, error-sED,
                    jSum[0], jSum[1], jMax[0], jMax[1], cSum[0], cSum[1],
                    cMax[0], cMax[1], aMax[0], jcSum[0], jcSum[1], jcMax[0], jcMax[1], aMax[1],
                    float(jSum[0]+ jSum[1])/alen, float(jMax[0]+ jMax[1])/alen, float(cSum[0]+ cSum[1])/alen,
                    float(cMax[0]+ cMax[1])/alen, float(jcSum[0]+ jcSum[1])/alen, float(jcMax[0]+ jcMax[1])/alen, val))

            array.append((error, out))

        print(misses, " misses, ", badNum, ' unexpected things')
        f.write(str(misses) + " misses, " + str(badNum) + ' unexpected things\n')

    with open(output, 'w') as f:
        f.write('sep=;\n')
        f.write('error; hits; misses; partials\n')
        for e, r in array:
            f.write('{}; {}; {}; {}; \n'.format(e, *r))
        if badNum > 0:
            f.write(str(misses) + " misses, " + str(badNum) + ' unexpected things\n')


if __name__ == '__main__':
    resultsTable(1000)

