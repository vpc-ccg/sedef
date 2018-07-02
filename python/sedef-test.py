# 786

import sedef

seq1 = 'ATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCATTTTCCTGAGAGGAAAGCTTTCCCACATTATACAGCTTCTGAAAGGGTTGCTTGACCCACAGATGTGAAGCTGAGGCTGAAGGAGACTGATGTGGTTTCTCCTCAGTTTCTCTGTGTGGCACCAGGTGGCAGCAGAGGTCAGCAAGGCAAACCCGAGCCCAGGGATGCGGGGTGGGGGCAGGTACATCCTCTCTTGAGCTACAGCAGATTAACTCTGTTCTGTTTCATTGTGGTTGTTTAGTTTGCGTTTTTTTTTCTCCAACTTTGTGCTTCATCGGGAAAAGCTTTGGATCACAATTCCCAGtgctgaagaaaaggccaaactctggaaaaaatttgaatattttgagccaaatgtgaggaccacaacctgtgagaacggaaaataaatcctgggaccccagactcactaagccaaagggaaaagccaagctgggaactggcttatgcaaacctgcttcccatctggttcctaaataagatagctattacacaaagacaaaaaagctacatccctgcctctacctccatcgcatgcaaaatgtgtattcagtgaacgctgaccaaagacagaagaatgcaaccatttgcctctgatttacccacacccattttttccacttcttcccctttccccaatacccgcacttttcccctttacttactgaggtccccagacaacctttgggaaaagcacggaccacagtttttcctgtggttctctgttcttttctcaggtgtgtccttaaccttgcaaatagatttcttgaaatgattgagactcaccttggttgtgttctttgattAGTgcctgtgacgcagcttcaggaggtcctgagaacgtgtgcacagtttagtcggcagaaacttagggaaatgtaagaccaccatcagcacataggagttctgcattggtttggtctgcattggtttggtctggaaggaggaaaattcaaagtaatggggcttacaggtcatagatagattcaaagattttctgattgtcaattggttgaaagaattattatctacagacctgctatcaatagaaaggagagtctgggttaagataagagactgtggagacc'
seq2 = 'ATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCATTTTCCTGAGAGGAAAGCTTTCCCACATTATTCAGCTTCTGAAAGGGTTGCTTGACCCACAGATGTGAAGCTGAGGCTGAAGGAGACTGATGTGGTTTCTCCTCAGTTTCTCTGTGCGGCACCAGGTGGCAGCAGAGGTCAGCAAGGCAAACCCGAGCCCGGGGATGCGGGGTGGGGGCAGCTACGTCCTCTCTTGAGCTACAGCAGATTCACTCTGTTCTGTTTCATTGTTGCTTAGTTTGCGTTTTGTTTCTCCAACTTTGTGCCTCATCAGGAAAAGCTTTGGATCACAATTCCCAGtgctgaagaaaaggccaaactctggaaaaaattttgaatattttgagccaaatgtgaggaccacaacctgtgagaacggaaaataaatcctgggaccccagactcactaagccaaagggaaaagccaagctgggaactggcttatgcaaacctgcttcccatctggttcctaaataagatagctattacacaaagataaaaaagctacatccctgcctctacctccctcgcatgtaaaatgtgtattcagtgaacactgaccaaagacagaagaatgcaaccatttgcctctgatttacccacacccattttttccacttcttcccctttccccaatacccgcacttttcccctttacttactgaggcccccagacaatctttgggaaaagcacggaccacagtttttcctgtggttctctgttcttttctcaggtgtgtccttaaccttgcaaatagatttcttgaaatgattgacactcaccttggttgtgttctttgatcagcgcctgtgacgcagcttcaggaggtcctgagaacgtgtgcacagtttagtcggcagaaacttagggaaacgtaagaccaccatcagtacgtaggagttgtgcattggtttggtctggaaggaggaaaattcaaagtaatggggcttacaggtcatagatagattcaaagattttctgattgtcaattgattgaaagaattattatctacagacctgctatcaatagaaaggagagtctgagttaagataagagactgtggagacc'

aln = sedef.PyAligner()
hits = aln.chain_align(seq1.upper(), seq2.upper())
print 'chains'
for h in hits:
	print '  {}..{} -> {}..{} === {} (sz={} g={} m={})'.format(h.query_start(), h.query_end(), h.ref_start(), h.ref_end(), 
		h.cigar(), h.alignment_size(), h.gaps(), h.mismatches())

hits = aln.jaccard_align(seq1.upper(), seq2.upper())
print 'jaccard'
for h in hits:
	print '  {}..{} -> {}..{} === {} (sz={} g={} m={})'.format(h.query_start(), h.query_end(), h.ref_start(), h.ref_end(), 
		h.cigar(), h.alignment_size(), h.gaps(), h.mismatches())

hits = aln.full_align(seq1.upper(), seq2.upper())
print 'full'
for h in hits:
	print '  {}..{} -> {}..{} === {} (sz={} g={} m={})'.format(h.query_start(), h.query_end(), h.ref_start(), h.ref_end(), 
		h.cigar(), h.alignment_size(), h.gaps(), h.mismatches())

