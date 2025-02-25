from genome_tookit import *

sequence = "AAAGAAAATGTGCTTCGAAATCGTA"
target_pattern = "TGT"
target_kmer = "AA"

gt = DNAPatternAnalysis(sequence, target_pattern, 1)
print("DNA Pattern Analysis:")
print(f"[1] Approximate matches found at index: {gt.approx_matches()}")
print(f"[2] Exact matches found at: {gt.find_all_occurrences_re()}")
print(f"[3] Introns found at: {gt.find_intron()}")
print(f"[4] K-mer Count: {gt.count_kmer(target_kmer)}")
print(f"[5] Most Frequent K-mer: {gt.frequent_kmers(3)}")
