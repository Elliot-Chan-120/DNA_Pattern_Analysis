# goals
# search for kmer with max hamming distance set to 2 in a sequence A
# generate complementary sequence B from the sequence A
# reverse the sequence B (generate reverse complement)
# search for kmer with max hamming distance set to 2 in sequence B
# collect display graph save results

from re import search, finditer

class DNAPatternAnalysis:
    """
    DNA Patterns: Approx. matching, All matches, Intron finder
    Default values: ATCG, AT, 2
    """
    def __init__(self, sequence= "ATCG", find_pattern= "AT", hamming_distance = 2):
        """initialization"""
        self.sequence = sequence.upper()
        self.find_pattern = find_pattern.upper()
        self.h_d = hamming_distance

    # DNA Pattern Recognition Toolkit Functions
    def h_d_zip(self, substring):
        """
        Hamming distance calculation w/ zip method
        \nGo to approx_matches, this just outputs True or False
        \nGets substring from approx_matches
        """
        zipped_dna = zip(substring, self.find_pattern)
        if len([(n1, n2) for n1, n2 in zipped_dna if n1 != n2]) <= self.h_d:
            return True
        else:
            return False

    def approx_matches(self):
        """
        Finds approx matches for the DNA seq and target pattern
        \n error margin = hamming distance -> how many mismatches we're willing to look past
        """
        matches = []
        pattern_length = len(self.find_pattern)

        for i in range(len(self.sequence) - pattern_length + 1):
            substring = self.sequence[i:i + pattern_length]
            if self.h_d_zip(substring):
                matches.append(i)
        if not matches:
            return None
        return matches


    def find_all_occurrences_re(self):
        """
        Finds all occurrences of the target pattern in the sequence
        :return:
        """
        pattern = f"(?="+self.find_pattern+")"
        mos = finditer(pattern, self.sequence)
        results = []
        for motif in mos:
            results.append(motif.span()[0])
        if len(results) > 0:
            return results
        else:
            return None

    def find_intron(self):
        intron = "GT[ATCG]{1,10}TACTAAC[ATCG]{1,10}AC"
        result = search(intron, self.sequence)
        if result != None:
            return result.span()[0]
        else:
            return None

    #K-MER ANALYSIS
    def count_kmer(self, kmer):
        """
        Counts repeating k-mers in a sequence. Includes overlapping k-mers
        **MUST INPUT**: kmer -> what specific kmer you're looking for
        """
        kmer_count = 0
        for position in range(len(self.sequence) - (len(kmer) - 1)):
            if self.sequence[position:position + len(kmer)] == kmer:
                kmer_count += 1
        return kmer_count

    def frequent_kmers(self, k_len):
        """
        Finds the most frequent k-mers of a given length in a DNA string
        **MUST INPUT**: k_len (int) = length of k-mers to search for
        returns list of most frequent k-mers in DNA string
        """
        # dict to store k-mer frequencies
        kmer_frequencies = {}

        # loop iterates through DNA string and extracts k-mers of given length and adds frequency
        for i in range(len(self.sequence) - k_len + 1):
            kmer = self.sequence[i:i+k_len]
            if kmer in kmer_frequencies:
                kmer_frequencies[kmer] += 1
            else:
                kmer_frequencies[kmer] = 1

        highest_frequency = max(kmer_frequencies.values())

        return [
            kmer for kmer, frequency in kmer_frequencies.items()
            if frequency == highest_frequency
        ]
