# Bioinformatics:: Ori Finder
# Ori-Finder, a software tool for finding replication origins in DNA sequences

# Processes through a sequence of DNA and determines how often a certain section of DNA pattern occurs within the inputted sequence.
def PatternCount(Pattern,Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

def HammingDistance(p, q):
    count = 0
    for i in range(0, len(p)):
        if p[i] != q[i]:
            count += 1
    return count

# Generates symbol array of genome based on symbol input. {Position : Count}.
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

# Determines the frequency of all sections of genetic code (length=k) from text and returns a dictionary of all possibilities and their frequencies.
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        count = 0
        for i in range(len(Text)-len(Pattern)+1):
            if Text[i:i+len(Pattern)] == Pattern:
                count = count+1
        freq[Pattern] = count
    return freq

# Creates a list containing all the most frequent k-mers (genetic code sections of k length) in Text.
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words

# Input:  A string Pattern
# Output: The reverse of Pattern
def Reverse(Pattern):
    return Pattern[len(Pattern)::-1]

# Creates the complementary genetic code sequence
def Complement(Pattern):
    basepairs = {"A":"T", "G":"C", "T":"A", "C":"G"}
    complement = ""
    for base in Pattern:
        complement += basepairs.get(base)
    return complement

# Creates the reverse complement of a sequence of code
def ReverseComplement(Pattern):   
    pattern = Reverse(Pattern)
    result = Complement(pattern)
    return result

# Finds all occurrences of a pattern in a genomic sequence Input: Strings Pattern and Genome.
# Output: All starting positions in Genome where Pattern appears as a substring.
def PatternMatching(Pattern, Genome):
    positions = [] 
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions

def ApproximatePatternMatching(Genome, Pattern, d):
    positions = []
    for i in range(len(Genome)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

# Tracking differences between the total number of occurrences of G and the total number of occurrences of C.
def SkewArray(Genome):
    skew = [0]
    score = {"A":0, "T":0, "C":-1, "G":1}
    for i in range(1,len(Genome)+1):
            skew.append(score[Genome[i-1]] + skew[i-1])
    return skew

# Skew array should achieve a minimum at the position where the reverse half-strand ends and the forward half-strand begins, which is the location of ori.  
def MinimumSkew(Genome):
    positions = [] # output variable
    skew = SkewArray(Genome)
    minimum = min(skew)
    for i in range(0,len(skew)):
        if skew[i] == minimum:
            positions.append(i)
    return positions

# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    count = 0
    for i in range(0, len(p)):
        if p[i] != q[i]:
            count += 1
    return count

