from Bio import Entrez, SeqIO
from collections import Counter

Entrez.email = "sravanvakkanti@gmail.com"


accession = "WP_000516135.1"

handle = Entrez.efetch(db="protein",
                       id=accession,
                       rettype="fasta",
                       retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

seq = str(record.seq)
pep = seq[:5]

mass = {
    'G':57,'A':71,'S':87,'P':97,'V':99,
    'T':101,'C':103,'I':113,'L':113,
    'N':114,'D':115,'K':128,'Q':128,
    'E':129,'M':131,'H':137,'F':147,
    'R':156,'Y':163,'W':186
}

aa = list(mass.keys())

def total_mass(p):
    return sum(mass[x] for x in p)

def cyclic_spec(p):
    pre = [0]
    for i in range(len(p)):
        pre.append(pre[i] + mass[p[i]])

    full = pre[-1]
    spec = [0]

    for i in range(len(p)):
        for j in range(i+1, len(p)+1):
            spec.append(pre[j] - pre[i])
            if i > 0 and j < len(p):
                spec.append(full - (pre[j] - pre[i]))

    return sorted(spec)

exp_spec = cyclic_spec(pep)
parent = total_mass(pep)

def score(p, spec):
    theo = cyclic_spec(p)
    return sum((Counter(theo) & Counter(spec)).values())

def expand(lst):
    return [x + y for x in lst for y in aa]

def leaderboard(spec, N):
    board = [""]
    leader = ""

    while board:
        board = expand(board)
        new = []

        for p in board:
            if total_mass(p) == parent:
                if score(p, spec) > score(leader, spec):
                    leader = p
            if total_mass(p) <= parent:
                new.append(p)

        board = sorted(new,
                       key=lambda x: score(x, spec),
                       reverse=True)[:N]

    return leader

best = leaderboard(exp_spec, 10)

print("Original Peptide:", pep)
print("Best Peptide Found:", best)
