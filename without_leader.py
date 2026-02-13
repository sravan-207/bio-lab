from Bio import Entrez, SeqIO

Entrez.email = "sravanvakkanti@gmail.com"


accession = "WP_000516135.1"

handle = Entrez.efetch(db="protein",
                       id=accession,
                       rettype="fasta",
                       retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

seq = str(record.seq)
print("Full Sequence:")
print(seq)

pep = seq[:5]
print("\nFirst 5 Amino Acids:", pep)

mass = {
    'G':57,'A':71,'S':87,'P':97,'V':99,
    'T':101,'C':103,'I':113,'L':113,
    'N':114,'D':115,'K':128,'Q':128,
    'E':129,'M':131,'H':137,'F':147,
    'R':156,'Y':163,'W':186
}

def total_mass(p):
    return sum(mass[x] for x in p)

parent = total_mass(pep)
print("Parent Mass:", parent)

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

s = cyclic_spec(pep)
print("Cyclic Spectrum:")
print(s)
