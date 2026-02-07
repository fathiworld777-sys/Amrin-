from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
print("All successful")

#Step : Sequence Quality Analysis

record = SeqIO.read("yfgA.fasta", "fasta")
print("The Sequence ID:", record.id)
print("The Sequence Description:", record.description)
print("Sequence:", record.seq)
print("Type of the seq record:", type(record))

print("The Length of the Amino acid:", len(record))

# Amino acid composition 

sequence = str(record.seq)

aa_count = {}

for aa in sequence:
    aa_count[aa] = aa_count.get(aa, 0) + 1

print("Amino acid composition:")
for aa, count in sorted(aa_count.items()):
    print(f"{aa}: {count}")

#Step 3: Sequence Filtering & Validation

if len(record.seq) < 50:
    print("Too short sequence")
else:
    print("Sequence is fine")


#Step 4: Homology Search (BLAST)
"""

result_handle = NCBIWWW.qblast("blastp", "nr", record.seq)

with open("blast_results.xml", "w") as out_b:
    out_b.write(result_handle.read())

"""
print("BLAST search complete. Results saved to blast_results.xml")

with open("blast_results.xml") as b:
    blast_record = NCBIXML.read(b)


top_hit = blast_record.alignments[0]
print("Top BLAST hit:")
print(top_hit.hit_def)


top_hit = blast_record.alignments[0]
top_hsp = top_hit.hsps[0]

#closest homolog

print("Closest homolog:")
print("Protein:", top_hit.hit_def)
print("Length:", top_hit.length)
print("E-value:", top_hsp.expect)
print("Identity:", top_hsp.identities, "/", top_hsp.align_length)
print("Percent identity:", (top_hsp.identities / top_hsp.align_length) * 100)

#Identify Conserved Regions

print("Query sequence:")
print(top_hsp.query)

print("Match:")
print(top_hsp.match)

print("Subject sequence:")
print(top_hsp.sbjct)

#Step 5: Functional Annotation

from Bio import SeqIO
record = SeqIO.read("yfgA.fasta", "fasta")
print(record.id)
print(record.description)
print(len(record.seq))
print(record.seq)

#Predict: Function, Biological role, Organism relevance

#1)The protein YfgA is involved in cell shape determination and regulation of cell wall synthesis in bacteria. 
# Specifically, it's a cytoskeletal protein that helps control the length of the long axis of the cell, maintaining its rod shape. 
# It's also thought to contribute to the regulation of penicillin-binding proteins

#2)Biological Role: YfgA acts in bacterial cell wall synthesis and cell shape maintenance. 
# It's part of the cytoskeletal system regulating cell morphology.

#3)Organism Relevance:In bacteria like E. coli, YfgA is important for maintaining cell shape and structural integrity,
#  crucial for survival and division. Disruption can lead to shape changes or lysis.

#Step 6: Biological Interpretation (Research Outcome)

#What does this sequence likely do?

#This sequence is for a hypothetical protein YfgA in E. coli. Likely function is unclear, 
# but it might be involved in cellular processes like metabolism or regulation, based on sequence motifs and homology.

#Why do they think so?

#Predictions are based on:
# Sequence motifs and domains (e.g., potential binding sites).
# Homology to other proteins with known functions.
# Gene context in E. coli genome (e.g., operon structure).
#But, experimental data is limited for YfgA, making it a hypothetical protein

#What evidence supports their claim?

#For YfgA, evidence is limited since it's hypothetical. Support for predicted functions might include:
# Sequence similarity to proteins with known domains.
# Gene neighborhood or operon context in E. coli.
# Structural predictions suggesting possible binding or catalytic roles.
#But, concrete experimental data (e.g., knockout studies, assays) is needed for confirmation.