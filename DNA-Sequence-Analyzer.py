from Bio import SeqIO
from Bio.Seq import Seq

# les Fichiers présents
import os
print("Fichiers présents :", os.listdir())

# Lire la séquence depuis un fichier FASTA
record = SeqIO.read("DNA_Seq.fasta", "fasta")
sequence = str(record.seq)
print("Séquence lue :", sequence)

# Longueur de la séquence
print("Longueur :", len(sequence))

# Contenu GC
gc_count = sequence.count("G") + sequence.count("C")
gc_content = (gc_count / len(sequence)) * 100
print("Contenu GC :", gc_content, "%")

# Transcription en ARN
dna_seq = Seq(sequence)
rna_seq = dna_seq.transcribe()
print("ARN :", rna_seq)

# Ajuster la longueur pour qu'elle soit multiple de 3
trimmed_seq = dna_seq[:len(dna_seq) - len(dna_seq) % 3]

# Traduction
protein_seq = trimmed_seq.translate()
print("Protéine :", protein_seq)

# Nettoyage pour éviter des acides aminés après le stop
protein_seq = trimmed_seq.translate(to_stop=True)
print("Protéine :", protein_seq)
