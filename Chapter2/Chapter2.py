#!/usr/bin/env python
# coding: utf-8

# # Next-Generation Sequencing (NGS)

# In fact, NGS technology consists of a series of methods consisting of preliminary preparation and fragmentation of the studied genome, sequencing, imaging, assembling the sequenced fragments, and data analysis.

# ## Accessing GenBank and moving around NCBI databases

# In[ ]:


from Bio import Entrez, Medline, SeqIO

Entrez.email = "put@your.email.here" 

#This gives you the list of available databases
handle = Entrez.einfo()
rec = Entrez.read(handle)
print(rec)


#Retrieving 20 records
handle = Entrez.esearch(db="nucleotide", term='CRT[Gene Name] AND "Plasmodium falciparum"[Organism]')
rec_list = Entrez.read(handle)

#Get more records
if rec_list['RetMax'] < rec_list['Count']:
    handle = Entrez.esearch(db="nucleotide", term='CRT[Gene Name] AND "Plasmodium falciparum"[Organism]',
                            retmax=rec_list['Count'])
    rec_list = Entrez.read(handle)

#Download all matching nucleotide sequences from GenBank
id_list = rec_list['IdList']
hdl = Entrez.efetch(db='nucleotide', id=id_list, rettype='gb')

#Parse the results
recs = list(SeqIO.parse(hdl, 'gb'))

#Focus on a record
for rec in recs:
    if rec.name == 'KM288867':
        break
print(rec.name)
print(rec.description)

#Extract some sequence information
for feature in rec.features:
    if feature.type == 'gene':
        print(feature.qualifiers['gene'])
    elif feature.type == 'exon':
        loc = feature.location
        print('Exon', loc.start, loc.end, loc.strand)
    else:
        print('not processed:\n%s' % feature)

#Annotations on the record
for name, value in rec.annotations.items():
    print('%s=%s' % (name, value))

#Access the fundamental piece of information(the sequence)
sequence = rec.seq
print(sequence)

#Take all reference annotations
refs = rec.annotations['references']
for ref in refs:
    if ref.pubmed_id != '':
        print(ref.pubmed_id)
        handle = Entrez.efetch(db="pubmed", id=[ref.pubmed_id],
                                rettype="medline", retmode="text")
        records = Medline.parse(handle)
        for med_rec in records:
            for k, v in med_rec.items():
                print('%s: %s' % (k, v))


# ## Performing basic sequence analysis

# In[ ]:


from Bio import Entrez, Seq, SeqIO
from Bio.Alphabet import IUPAC

Entrez.email = "put@your_email.here" 
hdl = Entrez.efetch(db='nucleotide', id=['NM_002299'], rettype='fasta')  # Lactase gene
seq = SeqIO.read(hdl, 'fasta')

#Saving our interest sequence
w_seq = seq[11:5795]
w_hdl = open('example.fasta', 'w')
SeqIO.write([w_seq], w_hdl, 'fasta')
w_hdl.close()

#Reading the sequence
recs = SeqIO.parse('example.fasta', 'fasta')
for rec in recs:
    seq = rec.seq
    print(rec.description)
    print(seq[:10])
    print(seq.alphabet)
    
#Change the alphabet of our sequence
seq = Seq.Seq(str(seq), IUPAC.unambiguous_dna)

#Transcribe the sequence
rna = seq.transcribe()

#Translate the gene
prot = seq.translate()


# ## Working with modern sequence formats

# In[ ]:


from collections import defaultdict

get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
import matplotlib.pyplot as plt

from Bio import SeqIO

recs = SeqIO.parse('SRR003265.filt.fastq', 'fastq')
rec = next(recs)
print(rec.id, rec.description, rec.seq)
print(rec.letter_annotations)


#Distribution of nucleotide reads
cnt = defaultdict(int)
for rec in recs:
    for letter in rec.seq:
        cnt[letter] += 1
tot = sum(cnt.values())
for letter, cnt in cnt.items():
    print('%s: %.2f %d' % (letter, 100. * cnt / tot, cnt))

#Plot the distribution of Ns
n_cnt = defaultdict(int)
for rec in recs:
    for i, letter in enumerate(rec.seq):
        pos = i + 1
        if letter == 'N':
            n_cnt[pos] += 1
seq_len = max(n_cnt.keys())
positions = range(1, seq_len + 1)
fig, ax = plt.subplots()
ax.plot(positions, [n_cnt[x] for x in positions])
ax.set_xlim(1, seq_len)

#Distribution of Phred scores
cnt_qual = defaultdict(int)
for rec in recs:
    for i, qual in enumerate(rec.letter_annotations['phred_quality']):
        if i < 25:
            continue
        cnt_qual[qual] += 1
tot = sum(cnt_qual.values())
for qual, cnt in cnt_qual.items():
    print('%d: %.2f %d' % (qual, 100. * cnt / tot, cnt))

