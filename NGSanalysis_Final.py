import pyfaidx
from Bio import SeqIO
import os
import csv
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

filenames = ["JDHS03B_S2_L001_","JDHS03D_S4_L001_",]

# Place fastq files into subfolder
# Modify these according to your primers and target site relative to your amplicon
FwdPrimer = "TCACTGTGTGGCTTCAGG"
RevPrimer = "TCGACCCAGCTGTCGGAG"
TargetStart = 28
TargetEnd = 95
#
# # Cycle through files
# for i in range(len(filenames)):

# # Trim off primers from R1
#     trimmed = (rec[TargetStart:TargetEnd] for rec in SeqIO.parse(filenames[i]+"R1.fastq", "fastq"))
#     count = SeqIO.write(trimmed, filenames[i]+"R1_trimmed.fastq", "fastq")
#
# # Make a reverse complement R1 file
#     trimmed_RC = (rec.reverse_complement() for rec in SeqIO.parse(filenames[i]+"R1_trimmed.fastq", "fastq"))
#     count = SeqIO.write(trimmed_RC, filenames[i]+"R1_trimmed_RC.fastq", "fastq")
#     print("Found %i trimmed RC sequences" % count)
#
# # Make a file of the reads in R1 whose reverse complement matches perfectly with part of R2. This removes illumina base-calling errors.
#     count = 0
#     matches = []
#     rcrecords = list(SeqIO.parse(filenames[i]+"R1_trimmed_RC.fastq", "fastq"))
#     for rec in SeqIO.parse(filenames[i]+"R2.fastq", "fastq"):
#         if rcrecords[count].seq in rec.seq:
#             matches.append(rec)
#         count += 1
#     count = SeqIO.write(matches, filenames[i]+"R1_matches.fastq", "fastq")
#     print("Found %i perfectly matching sequences" % count)
#
# # make file containing R1 reads that have either the forward or reverse primer (to remove non target sequences)
#     primer_reads = []
#     for rec in SeqIO.parse(filenames[i]+"R1_matches.fastq", "fastq"):
#         if FwdPrimer in rec.seq:
#                 primer_reads.append(rec)
#         elif RevPrimer in rec.seq:
#                 primer_reads.append(rec)
#     count = SeqIO.write(primer_reads, filenames[i]+"R1_primed_matches.fastq", "fastq")
#     print("Found %i matched sequences with primers" % count)
#
# # make file with all reads oriented in same direction and with the primers trimmed
#     count=0
#     oriented = []
#     for rec in SeqIO.parse(filenames[i]+"R1_matches.fastq", "fastq"):
#         if FwdPrimer in rec.seq:
#                 oriented.append(rec[TargetStart:TargetEnd])
#         elif RevPrimer in rec.seq:
#                 rec = rec[TargetStart:TargetEnd].reverse_complement()
#                 oriented.append(rec)
#     count = SeqIO.write(oriented, filenames[i]+"R1_primed_matches_oriented_trimmed.fastq", "fastq")
#     print("Found %i matched sequences with primers_oriented" % count)
#
# # alignment
#     os.system("bwa index 'pSH0201.fa'")
#     os.system("bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' 'pSH0201.fa' "+filenames[i]+"'R1_primed_matches_oriented_trimmed.fastq' > "+filenames[i]+"'.sam'")
#     os.system("samtools fixmate -O bam "+filenames[i]+"'.sam' "+filenames[i]+"'_fixmate.bam'")
#     os.system("samtools sort -O bam -o "+filenames[i]+"'_sorted.bam' -T '/tmp/201_651_temp' "+filenames[i]+"'_fixmate.bam'")
#     os.system("samtools index "+filenames[i]+"'_sorted.bam'")
#     os.system("samtools mpileup -f 'pSH0201.fa' -q 1 -B "+filenames[i]+"'_sorted.bam' > "+filenames[i]+"'intermediate.mpileup'")
#
# # variant analysis
#     os.system("java -jar VarScan.v2.3.9.jar  mpileup2snp "+filenames[i]+"'intermediate.mpileup' --pileup 1 \
#     --min-coverage 1 \
#     --min-reads2 1 \
#     --variants 1 \
#     --p-value 0.99 \
#     --min-var-freq 0.0003 \
#     --output-vcf 1 > "+filenames[i]+"'snp3.vcf'")
#     os.system("bcftools query -f '%REF\t%ALT[\t%FREQ]\t%POS\n' "+filenames[i]+"'snp3.vcf' > "+filenames[i]+"'snp3.txt'")
#     with open(filenames[i]+'snp3.txt', 'r') as infile, open(filenames[i]+'snp3formatted.txt', 'w') as outfile:
#         temp = infile.read().replace("%", "")
#         outfile.write(temp)

#plot
freq1=[]
pos1=[]
f = open('JDHS03A_S1_L001_snp3formatted.txt')
reader = csv.reader(f,delimiter='\t')
for row in reader:
	freq1.append(row[2])
	pos1.append(row[3])
freq2=[]
pos2=[]
f = open('JDHS03B_S2_L001_snp3formatted.txt')
reader = csv.reader(f,delimiter='\t')
for row in reader:
	freq2.append(row[2])
	pos2.append(row[3])
freq3=[]
pos3=[]
f = open('JDHS03C_S3_L001_snp3formatted.txt')
reader = csv.reader(f,delimiter='\t')
for row in reader:
	freq3.append(row[2])
	pos3.append(row[3])
freq4=[]
pos4=[]
f = open('JDHS03D_S4_L001_snp3formatted.txt')
reader = csv.reader(f,delimiter='\t')
for row in reader:
	freq4.append(row[2])
	pos4.append(row[3])
freq5=[]
pos5=[]
f = open('JDHS03E_S5_L001_snp3formatted.txt')
reader = csv.reader(f,delimiter='\t')
for row in reader:
	freq5.append(row[2])
	pos5.append(row[3])



pos1=map(int,pos1)
pos1[:]=[-1*(x - 583.5) for x in pos1]
pos2=map(int,pos2)
pos2[:]=[-1*(x - 583.5) for x in pos2]
pos3=map(int,pos3)
pos3[:]=[-1*(x - 583.5) for x in pos3]
pos4=map(int,pos4)
pos4[:]=[-1*(x - 583.5) for x in pos4]
pos5=map(int,pos5)
pos5[:]=[-1*(x - 583.5) for x in pos5]

fig1, ax = plt.subplots(figsize=(9,2.7), dpi=100)

fig1=matplotlib.pyplot.scatter(pos2,freq2, color='blue',label = r'$\bf{nCas9-PolI3M on target 1}$')
fig1=matplotlib.pyplot.scatter(pos4,freq4, color='orange', label = r'$\bf{nCas9-PolI3M on target 2}$')
fig1=matplotlib.pyplot.scatter(pos3,freq3, color='green', label = r'$\bf{nCas9 + PolI3M on target}$')
fig1=matplotlib.pyplot.scatter(pos5,freq5, color='purple', label = r'$\bf{nCas9-PolI3M off target}$')
fig1=matplotlib.pyplot.scatter(pos1,freq1, color='red', label = r'$\bf{nCas9 on target}$')
# plt.legend(loc='upper left',frameon=False, fontsize=12)
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.set_xlabel('Distance from nick (bp)',fontsize=12,fontweight='bold')
# ax.set_ylabel('Variant frequency (%)',fontsize=12,fontweight='bold')
# ax.tick_params(axis='x', labelsize=12)
# ax.tick_params(axis='y', labelsize=12)
ax.set_yscale('log')
# plt.ylim([0.01,1.5])
plt.ylim([0.04,2])
plt.xlim([-8,40])
plt.tight_layout()
matplotlib.pyplot.show()
