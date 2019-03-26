
"""
This program reads in a BLASTn output file called 'align'.xml
It displays alignment as well as BLASTn results.
"""

# >>> pip install xmltodict
import xmltodict

with open('./align.xml', 'r') as infile:
    raw_xml = infile.read()

xml = xmltodict.parse(raw_xml)
alignment = xml['BlastOutput']['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit']
hit_def = alignment['Hit_def']
seq_a = alignment['Hit_hsps']['Hsp']['Hsp_qseq']
seq_b = alignment['Hit_hsps']['Hsp']['Hsp_hseq']
mid_line = alignment['Hit_hsps']['Hsp']['Hsp_midline']
score = int(alignment['Hit_hsps']['Hsp']['Hsp_score'])
e = float(alignment['Hit_hsps']['Hsp']['Hsp_evalue'])
identities = int(alignment['Hit_hsps']['Hsp']['Hsp_identity'])
gaps = int(alignment['Hit_hsps']['Hsp']['Hsp_gaps'])

chunk_size = 120
print('\tBLASTN Alignment:\n')
for i in range(int(len(seq_a) / chunk_size) + 1):
    start = i * chunk_size
    stop = start + chunk_size
    index = start + 1
    print(index, '\t', seq_a[start:stop])
    print(
        ''.join([' ' for num in str(index)]), '\t', 
        mid_line[start:stop]
    )
    print(index, '\t', seq_b[start:stop], '\n')

print(hit_def)
print(f"Score: {score} bits")
print(f"Expect: {e}")
print(f"Identities: {identities}/{len(seq_a)} - {int(identities/len(seq_a) * 100)}%")
print(f"Gaps {gaps}/{len(seq_a)} - {int(gaps / len(seq_a))}%")