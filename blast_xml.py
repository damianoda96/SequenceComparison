
"""
Interface for BLAST-XML files.
"""

# >>> pip install xmltodict
import xmltodict

class BlastHit():
    def __init__(self, hit):
        self.seq_a = hit['Hsp_qseq']
        self.seq_b = hit['Hsp_hseq']
        self.mid_line = hit['Hsp_midline']
        self.score = int(hit['Hsp_score'])
        self.e = float(hit['Hsp_evalue'])
        self.identities = int(hit['Hsp_identity'])
        self.gaps = int(hit['Hsp_gaps'])
    
    def display(self):
        chunk_size = 120
        for i in range(int(len(self.seq_a) / chunk_size) + 1):
            start = i * chunk_size
            stop = start + chunk_size
            index = start + 1
            print(index, '\t', self.seq_a[start:stop])
            print(
                ''.join([' ' for num in str(index)]), '\t', 
                self.mid_line[start:stop]
            )
            print(index, '\t', self.seq_b[start:stop], '\n')

        print(f"Score: {self.score} bits")
        print(f"Expect: {self.e}")
        print(f"Identities: {self.identities}/{len(self.seq_a)} - {int(self.identities/len(self.seq_a) * 100)}%")
        print(f"Gaps {self.gaps}/{len(self.seq_a)} - {int(self.gaps / len(self.seq_a))}%")


class BlastXML():
    def __init__(self, fname):
        with open(fname, 'r') as infile:
            raw_xml = infile.read()

        self.xml = xmltodict.parse(raw_xml)
        self.title = self.xml['BlastOutput']['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit']['Hit_def']
        self.hits = []
        unparsed_hits = self.xml['BlastOutput']['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit']['Hit_hsps']['Hsp']
        if type(unparsed_hits) == list:
            for hit in unparsed_hits:
                self.hits.append(BlastHit(hit))
        else:
            self.hits.append(BlastHit(unparsed_hits))
            

    def display(self):
        for i, (hit) in enumerate(self.hits):
            print('-' * 20)
            print(f"Hit #{i+1}\n")
            hit.display()
        print('-' * 20)

if __name__ == '__main__':
    print('BLAST-N RESULTS:')
    blast = BlastXML('./blastn.xml')
    blast.display()

    print('MEGA-BLAST RESULTS:')
    blast = BlastXML('./megablast.xml')
    blast.display()
