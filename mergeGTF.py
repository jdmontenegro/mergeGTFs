import pandas as pd
import pyranges as pr
import igraph as ig
import argparse
from collections import defaultdict

class GTFMerger:
    def __init__(self, gtf_files, output_gtf, singleton_gtf, min_shared_locus):
        self.gtf_files = gtf_files
        self.output_gtf = output_gtf
        self.singleton_gtf = singleton_gtf
        self.min_shared_locus = min_shared_locus
        self.exons = []  # List to store parsed exons
        self.exon_dict = {}  # To track unique exons
        self.graph = ig.Graph(directed=True)
        self.node_to_exon = {}
        self.transcript_paths = defaultdict(list)
        self.gene_id_map = {}
        self.path_dict = {}

    def parse_gtf(self):
        for file in self.gtf_files:
            df = pd.read_csv(file, sep='\t', comment='#', header=None,
                             names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'])
            df = df[df['feature'] == 'exon']
            
            for _, row in df.iterrows():
                attrs = self.parse_attributes(row['attributes'])
                exon_key = (row['chr'], row['start'], row['end'], row['strand'])
                transcript_id = attrs.get('transcript_id', 'UNKNOWN')
                
                if exon_key in self.exon_dict:
                    self.exon_dict[exon_key]['transcript_ids'].add(transcript_id)
                else:
                    self.exon_dict[exon_key] = {
                        'chr': row['chr'], 'start': row['start'], 'end': row['end'],
                        'strand': row['strand'], 'transcript_ids': {transcript_id},
                        'locus_id': None
                    }
        
        self.exons = self.merge_start_end_codons(list(self.exon_dict.values()))

    def parse_attributes(self, attr_str):
        attr_dict = {}
        for attr in attr_str.split(';'):
            attr = attr.strip()
            if attr:
                parts = attr.split(' ')
                if len(parts) >= 2:
                    key, value = parts[0], ' '.join(parts[1:])
                    attr_dict[key] = value.strip('"')
        return attr_dict

    def merge_start_end_codons(self, exons):
        merged_exons = []
        exons.sort(key=lambda x: (x['chr'], x['start'], x['end'], x['strand']))
        
        i = 0
        while i < len(exons):
            current_exon = exons[i]
            j = i + 1
            while j < len(exons):
                next_exon = exons[j]
                
                if (current_exon['chr'] == next_exon['chr'] and
                    current_exon['strand'] == next_exon['strand']):
                    
                    if (current_exon['start'] < next_exon['start'] and
                        current_exon['end'] == next_exon['end']):
                        # Merge start codons
                        current_exon['start'] = min(current_exon['start'], next_exon['start'])
                        current_exon['transcript_ids'].update(next_exon['transcript_ids'])
                        j += 1
                    elif (current_exon['end'] > next_exon['end'] and
                          current_exon['start'] == next_exon['start']):
                        # Merge end codons
                        current_exon['end'] = max(current_exon['end'], next_exon['end'])
                        current_exon['transcript_ids'].update(next_exon['transcript_ids'])
                        j += 1
                    else:
                        break
                else:
                    break
            
            merged_exons.append(current_exon)
            i = j
        
        return merged_exons

    def assign_locus_ids(self):
        sorted_exons = sorted(self.exons, key=lambda x: (x['chr'], x['start'], x['end'], x['strand']))
        df = pd.DataFrame(sorted_exons)
        df.rename(columns={'chr': 'Chromosome', 'start': 'Start', 'end': 'End'}, inplace=True)
        loci = pr.PyRanges(df)
        cluster_ids = loci.cluster().df['Cluster'].tolist()
        for i, exon in enumerate(sorted_exons):
            exon['locus_id'] = cluster_ids[i]

    def build_graph(self):
        node_idx = {}
        edges = []
        
        for i, exon in enumerate(self.exons):
            node_idx[(exon['chr'], exon['start'], exon['end'], exon['strand'])] = i
            self.graph.add_vertex(name=str(i), label=f"{exon['chr']}:{exon['start']}-{exon['end']}({exon['strand']})")
            self.node_to_exon[i] = exon
            for tid in exon['transcript_ids']:
                self.transcript_paths[tid].append(i)
        
        for tid, path in self.transcript_paths.items():
            path.sort(key=lambda x: (self.node_to_exon[x]['chr'], self.node_to_exon[x]['start']))
            path_tuple = tuple(path)
            
            if path_tuple in self.path_dict:
                continue  # Skip duplicate paths
            
            self.path_dict[path_tuple] = tid
            for j in range(len(path) - 1):
                edges.append((path[j], path[j + 1], tid))
                
        self.graph.add_edges([(e[0], e[1]) for e in edges])
        self.graph.es['transcript_id'] = [e[2] for e in edges]

    def assign_gene_ids(self):
        gene_clusters = self.graph.connected_components(mode='weak')
        
        for gid, cluster in enumerate(gene_clusters):
            gene_id = f'GENE_{gid}'
            for node in cluster:
                self.gene_id_map[node] = gene_id
    def write_gtf(self):
        with open(self.output_gtf, 'w') as f:
            for path, transcript_id in sorted(self.path_dict.items(), key=lambda x: (self.node_to_exon[x[0][0]]['chr'], self.node_to_exon[x[0][0]]['start'])):
                gene_id = self.gene_id_map[path[0]]
                f.write(f"{self.node_to_exon[path[0]]['chr']}\tMerged\ttranscript\t{self.node_to_exon[path[0]]['start']}\t{self.node_to_exon[path[-1]]['end']}\t.\t{self.node_to_exon[path[0]]['strand']}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                for node in path:
                    exon = self.node_to_exon[node]
                    f.write(f"{exon['chr']}\tMerged\texon\t{exon['start']}\t{exon['end']}\t.\t{exon['strand']}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")

    def run(self):
        self.parse_gtf()
        self.assign_locus_ids()
        self.build_graph()
        self.assign_gene_ids()
        self.write_gtf()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge GTF files and assign gene and locus IDs.")
    parser.add_argument("gtf_files", nargs='+', help="Input GTF files.")
    parser.add_argument("--output", default="merged.gtf", help="Output GTF file (default: merged.gtf)")
    parser.add_argument("--singleton", default="singleton_exons.gtf", help="Output file for singleton exons (default: singleton_exons.gtf)")
    parser.add_argument("--min_shared_locus", type=int, default=2, help="Minimum shared locus IDs before assigning the same gene ID (default: 2)")
    
    args = parser.parse_args()
    merger = GTFMerger(args.gtf_files, args.output, args.singleton, args.min_shared_locus)
    merger.run()

