import sys
import os 
import statistics

from Bio import SeqIO

def parsing():
    if len(sys.argv) != 2:
        print("Usage: python bakta.py <file.gbff>")
        sys.exit(1)

    file = sys.argv[1]

    if not os.path.exists(file):
        print("path/file doesnt exist")
        sys.exit(1)

    return file

def file_process(file):
    file_list = list(SeqIO.parse(file,'genbank'))
    genes = []
    for record in file_list:
        for feature in record.features:
            if feature.type == 'CDS':
                genes.append(feature)

    unknown = 0
    known = 0
    lengths = []

    for gene in genes:
        length = gene.location.end - gene.location.start
        lengths.append(length)
       

       #not to ommit some protein qualities 
        if 'product' in gene.qualifiers:
            product = gene.qualifiers['product'][0].lower()
            if any(x in product for x in ['hypothetical', 'putative', 'uncharacterized', 'predicted']):
                unknown += 1
            else:
                known += 1
        else:
            unknown += 1
    
    median_len = statistics.median(lengths)
    
    return {
        'genes': genes,
        'unknown': unknown,
        'known': known,
        'median_len': median_len
    }

def print_results(data):
    print("-" * 30)
    print(f"Total genes: {len(data['genes'])}")
    print(f"Genes with unknown function: {data['unknown']}")
    print(f"Genes with known function: {data['known']}")
    print(f"Median gene length: {data['median_len']:.0f} bp")

def main():
    print("\n" + "=" * 30)
    print("ANNOTATION INFO")
    print("=" * 30)
    
    file = parsing()
    
    print(f"\nAnalyzing: {os.path.basename(file)}")
    analysis = file_process(file)
    print_results(analysis)

if __name__ == "__main__":
    main()
