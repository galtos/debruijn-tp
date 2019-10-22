import argparse
import networkx as nx

def read_fastq(path_fastq):
    """
    read fastq file and output a list of the reads
    """
    with open(path_fastq,'r') as filin:
        for ligne in filin:
            yield next(filin)
            next(filin)
            next(filin)

def cut_kmer(seq,k_mer):
    """
    take a sequence and a k-mer size and return an iterator of k-mer
    """
    for i in range(len(seq)-k_mer):
        yield seq[i:(i+k_mer)]

def build_kmer_dict(path_fastq, k_mer):
    """
    take a fastq file in argument and a k mer size and return a dictionary with the
    k-mer as keys and the number of occurences as values
    """
    dict_k_mer = {}
    for seq_read in read_fastq(path_fastq):
        for iterator_k_mer in cut_kmer(seq_read, k_mer):
            if iterator_k_mer in dict_k_mer:
                dict_k_mer[iterator_k_mer] +=1
            else:
                dict_k_mer[iterator_k_mer] = 1
    return(dict_k_mer)

def build_graph(dict_k_mer):
    """
    creat a networkx graph
    """
    G=nx.DiGraph()
    for k_mer in dict_k_mer:
        G.add_node((k_mer[:-1],k_mer[1:]), weight=dict_k_mer[k_mer])
    print(G.nodes)
  
def get_starting_nodes():
    pass


def std():
    pass


def get_sink_nodes():
    pass


def path_average_weight():
    pass


def remove_paths():
    pass


def select_best_path():
    pass


def save_contigs():
    pass


def get_contigs():
    pass


def solve_bubble():
    pass


def simplify_bubbles():
    pass


def solve_entry_tips():
    pass


def solve_out_tips():
    pass
    
def main():
    """
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", nargs="?",\
        help="fichier fastq single end")
    parser.add_argument("-k", nargs="?",\
        help="taille des kmer (optionnel - default 21)",)
    parser.add_argument("-o", nargs="?",\
        help="fichier config")

    args = parser.parse_args()
    
    dict_k_mer = build_kmer_dict(args.i, 3)
    build_graph(dict_k_mer)
if __name__ == "__main__":
    main()
