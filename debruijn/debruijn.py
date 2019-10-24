import os
import argparse
import statistics
import random
import networkx as nx

def read_fastq(path_fastq):
    """
    read fastq file and output a list of the reads
    """
    with open(path_fastq,'r') as filin:
        for ligne in filin:
            yield next(filin)[0:-1]
            next(filin)
            next(filin)

def cut_kmer(seq, k_mer):
    """
    take a sequence and a k-mer size and return an iterator of k-mer
    """
    for i in range(len(seq)-k_mer+1):
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
        G.add_edge(k_mer[:-1],k_mer[1:], weight=dict_k_mer[k_mer])
    return(G)
  
def get_starting_nodes(G):
    start_nodes = []
    for k_mer in G:
        if list(G.predecessors(k_mer)) == []:
            start_nodes.append(k_mer)
    return start_nodes

def std(list_value):
    """
    compute the standard deviation
    """
    return statistics.stdev(list_value)


def get_sink_nodes(G):
    end_nodes = []
    for k_mer in G:
        if list(G.successors(k_mer)) == []:
            end_nodes.append(k_mer)
    return end_nodes


def path_average_weight(G, chemin):
    tot_weight = []
    list_chemin = list(chemin)
    print(list_chemin)
    for i in range(0,len(list_chemin)-1):
        node_1 = list_chemin[i]
        node_2 = list_chemin[i+1]
        tot_weight.append(G[node_1][node_2]["weight"])
    return sum(tot_weight)/len(tot_weight)


def remove_paths(graph, list_chemin, delete_entry_node, delete_sink_node):
    list_chemin = list(list_chemin)
    list_path_clean = []       
    start_nodes = list(get_starting_nodes(graph))
    end_nodes = list(get_sink_nodes(graph))
    for chemin in list_chemin:
        for node in chemin:
            if delete_entry_node == True and delete_sink_node == False:
                print(1, node, start_nodes)
                if node in start_nodes:
                    graph.remove_node(node)
                    print(node,1)
            if delete_entry_node == False and delete_sink_node == True:
                if node in end_nodes:
                    graph.remove_node(node)
            if delete_entry_node == False and delete_sink_node == False:
                if node not in start_nodes and node not in end_nodes:
                    graph.remove_node(node)
            if delete_entry_node == True and delete_sink_node == True:
                graph.remove_node(node)
    return graph


def select_best_path(graph, list_chemin, list_size, list_weight, delete_entry_node = False, delete_sink_node = False):
    #frequent
    
    max_weight = []
    index_weight = []
    i=0
    for weight in list_weight:
        if max_weight == []:
            max_weight.append(weight)
            index_weight.append(i)
        else:
            print(weight, "--")
            if weight > max_weight[0]:
                max_weight = []
                index_weight = [] 
                max_weight.append(weight)
                index_weight.append(i)
            elif weight == max_weight[0]:
                max_weight.append(weight)
                index_weight.append(i)
        i +=1
    best_chemin = list_chemin[index_weight[0]]
    max_size = []
    index_size = []   
    if len(max_weight) > 1:

        for index in index_weight:
            if max_size == []:
                max_size.append(list_size[index])
            else:
                if list_size[index] > max_size[0]:
                    max_size = []
                    index_size = [] 
                    max_size.append(list_size[index])
                    index_size.append(index)
                elif list_size[index] == max_size[0]:
                    max_size.append(list_size[index])
                    index_size.append(index)
        best_chemin = list_chemin[index_size[0]]
        print("sup à 1")
    if len(max_size) > 1:
        best_chemin = list_chemin[index_size[randint(0,len(max_size))]]
    for chemin in list_chemin:
        if chemin not in best_chemin:
            remove_paths(graph, chemin, delete_entry_node, delete_sink_node)
    print(max_weight, index_weight)
    return graph

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))
    
def save_contigs(contigs, path_output):
    file_output = open(path_output,"w")
    n_contig = 0
    for contig in contigs:
        file_output.write(">contig_"+str(n_contig)+" len="+str(contig[1])+"\n")
        file_output.write(fill(contig[0]))
        file_output.write("\n")
        n_contig += 1
    file_output.close()
    pass


def get_contigs(G, start_nodes, end_nodes):
    contigs = []
    for start in start_nodes:
        for end in end_nodes:
            if list(nx.all_simple_paths(G, start, end)) != []:
                contig_seq = nx.shortest_path(G, start, end)
                contig = []
                for i in range(len(contig_seq)-1):
                    contig.append(contig_seq[i][0])
                    flag = i
                contig.append(contig_seq[flag+1])
                join_contig = ''.join(contig)
                size_contig = len(join_contig)
                contigs.append((join_contig,size_contig))
    return contigs


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
    
    random.seed(9001)
    
    dict_k_mer = build_kmer_dict(args.i, int(args.k))
    G = build_graph(dict_k_mer)
    sart_nodes = get_starting_nodes(G)
    end_nodes = get_sink_nodes(G)
    contigs = get_contigs(G, sart_nodes, end_nodes)
    save_contigs(contigs, "test.txt")
    #path_average_weight(G,)
    
    graph_1 = nx.DiGraph()
    graph_1.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    graph_1 = select_best_path(graph_1, [[1,2], [3,2]], [1, 1], [5, 10], delete_entry_node=True)



if __name__ == "__main__":
    main()
