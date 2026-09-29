from collections import defaultdict, Counter
import sys

def get_seg_and_name(L):
    (_, _, seg, name) = L.split("|")
    return (seg, name)

def process_cluster(cluster_name, protein_list, clusters):
    # Handle proteins
    for (C_Segment, C_Protein), counter in protein_list.items():
        if len(counter) > 0:
           clusters.write(f"{cluster_name}\t{dict(counter)}\n")


clusters = open(sys.argv[1], "w")

Current_cluster = None
protein_list = defaultdict(Counter)

with open(sys.argv[2]) as f:
    for line in f:
        cluster, member = line.strip().split("\t")

        # If cluster changes, process previous one
        if Current_cluster is not None and Current_cluster != cluster:
            process_cluster(Current_cluster, protein_list,clusters)
            protein_list = defaultdict(Counter)

        Current_cluster = cluster

        # Extract info
        S_clust, N_clust = get_seg_and_name(cluster)
        S_member, N_member = get_seg_and_name(member)

        # Count occurrences
        protein_list[(S_clust, N_clust)][f"{S_member}:{N_member}"] += 1

if Current_cluster is not None:
    process_cluster(Current_cluster, protein_list,clusters)
                    

clusters.close()
