from collections import defaultdict, Counter

def get_seg_and_name(L):
    (_, _, seg, name) = L.split("|")
    return (seg, name)

def process_cluster(cluster_name, protein_list, seg_list, seg_miss, protein_miss):
    # Handle proteins
    for (C_Segment, C_Protein), counter in protein_list.items():
        if len(counter) > 0:
           seg_miss.write(f"{cluster_name}\t{dict(counter)}\n")

    # Handle segments
    for (C_Segment, C_Protein), counter in seg_list.items():
        for segment, count in counter.items():
            protein_miss.write(f"{cluster_name}\t{dict(counter)}\n")


seg_miss = open("analysis/seg_miss.txt", "w")
protein_miss = open("analysis/protein_miss.txt", "w")

Current_cluster = None
protein_list = defaultdict(Counter)
seg_list = defaultdict(Counter)

with open("analysis/clusters.tsv") as f:
    for line in f:
        cluster, member = line.strip().split("\t")

        # If cluster changes, process previous one
        if Current_cluster is not None and Current_cluster != cluster:
            process_cluster(Current_cluster, protein_list, seg_list,
                            seg_miss, protein_miss)

            # Reset for next cluster
            protein_list = defaultdict(Counter)
            seg_list = defaultdict(Counter)

        Current_cluster = cluster

        # Extract info
        S_clust, N_clust = get_seg_and_name(cluster)
        S_member, N_member = get_seg_and_name(member)

        # Count occurrences
        protein_list[(S_clust, N_clust)][f"{S_member}:{N_member}"] += 1
        seg_list[(S_clust, N_clust)][f"{S_member}:{N_member}"] += 1

if Current_cluster is not None:
    process_cluster(Current_cluster, protein_list, seg_list,
                    seg_miss, protein_miss)

seg_miss.close()
protein_miss.close()
