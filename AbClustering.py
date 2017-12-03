import sys
import editdistance

def define_clones(clonal_ab_library, sim_cutoff):
    from IgLibrary import Clone
    # group by V, J
    # then subdivide into clones based on CDR3
    # then further subdivide into clones based on full sequence

    vj_groups = {}
    for ab in clonal_ab_library.entries:
        top_v, top_d, top_j = ab.extract_vdj()
        if (top_v, top_j) in vj_groups:
            vj_groups[(top_v,top_j)].append(ab)
        else:
            vj_groups[(top_v,top_j)] = [ab]

    ab_clones = []
    for vj_call, ab_group in vj_groups.items():
        # cluster based on cdr3
        vj_clusters = cluster_seqs(ab_group, sim_cutoff, 'cdr3_seq')

        vj_final_clusters = []
        # refine clusters based on full sequence
        for vj_clust in vj_clusters:
            vj_refined_cluster = cluster_seqs(vj_clust, sim_cutoff, 'ungapped_seq')
            vj_final_clusters += vj_refined_cluster

        for i, clust in enumerate(vj_final_clusters):
            clone_id = "_".join([vj_call[0], vj_call[1], "clone", str(i+1)])
            ab_clones.append(Clone(clone_id, clust))

    ab_clones.sort(key=lambda x: x.clone_size, reverse=True)

    return ab_clones


def cluster_seqs(ab_sequences, sim_cutoff, clust_field='cdr3_nt'):  # or seq
    clusters = []
    for s in ab_sequences:
        cluster_id = calc_similarity(getattr(s, clust_field), clusters, sim_cutoff, clust_field)
        if cluster_id is not None:
            clusters[cluster_id].append(s)
        else:
            clusters.append([s])

    return clusters


def calc_similarity(sequence, clusters, sim_cutoff, clust_field):
    for c, c_seqs in enumerate(clusters):
        dist = editdistance.eval(sequence, getattr(c_seqs[0], clust_field))
        similarity = dist/len(sequence)
        if similarity <= sim_cutoff:
            return c

    return None
