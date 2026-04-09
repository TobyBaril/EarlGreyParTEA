"""cluster_utils.py — shared cd-hit-est .clstr parsing for EarlGreyParTEA.

This module is imported by saturation_plot.py, shared_unique_plot.py, and
busco_te_qc_plot.py.  Centralising the parser avoids duplication and ensures
all scripts apply identical cluster-attribution logic.
"""


def parse_clstr(clstr_file, species_list):
    """Parse a cd-hit-est .clstr file and attribute clusters to source genomes.

    Each sequence in the combined library was prefixed with its source before
    clustering:
        {species}_         for de-novo TE families
        REPMASKER_{spec}_  for RepeatMasker database families
        CUSTOM_            for user-supplied custom library sequences

    Species names are matched longest-first so that a name that is a prefix of
    another (e.g. 'Sp1' vs 'Sp10') is not misattributed.

    Parameters
    ----------
    clstr_file : str
        Path to the .clstr file produced by cd-hit-est.
    species_list : list[str]
        Species names exactly as used for sequence prefixing.

    Returns
    -------
    species_to_clusters : dict[str, set[int]]
        Maps each species name to the set of cluster IDs it contributed at
        least one sequence to.
    existing_clusters : set[int]
        Cluster IDs that contain at least one REPMASKER_ or CUSTOM_ sequence.
        These form the x=0 baseline in the saturation curve.
    cluster_to_species : dict[int, frozenset[str]]
        Maps each cluster ID to the frozenset of species that contributed to
        it.  Used by shared_unique_plot.py.
    rep_name_to_cluster : dict[str, int]
        Maps the *representative* sequence NAME (the entry flagged with '*' in
        the .clstr file, stripped of any trailing '#class/subclass') to its
        cluster ID.  Used by shared_unique_plot.py when matching GFF NAME=
        attributes back to clusters.
    """
    sorted_species = sorted(species_list, key=len, reverse=True)
    species_to_clusters = {sp: set() for sp in species_list}
    existing_clusters = set()
    cluster_members = {}      # cluster_id -> list of (header, is_rep)
    current_cluster = None

    with open(clstr_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">Cluster"):
                current_cluster = int(line.split()[1])
                cluster_members[current_cluster] = []
            elif current_cluster is not None and ">" in line:
                # Extract header: text between ">" and "..."
                header = line.split(">", 1)[1].split("...")[0]
                is_rep = line.strip().endswith("*")
                cluster_members[current_cluster].append((header, is_rep))

    # Build all derived dicts in a single pass over parsed members
    rep_name_to_cluster = {}
    cluster_to_species_raw = {}

    for cid, members in cluster_members.items():
        species_in_cluster = set()
        for header, is_rep in members:
            if header.startswith("REPMASKER_") or header.startswith("CUSTOM_"):
                existing_clusters.add(cid)
                if is_rep:
                    # Strip class suffix if present (NAME#class/subclass)
                    rep_bare = header.split("#")[0]
                    rep_name_to_cluster[rep_bare] = cid
                    rep_name_to_cluster[header] = cid
            else:
                for sp in sorted_species:
                    if header.startswith(sp + "_"):
                        species_in_cluster.add(sp)
                        species_to_clusters[sp].add(cid)
                        if is_rep:
                            # The NAME= in the GFF will be the header without
                            # the {species}_ prefix, but in cluster_all_species
                            # we keep the full prefixed name.  Store both the
                            # full name and the bare name (after stripping
                            # {species}_ prefix) so GFF matching succeeds.
                            rep_bare = header.split("#")[0]
                            rep_name_to_cluster[rep_bare] = cid
                            rep_name_to_cluster[header] = cid
                            # Also store with stripped species prefix so that
                            # GFF NAME= values (which do NOT carry the prefix)
                            # can be resolved.
                            bare_no_prefix = header[len(sp) + 1:].split("#")[0]
                            rep_name_to_cluster[bare_no_prefix] = cid
                        break
        cluster_to_species_raw[cid] = species_in_cluster

    cluster_to_species = {cid: frozenset(sp) for cid, sp in cluster_to_species_raw.items()}

    return species_to_clusters, existing_clusters, cluster_to_species, rep_name_to_cluster
