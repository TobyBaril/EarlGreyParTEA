import os
import sys
import subprocess
import urllib.request
from pathlib import Path
import yaml

def running_tea(stage="Starting Earl Grey"):
    """Display ASCII art tea cup with stage name"""
    tea_art = rf"""    
          )  (
         (   ) )
         ) ( (
       _______)_
    .-'---------|  
   ( C|/\/\/\/\/|
    '-./\/\/\/\/|
      '_________'
       '-------'
    <<< {stage} >>>"""
    print(tea_art)

def convert_seconds(seconds):
    """Convert seconds to HH:MM:SS.ss format"""
    seconds = int(seconds)
    h = seconds // 3600
    m = seconds % 3600 // 60
    s = seconds % 60
    return f"{h:02d}:{m:02d}:{s:05.2f}"


# --------------------------------------------------
# Message helpers (consistent style)
# --------------------------------------------------
def msg_header(title):
    print("\n" + "=" * 60)
    print(f"{title}")
    print("=" * 60)

def msg_info(text):
    print(f"[INFO] {text}")

def msg_warn(text):
    print(f"[WARNING] {text}")

def msg_error(text):
    sys.exit(f"\n[ERROR] {text}\n")

# --------------------------------------------------
# Pipeline mode visualization
# --------------------------------------------------
def show_pipeline_mode_visualization(pipeline_mode, config):
    """Display visual table showing which rules are active in current pipeline mode
    
    Args:
        pipeline_mode: One of 'full', 'libconstruct', 'annotate'
        config: Full config dict to resolve conditional rules
    """
    
    # Extract config parameters that affect rule execution
    run_heliano = config.get('run_heliano', False)
    softmask = config.get('softmask', False)
    repeatmasker_species = config.get('repeatmasker_species', '')
    custom_library = config.get('custom_library', '')
    skip_clustering = config.get('skip_clustering', False)
    run_shared_unique = config.get('run_shared_unique', False)
    run_busco_phylo = config.get('run_busco_phylo', False)
    
    # Define rule groups and their activity per mode
    # ✓ = active, ✗ = inactive, ○ = conditional (depends on config)
    rules = {
        'Library Construction': {
            'prep_genome': {'libconstruct': '✓', 'annotate': '✓', 'full': '✓'},
            'repeatmasker_initial': {'libconstruct': '○', 'annotate': '✗', 'full': '○'},
            'build_db': {'libconstruct': '✓', 'annotate': '✗', 'full': '✓'},
            'repeatmodeler': {'libconstruct': '✓', 'annotate': '✗', 'full': '✓'},
            'testrainer': {'libconstruct': '✓', 'annotate': '✗', 'full': '✓'},
            'clustering': {'libconstruct': '✓', 'annotate': '✗', 'full': '✓'},
            'split_chimeras': {'libconstruct': '○', 'annotate': '✗', 'full': '○'},
            'saturation_plot': {'libconstruct': '✓', 'annotate': '✗', 'full': '✓'},
        },
        'Annotation': {
            'symlink_library': {'libconstruct': '✗', 'annotate': '✓', 'full': '✗'},
            'repeatmasker_annotation': {'libconstruct': '✗', 'annotate': '✓', 'full': '✓'},
            'heliano_detection': {'libconstruct': '✗', 'annotate': '○', 'full': '○'},
            'merge_repeats': {'libconstruct': '✗', 'annotate': '✓', 'full': '✓'},
            'calculate_divergence': {'libconstruct': '✗', 'annotate': '✓', 'full': '✓'},
            'generate_charts': {'libconstruct': '✗', 'annotate': '✓', 'full': '✓'},
            'softmasked_genome': {'libconstruct': '✗', 'annotate': '○', 'full': '○'},
        },
        'Optional: Shared/Unique TE Content': {
            'shared_unique_plot (cluster)':        {'libconstruct': '✗', 'annotate': '✗', 'full': '○'},
            'shared_unique_pa_plot (presence/abs)':{'libconstruct': '✗', 'annotate': '○', 'full': '✗'},
            'shared_unique_*_phylo':               {'libconstruct': '✗', 'annotate': '○', 'full': '○'},
        },
        'Optional: BUSCO Phylogenomics': {
            'run_busco':               {'libconstruct': '○', 'annotate': '○', 'full': '○'},
            'busco_summary_table':     {'libconstruct': '○', 'annotate': '○', 'full': '○'},
            'extract_busco_aa':        {'libconstruct': '○', 'annotate': '○', 'full': '○'},
            'align_busco_gene':        {'libconstruct': '○', 'annotate': '○', 'full': '○'},
            'create_supermatrix':      {'libconstruct': '○', 'annotate': '○', 'full': '○'},
            'run_fasttree':            {'libconstruct': '○', 'annotate': '○', 'full': '○'},
            'busco_completeness_phylo':{'libconstruct': '○', 'annotate': '○', 'full': '○'},
            'busco_te_qc':             {'libconstruct': '○', 'annotate': '✗', 'full': '○'},
        },
    }
    
    print("\n" + "="*60)
    print(f"PIPELINE MODE: {pipeline_mode.upper()}")
    print("="*60)
    
    for group, group_rules in rules.items():
        print(f"\n{group}:")
        for rule_name, activity in group_rules.items():
            status = activity[pipeline_mode]
            # Update conditional rules based on actual config
            if status == '○':
                if rule_name == 'heliano_detection':
                    status = '✓' if run_heliano else '✗'
                elif rule_name == 'softmasked_genome':
                    status = '✓' if softmask else '✗'
                elif rule_name == 'repeatmasker_initial':
                    status = '✓' if (repeatmasker_species or custom_library) else '✗'
                elif rule_name == 'split_chimeras':
                    status = '✓' if (config.get('split_chimeras', False) and not skip_clustering) else '✗'
                elif rule_name in (
                    'shared_unique_plot (cluster)',
                    'shared_unique_pa_plot (presence/abs)',
                    'shared_unique_*_phylo',
                ):
                    status = '✓' if run_shared_unique else '✗'
                elif rule_name in (
                    'run_busco', 'busco_summary_table', 'extract_busco_aa',
                    'align_busco_gene', 'create_supermatrix', 'run_fasttree',
                    'busco_completeness_phylo', 'busco_te_qc',
                ):
                    status = '✓' if run_busco_phylo else '✗'
            
            # Add note for clustering behavior
            note = ""
            if rule_name == 'clustering' and status == '✓' and skip_clustering:
                note = " (concatenate mode)"
            # Add note for presence/absence shared_unique in annotate mode
            if rule_name == 'shared_unique_*_phylo' and status == '✓':
                note = " (requires run_busco_phylo: true)"
            
            print(f"  {status} {rule_name}{note}")
    
    print("\n" + "="*60)
    print("Legend: ✓ = active  ✗ = inactive  ○ = enabled via config")
    print("="*60 + "\n")

# --------------------------------------------------
# Validation + defaults with custom messages
# --------------------------------------------------
def validate_parameters(config, outfile = None):

    msg_header("Earl Grey configuration check")

    # ---- Required parameters ----
    required = ['genome', 'species', 'output_dir']
    for param in required:
        if not config.get(param):
            msg_error(f"Required parameter '{param}' not specified in config file")

    msg_info("Required parameters detected")

    # ---- Defaults + messages ----
    defaults = {
        'iterations': (10, "De Novo Sequences will be extended through a maximum of {} iterations"),
        'max_consensus_seqs': (20, "{} sequences will be used in BEAT consensus generation"),
        'skip_clustering': (False, None),
        'clustering_identity': (0.8, None),
        'clustering_coverage': (0.8, None),
        'clustering_coverage_long': (0.0, None),
        'clustering_length_diff': (0.5, None),
        'split_chimeras': (False, None),
        'chimera_overlap_min': (50, None),
        'chimera_min_members': (3, None),
        'chimera_min_component_span': (0.1, None),
        'softmask': (False, None),
        'margin': (False, None),
        'flank': (1000, "Blast, extend, align, trim process will add {}bp to each end in each iteration"),
        'min_consensus_seqs': (3, "Blast, extend, align, trim process will require {} sequences to generate a new consensus sequence"),
        'run_heliano': (False, None),
        'repeatmasker_species': ("", None),
        'custom_library': ("", None),
        'saturation_permutations': (100, None),
        'run_shared_unique': (False, None),
        'run_busco_phylo': (False, None),
        'busco_lineage': ("", None),
        'busco_prefix': ("busco", None),
        'busco_min_occupancy': (0.5, None),
        'slurm_partition': ("", None),
        'slurm_account': ("", None),
        'slurm_extra': ("", None),
        'lsf_queue': ("", None),
        'lsf_project': ("", None),
        'lsf_extra': ("", None),
    }

    msg_header("Parameter values")

    for param, (default_val, message_template) in defaults.items():
        if param not in config or config.get(param) is None or config.get(param) == "":
            config[param] = default_val
            if message_template:
                msg_info(message_template.format(default_val))
        else:
            if message_template:
                msg_info(message_template.format(config[param]))

    # --------------------------------------------------
    # Verbose user messages 
    # --------------------------------------------------
    msg_header("Pipeline behaviour")

    # RepeatMasker
    repspec = config.get('repeatmasker_species', "")
    custom_lib = config.get('custom_library', "")
    
    if repspec and custom_lib:
        msg_error("Both RepeatMasker species and custom library specified - only one can be used at a time")
    elif repspec:
        msg_info(f"Running with initial RepeatMasker masking using species: {repspec}")
    elif custom_lib:
        msg_info(f"Running with initial RepeatMasker masking using custom library: {os.path.basename(custom_lib)}")
    else:
        msg_info("RepeatMasker species/library not specified, running Earl Grey without an initial mask with known repeats")

    # Clustering
    skip_clustering = config.get('skip_clustering', False)
    if skip_clustering:
        msg_info("TE consensus sequences will NOT be clustered (libraries will be concatenated)")
    else:
        cluster_id = config.get('clustering_identity', 0.8)
        cluster_cov = config.get('clustering_coverage', 0.8)
        cluster_cov_long = config.get('clustering_coverage_long', 0.0)
        cluster_len = config.get('clustering_length_diff', 0.5)
        if not (0.0 < cluster_id <= 1.0):
            msg_error(f"clustering_identity={cluster_id} is out of range. Must be between 0 (exclusive) and 1 (inclusive).")
        if not (0.0 <= cluster_cov <= 1.0):
            msg_error(f"clustering_coverage={cluster_cov} is out of range. Must be between 0.0 and 1.0.")
        if not (0.0 <= cluster_cov_long <= 1.0):
            msg_error(f"clustering_coverage_long={cluster_cov_long} is out of range. Must be between 0.0 and 1.0.")
        if not (0.0 <= cluster_len <= 1.0):
            msg_error(f"clustering_length_diff={cluster_len} is out of range. Must be between 0.0 and 1.0.")
        aL_note = f", aL: {cluster_cov_long}" if cluster_cov_long > 0.0 else " (aL: disabled)"
        msg_info(f"TE consensus sequences will be clustered (identity: {cluster_id}, aS: {cluster_cov}, length_diff: {cluster_len}{aL_note})")
        if cluster_cov_long > 0.0:
            msg_info(f"Long-sequence coverage filter active: sequences must contribute >= {cluster_cov_long:.0%} of the longer sequence to the alignment")
        msg_warn("Clustering may affect subfamilies and create chimeras")
        # Chimera splitting
        if config.get('split_chimeras', False):
            ovlp = config.get('chimera_overlap_min', 50)
            min_mem = config.get('chimera_min_members', 3)
            span = config.get('chimera_min_component_span', 0.1)
            if not (ovlp > 0):
                msg_error(f"chimera_overlap_min={ovlp} must be a positive integer.")
            if not (min_mem >= 1):
                msg_error(f"chimera_min_members={min_mem} must be >= 1.")
            if not (0.0 < span <= 1.0):
                msg_error(f"chimera_min_component_span={span} is out of range. Must be between 0 (exclusive) and 1 (inclusive).")
            msg_info(f"Chimera detection enabled (overlap_min: {ovlp} nt, min_members: {min_mem}, min_component_span: {span})")
        else:
            msg_info("Chimera detection disabled (split_chimeras: false)")

    # Saturation plot
    pipeline_mode_for_sat = config.get('pipeline_mode', 'full')
    if pipeline_mode_for_sat in ('full', 'libconstruct'):
        sat_perms = config.get('saturation_permutations', 100)
        if skip_clustering:
            msg_info(f"TE family saturation plot will be generated ({sat_perms} permutations, fallback mode: raw sequence counts)")
            msg_warn("Saturation plot: cross-genome deduplication unavailable when skip_clustering=True")
        else:
            msg_info(f"TE family saturation plot will be generated ({sat_perms} permutations)")

    # SoftMask
    softmask = config.get('softmask', False)
    if softmask is True or softmask == 'yes':
        msg_info("Softmasked genome will be generated")
    else:
        msg_info("Softmasked genome will not be generated")

    # Margin
    margin = config.get('margin', False)
    if margin is True or margin == 'yes':
        msg_info("Short TE sequences (<100bp) will be removed")
    else:
        msg_info("Short TE sequences (<100bp) will not be removed")

    # Helitrons
    run_heliano = config.get('run_heliano', False)
    if run_heliano is True or run_heliano == 'yes':
        msg_info("HELITRON detection will be run using HELIANO")
    else:
        msg_info("HELITRON detection will not be run")

    # Pipeline mode
    pipeline_mode = config.get('pipeline_mode', 'full')
    if pipeline_mode not in ['full', 'libconstruct', 'annotate']:
        msg_error(f"Invalid pipeline_mode '{pipeline_mode}'. Must be 'full', 'libconstruct', or 'annotate'")
    
    if pipeline_mode == 'full':
        msg_info("Pipeline mode: FULL - Complete library construction and annotation")
    elif pipeline_mode == 'libconstruct':
        msg_info("Pipeline mode: LIBCONSTRUCT - Library construction only (stops after clustering)")
    elif pipeline_mode == 'annotate':
        msg_info("Pipeline mode: ANNOTATE - Annotation only with user-supplied library")
        annotation_lib = config.get('annotation_library', '')
        if not annotation_lib:
            msg_error("Pipeline mode 'annotate' requires 'annotation_library' parameter specifying path to library fasta")
        if not Path(annotation_lib).exists():
            msg_error(f"Annotation library not found: {annotation_lib}")
        msg_info(f"Using annotation library: {os.path.basename(annotation_lib)}")

    # Optional module: shared/unique TE content
    run_shared_unique = config.get('run_shared_unique', False)
    if run_shared_unique:
        if pipeline_mode == 'libconstruct':
            msg_error(
                "run_shared_unique=true requires annotation outputs (GFF files) "
                "which are not produced in 'libconstruct' mode. "
                "Use 'full' or 'annotate' mode instead."
            )
        elif pipeline_mode == 'annotate':
            msg_warn(
                "run_shared_unique=true with pipeline_mode='annotate': "
                "presence/absence strategy will be used (no .clstr file). "
                "Homologous families annotated under different names will "
                "appear as unique in each species."
            )
            msg_info("Shared/unique TE content analysis: ENABLED (presence/absence mode)")
        else:
            msg_info("Shared/unique TE content analysis: ENABLED (cluster-based mode)")
    else:
        msg_info("Shared/unique TE content analysis: disabled (set run_shared_unique: true to enable)")

    # Optional module: BUSCO phylogenomics
    run_busco_phylo = config.get('run_busco_phylo', False)
    if run_busco_phylo:
        busco_lineage = config.get('busco_lineage', '')
        if not busco_lineage:
            msg_error(
                "run_busco_phylo=true requires 'busco_lineage' to be set "
                "(e.g. busco_lineage: fungi_odb10). "
                "See https://busco.ezlab.org/busco_userguide.html for available lineages."
            )
        min_occ = config.get('busco_min_occupancy', 0.5)
        if not (0.0 < min_occ <= 1.0):
            msg_error(
                f"busco_min_occupancy={min_occ} is out of range. "
                "Must be a fraction between 0 (exclusive) and 1 (inclusive)."
            )
        msg_info(f"BUSCO phylogenomics: ENABLED (lineage={busco_lineage}, "
                 f"min_occupancy={min_occ:.0%})")
    else:
        msg_info("BUSCO phylogenomics: disabled (set run_busco_phylo: true to enable)")

    # Show visual pipeline mode overview
    show_pipeline_mode_visualization(pipeline_mode, config)

    # --------------------------------------------------
    # Structural validation 
    # --------------------------------------------------
    msg_header("Input validation")

    if not isinstance(config['genome'], dict):
        msg_error("'genome' must be a dictionary: species → fasta")

    for sp, path in config['genome'].items():
        if not Path(path).exists():
            msg_error(f"Genome file for species '{sp}' not found: {path}")

    if not isinstance(config['species'], list):
        msg_error("'species' must be a list")

    for sp in config['species']:
        if sp not in config['genome']:
            msg_error(f"Species '{sp}' listed but no genome provided")

    # Auto-detect script_dir if not specified (for conda/mamba installations)
    if not config.get("script_dir") or config["script_dir"] == "":
        # Try to find script directory in standard conda installation locations
        possible_locations = []
        
        # Check CONDA_PREFIX environment variable (version-agnostic)
        if os.environ.get("CONDA_PREFIX"):
            conda_prefix = os.environ["CONDA_PREFIX"]
            # Look for any version of earlgrey in share directory
            share_dir = os.path.join(conda_prefix, "share")
            if os.path.isdir(share_dir):
                # Find all earlgrey-* directories
                import glob
                for earlgrey_dir in glob.glob(os.path.join(share_dir, "earlgrey-*")):
                    scripts_path = os.path.join(earlgrey_dir, "scripts")
                    if os.path.isdir(scripts_path):
                        possible_locations.append(scripts_path)
                # Also check generic earlgrey directory (no version)
                generic_path = os.path.join(share_dir, "earlgrey", "scripts")
                if os.path.isdir(generic_path):
                    possible_locations.append(generic_path)
        
        # Check relative to this script's location
        script_location = Path(__file__).parent.parent / "scripts"
        possible_locations.append(str(script_location))
        
        # Try each location, prioritizing newer versions
        for location in possible_locations:
            if os.path.isdir(location):
                testrainer_check = os.path.join(location, "TEstrainer")
                if os.path.isdir(testrainer_check):
                    config["script_dir"] = location
                    msg_info(f"Auto-detected script directory: {location}")
                    break
        
        if not config.get("script_dir") or config["script_dir"] == "":
            msg_error("Script directory not specified and could not be auto-detected. Please set 'script_dir' in your config file or ensure EarlGrey is properly installed.")
    
    if not os.path.isdir(config["script_dir"]):
        msg_error(f"Script directory not found: {config['script_dir']}")
    
    testrainer_dir = os.path.join(config["script_dir"], "TEstrainer")
    if not os.path.isdir(testrainer_dir):
        msg_error(f"TEstrainer module not found in {config['script_dir']}. Please ensure EarlGrey is properly installed.")

    # --------------------------------------------------
    # EarlGrey configuration check
    # --------------------------------------------------
    msg_header("EarlGrey configuration check")

    try:
        import glob as _glob
        # Locate RepeatMasker and derive conda prefix
        rm_result = subprocess.run(["which", "RepeatMasker"],
                                   capture_output=True, text=True, check=True)
        rm_path = rm_result.stdout.strip()
        conda_prefix_rm = rm_path.replace("/bin/RepeatMasker", "")

        # Locate famdb directory — supports both:
        #   FamDB 3.0.0+ standalone (share/famdb-*/Libraries/famdb/) — Dfam 4.0+
        #   FamDB < 3.0.0 embedded  (share/RepeatMasker/Libraries/famdb/) — Dfam 3.9
        famdb_shares = sorted(_glob.glob(
            os.path.join(conda_prefix_rm, "share", "famdb-*")
        ))
        if famdb_shares:
            library_path = os.path.join(famdb_shares[-1], "Libraries", "famdb")
        else:
            library_path = os.path.join(
                conda_prefix_rm, "share", "RepeatMasker", "Libraries", "famdb"
            )

        if not os.path.isdir(library_path):
            msg_error(
                f"Dfam famdb directory not found: {library_path}\n"
                "Please ensure EarlGrey (>=7.3.0) is properly installed and "
                "Dfam partitions have been downloaded via: download_dfam.py"
            )

        # Fast-path: legacy completion marker (written by configure_dfam.sh in <=7.2 envs)
        complete_marker = os.path.join(library_path, ".earlgrey.config.complete")
        if os.path.exists(complete_marker):
            msg_info("EarlGrey configuration marker found — Dfam library configured")
        else:
            # Check for h5 partition files:
            #   Dfam 4.0 names: dfam40.*.h5
            #   Dfam 3.9 names: dfam39_full.*.h5
            h5_files = [
                f for f in os.listdir(library_path)
                if f.endswith(".h5") and ("dfam40" in f or "dfam39_full" in f)
            ]
            if len(h5_files) < 1:
                _print_dfam_setup_instructions(library_path, famdb_shares)
                msg_error(
                    "\nDfam library partitions not found in:\n"
                    f"  {library_path}\n"
                    "Download them using download_dfam.py (see instructions above)."
                )
            else:
                msg_info(f"Found {len(h5_files)} Dfam partition file(s) in {library_path}")

    except subprocess.CalledProcessError:
        msg_warn("Could not locate RepeatMasker installation — skipping library check")
    except Exception as e:
        msg_warn(f"Error checking EarlGrey configuration: {e}")

    # --------------------------------------------------
    # Absolutize all paths so rules work correctly after `cd`
    # --------------------------------------------------
    config['output_dir'] = os.path.abspath(config['output_dir'])
    config['genome'] = {sp: os.path.abspath(p) for sp, p in config['genome'].items()}
    if config.get('custom_library'):
        config['custom_library'] = os.path.abspath(config['custom_library'])
    if config.get('annotation_library'):
        config['annotation_library'] = os.path.abspath(config['annotation_library'])

    # --------------------------------------------------
    # Output setup
    # --------------------------------------------------
    msg_header("Output setup")

    outdir = Path(config['output_dir'])
    outdir.mkdir(parents=True, exist_ok=True)
    msg_info(f"Output directory: {outdir}")

    for sp in config['species']:
        sp_dir = outdir / f"{sp}_EarlGrey"
        sp_dir.mkdir(parents=True, exist_ok=True)
        msg_info(f"Created directory: {sp_dir}")

    # Save validated config (optional but recommended)
    if outfile:
        with open(outfile, "w") as f:
            yaml.safe_dump(config, f, sort_keys=False)
            msg_info(f"Validated configuration saved to: {outfile}")


    msg_header("Configuration complete")
    msg_info("All checks passed. Workflow ready to start.\n")
    
    print("\nPlease cite the following paper when using this software:")
    print("Baril, T., Galbraith, J. and Hayward, A., 2024. Earl Grey: a fully automated user-friendly transposable element annotation and analysis pipeline. Molecular Biology and Evolution, 41(4), p.msae068. \n")
    

    return config


def make_directories(directory, species, RepSpec=None, startCust=None, run_heliano=None):
    outdir = os.path.join(directory, f"{species}_EarlGrey")
    os.makedirs(outdir, exist_ok=True)

    if RepSpec or startCust:
        os.makedirs(os.path.join(outdir, f"{species}_RepeatMasker"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_Database"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_RepeatModeler"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_strainer"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_RepeatMasker_Against_Custom_Library"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_RepeatLandscape"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_mergedRepeats"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_summaryFiles"), exist_ok=True)
    if run_heliano:
        os.makedirs(os.path.join(outdir, f"{species}_heliano"), exist_ok=True)
    return outdir

def _print_dfam_setup_instructions(library_path, famdb_shares):
    """Print Dfam setup instructions appropriate for the detected environment."""
    print("\n" + "=" * 60)
    if famdb_shares:
        # New FamDB 3.0.0+ / Dfam 4.0 — standalone conda package
        print("Dfam 4.0 library setup required")
        print("=" * 60)
        print(
            "EarlGrey >=7.3.0 uses FamDB 3.0.0 / Dfam 4.0.\n"
            "RepeatMasker is pre-configured by the conda package — "
            "no 'perl ./configure' is needed.\n\n"
            "Download the Dfam 4.0 partitions using the interactive tool:\n\n"
            "    download_dfam.py\n\n"
            "This will guide you through selecting which partitions to download\n"
            "(curated consensus, HMMs, etc.) and write them to:\n"
            f"    {library_path}\n\n"
            "After download completes, re-run the pipeline."
        )
    else:
        # Legacy FamDB / Dfam 3.9 — embedded in RepeatMasker
        print("Dfam 3.9 library setup required")
        print("=" * 60)
        script_path = os.path.join(os.getcwd(), "configure_dfam39.sh")
        print(
            "Only the minimal Dfam partition 0 is present.\n"
            "Download the full set (0-16) and reconfigure RepeatMasker:\n\n"
            f"    cd {library_path}\n"
            "    curl -o 'dfam39_full.#1.h5.gz' \\\n"
            "      'https://dfam.org/releases/current/families/FamDB/"
            "dfam39_full.[0-16].h5.gz'\n"
            "    gunzip -f *.gz\n\n"
            f"A script template has been written to: {script_path}"
        )
        _generate_dfam39_config_script_legacy(library_path)
    print("=" * 60 + "\n")


def _generate_dfam39_config_script_legacy(library_path):
    """Write a legacy Dfam 3.9 configuration script for reference."""
    script_path = os.path.join(os.getcwd(), "configure_dfam39.sh")
    rm_share = library_path.replace("/Libraries/famdb", "")
    bin_dir = rm_share.replace("/share/RepeatMasker", "/bin")
    script_content = (
        "#!/bin/bash\n"
        "# Legacy Dfam 3.9 configuration script\n"
        f"cd {library_path}/\n"
        "curl -o 'dfam39_full.#1.h5.gz' "
        "'https://dfam.org/releases/current/families/FamDB/dfam39_full.[0-16].h5.gz'\n"
        "gunzip *.gz\n"
        f"mv {library_path}/min_init.0.h5 {library_path}/min_init.0.h5.bak\n"
        f"cd {rm_share}\n"
        "perl ./configure \\\n"
        f"    -libdir {library_path.replace('/famdb', '')} \\\n"
        f"    -trf_prgm {bin_dir}/trf \\\n"
        f"    -rmblast_dir {bin_dir} \\\n"
        f"    -hmmer_dir {bin_dir} \\\n"
        "    -default_search_engine rmblast\n"
        f"touch {library_path}/.earlgrey.config.complete\n"
    )
    with open(script_path, "w") as f:
        f.write(script_content)

