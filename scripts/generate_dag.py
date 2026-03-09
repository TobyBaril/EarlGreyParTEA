#!/usr/bin/env python3
"""
Generate DAG visualization for EarlGrey Pangenome Pipeline
This script creates a visual representation of the workflow directed acyclic graph
"""

import subprocess
import sys
import os
import shutil
from pathlib import Path


def generate_dag(
    snakefile,
    configfile,
    output_dir,
    cores=1,
    format="svg",
    filename_prefix="pipeline_dag"
):
    """
    Generate DAG visualization using Snakemake's built-in --dag option
    
    Args:
        snakefile: Path to Snakefile
        configfile: Path to config file
        output_dir: Directory to save DAG files
        cores: Number of cores for DAG generation (default 1)
        format: Output format - 'svg', 'png', or 'pdf' (default 'svg')
        filename_prefix: Prefix for output files (default 'pipeline_dag')
    """
    
    # Check if graphviz/dot is available
    if not shutil.which('dot'):
        print("[WARNING] Graphviz 'dot' command not found. Cannot generate DAG visualization.")
        print("          Install with: conda install graphviz  or  sudo apt install graphviz")
        return False
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Generate DOT format DAG
    dot_file = Path(output_dir) / f"{filename_prefix}.dot"
    image_file = Path(output_dir) / f"{filename_prefix}.{format}"
    
    try:
        # Get DAG in DOT format from Snakemake
        print(f"[INFO] Generating workflow DAG visualization...")
        cmd = [
            "snakemake",
            "--snakefile", str(snakefile),
            "--configfile", str(configfile),
            "--dag",
            "--cores", str(cores),
            "--quiet", "all"  # Suppress all output except the DAG
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False  # Don't fail if there are warnings
        )
        
        # Extract only the DOT graph from stdout (starts with "digraph")
        stdout_lines = result.stdout.split('\n')
        dot_content = []
        in_graph = False
        
        for line in stdout_lines:
            if line.strip().startswith('digraph'):
                in_graph = True
            if in_graph:
                dot_content.append(line)
        
        if not dot_content:
            print("[WARNING] No DAG graph generated. Check if workflow is valid.")
            return False
        
        # Save DOT file
        with open(dot_file, 'w') as f:
            f.write('\n'.join(dot_content))
        print(f"[INFO] Saved DAG in DOT format: {dot_file}")
        
        # Convert to image using Graphviz
        dot_cmd = ["dot", f"-T{format}", str(dot_file), "-o", str(image_file)]
        subprocess.run(dot_cmd, check=True, capture_output=True)
        print(f"[INFO] Saved DAG visualization: {image_file}")
        
        # Also generate a simplified rule graph (higher level view)
        rulegraph_dot = Path(output_dir) / f"{filename_prefix}_rulegraph.dot"
        rulegraph_image = Path(output_dir) / f"{filename_prefix}_rulegraph.{format}"
        
        cmd_rulegraph = [
            "snakemake",
            "--snakefile", str(snakefile),
            "--configfile", str(configfile),
            "--rulegraph",
            "--cores", str(cores),
            "--quiet", "all"
        ]
        
        result_rulegraph = subprocess.run(
            cmd_rulegraph,
            capture_output=True,
            text=True,
            check=False
        )
        
        # Extract only the DOT graph from stdout
        rulegraph_lines = result_rulegraph.stdout.split('\n')
        rulegraph_content = []
        in_graph = False
        
        for line in rulegraph_lines:
            if line.strip().startswith('digraph'):
                in_graph = True
            if in_graph:
                rulegraph_content.append(line)
        
        if rulegraph_content:
            with open(rulegraph_dot, 'w') as f:
                f.write('\n'.join(rulegraph_content))
            
            dot_cmd_rulegraph = ["dot", f"-T{format}", str(rulegraph_dot), "-o", str(rulegraph_image)]
            subprocess.run(dot_cmd_rulegraph, check=True, capture_output=True)
            print(f"[INFO] Saved rule graph visualization: {rulegraph_image}")
        
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"[WARNING] Failed to generate DAG: {e}")
        if e.stderr:
            print(f"          Error: {e.stderr}")
        return False
    except Exception as e:
        print(f"[WARNING] Unexpected error generating DAG: {e}")
        return False


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: generate_dag.py <snakefile> <configfile> <output_dir> [format] [cores]")
        print("  format: svg (default), png, or pdf")
        print("  cores: number of cores (default 1)")
        sys.exit(1)
    
    snakefile = sys.argv[1]
    configfile = sys.argv[2]
    output_dir = sys.argv[3]
    format = sys.argv[4] if len(sys.argv) > 4 else "svg"
    cores = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    
    success = generate_dag(snakefile, configfile, output_dir, cores, format)
    sys.exit(0 if success else 1)
