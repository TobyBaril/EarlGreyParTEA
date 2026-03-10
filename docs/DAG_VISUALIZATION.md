# Workflow Visualization

This directory contains automatically generated DAG (Directed Acyclic Graph) visualizations of your EarlGrey pangenome pipeline execution, plus an example rulegraph for documentation.

## Example Rulegraph

**[example_rulegraph.svg](example_rulegraph.svg)** - A pre-generated rulegraph showing the library construction workflow. This serves as a visual reference of how ParTEA orchestrates parallel TE library construction across multiple genomes.

**Reading the rulegraph:**
- **Boxes** represent pipeline rules (steps)
- **Arrows** show dependencies (what needs to run before what)
- **Parallel branches** indicate rules that run simultaneously (per-genome operations)
- **Convergence points** show where results are combined (e.g., clustering)

## Files Generated During Pipeline Runs

### Full DAG (`dag_{mode}_mode.svg`)
Shows the complete directed acyclic graph of all jobs that will be executed, with wildcards resolved to actual values (species names, file paths, etc.). This is useful for:
- Understanding the exact sequence of operations
- Debugging job dependencies
- Seeing parallelization opportunities

### Rule Graph (`dag_{mode}_mode_rulegraph.svg`)
Shows a simplified high-level view of rule dependencies without resolving wildcards. This is useful for:
- Understanding the overall workflow structure
- Documentation and presentations
- Quick overview of the pipeline logic

### DOT Files (`.dot`)
Graph descriptions in DOT format (Graphviz language). These can be:
- Edited manually to customize visualization
- Converted to other formats (PNG, PDF, etc.)
- Processed by other graph tools

## Viewing the Visualizations

### SVG Files (default)
- Open in any modern web browser
- View in image viewers that support SVG
- Scale to any size without quality loss

### Converting to Other Formats
If you need PNG or PDF versions:
```bash
# Convert to PNG
dot -Tpng dag_full_mode.dot -o dag_full_mode.png

# Convert to PDF
dot -Tpdf dag_full_mode.dot -o dag_full_mode.pdf
```

## Configuration

Control DAG generation in your config file:
```yaml
# Workflow visualization (requires graphviz/dot installed)
generate_dag: True  # Set to False to skip DAG generation
dag_format: "svg"  # Options: "svg", "png", "pdf"
```

## Interpreting the Graphs

### Nodes (boxes)
- Each node represents a job or rule
- Node labels show the rule name and parameters
- Different pipeline modes show different active rules

### Edges (arrows)
- Arrows show dependencies between jobs
- A job depends on all jobs that point to it
- Jobs without dependencies can run in parallel

### Colors (in rule graph)
- Different colors indicate different rule types
- Library construction rules vs annotation rules
- Can help identify workflow stages

## Troubleshooting

**No visualizations generated?**
- Check if graphviz is installed: `dot -V`
- Install: `conda install graphviz` or `sudo apt install graphviz`
- Check console output for warning messages

**Graph too large or complex?**
- Use the rule graph (simplified view) instead of full DAG
- Consider breaking your analysis into smaller batches
- Convert to PNG/PDF at lower resolution if SVG rendering is slow

## Learn More

- Graphviz documentation: https://graphviz.org/
- Snakemake DAG visualization: https://snakemake.readthedocs.io/en/stable/executing/cli.html#visualization
