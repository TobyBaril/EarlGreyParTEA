# EarlGrey ParTEA Packaging Guide

## Overview

**Repository:** https://github.com/TobyBaril/EarlGreyParTEA (separate from main EarlGrey repo)

EarlGrey ParTEA should be packaged as a **separate conda/mamba package** that depends on EarlGrey. This design:
- Keeps the pangenome pipeline separate from single-genome EarlGrey
- Allows independent versioning and updates
- Ensures compatibility through dependency management
- Works robustly across EarlGrey version updates
- Maintains separate GitHub repositories for clear organization

## Package Structure

### Package 1: earlgrey (existing package)
```
$CONDA_PREFIX/
├── bin/
│   └── earlGrey                    # Original EarlGrey command
└── share/
    └── earlgrey-X.Y.Z-N/           # Version-specific directory
        ├── scripts/
        │   ├── TEstrainer/
        │   ├── repeatCraft/
        │   └── ... (all EarlGrey scripts)
        └── ... (other EarlGrey files)
```

### Package 2: earlgrey-partea (new separate package)
```
$CONDA_PREFIX/
├── bin/
│   ├── earlGreyParTEA              # Main wrapper script
│   ├── earlGreyParTEA_LibConstruct # Library construction wrapper
│   └── earlGreyParTEA_AnnotationOnly # Annotation-only wrapper
└── share/
    └── earlgrey-partea-X.Y.Z-N/   # Version-specific directory
        ├── Snakefile
        ├── config/
        │   └── config.yaml
        ├── rules/
        │   ├── annotate.smk
        │   ├── annotate_simple.smk
        │   ├── clustering.smk
        │   └── lib_construct.smk
        └── scripts/
            ├── on_start_functions.py
            ├── generate_dag.py
            └── ... (other scripts)
```

## Conda Recipe Example

### meta.yaml for earlgrey-partea
```yaml
{% set name = "earlgrey-partea" %}
{% set version = "1.0.0" %}

package:
  name: {{ name|lower }}
  version: {{ version }}

source:
  url: https://github.com/TobyBaril/EarlGreyParTEA/archive/refs/tags/v{{ version }}.tar.gz
  sha256: YOUR_SHA256_HERE

build:
  number: 0
  noarch: python

requirements:
  host:
    - python >=3.9
  run:
    - python >=3.9
    - earlgrey >=7.0.3        # Dependency on base EarlGrey package
    - snakemake >=8.0
    - cd-hit
    - graphviz               # For DAG visualization

test:
  commands:
    - earlGreyParTEA --help
    - earlGreyParTEA_LibConstruct --help
    - earlGreyParTEA_AnnotationOnly --help

about:
  home: https://github.com/TobyBaril/EarlGreyParTEA
  license: OSI-approved BSD License
  summary: 'Pangenome transposable element annotation pipeline using EarlGrey'
  description: |
    EarlGrey ParTEA extends EarlGrey for pangenome-scale transposable element
    annotation. It processes multiple genomes in parallel, performs cross-species
    clustering, and generates comprehensive repeat annotations.

extra:
  recipe-maintainers:
    - TobyBaril
```

### build.sh for earlgrey-partea
```bash
#!/bin/bash

# Install wrapper scripts to bin/
mkdir -p $PREFIX/bin
install -m 755 earlGreyParTEA $PREFIX/bin/
install -m 755 earlGreyParTEA_LibConstruct $PREFIX/bin/
install -m 755 earlGreyParTEA_AnnotationOnly $PREFIX/bin/

# Install Snakemake workflow to share directory
mkdir -p $PREFIX/share/${PKG_NAME}-${PKG_VERSION}
cp Snakefile $PREFIX/share/${PKG_NAME}-${PKG_VERSION}/
cp -r rules $PREFIX/share/${PKG_NAME}-${PKG_VERSION}/
cp -r scripts $PREFIX/share/${PKG_NAME}-${PKG_VERSION}/
cp -r config $PREFIX/share/${PKG_NAME}-${PKG_VERSION}/
```

## Version-Agnostic Path Detection

### How It Works

The wrapper scripts and Python code use **glob patterns** to find installed files regardless of version:

#### In Wrapper Scripts (bash)
```bash
# Check for versioned earlgrey installations
for earlgrey_dir in "$CONDA_PREFIX/share/earlgrey"-*; do
    if [ -d "$earlgrey_dir" ]; then
        POSSIBLE_LOCATIONS+=("$earlgrey_dir/pangenome/Snakefile")
    fi
done
```

This finds:
- `share/earlgrey-7.0.3-0/` ✓
- `share/earlgrey-7.1.0-0/` ✓
- `share/earlgrey-8.0.0-1/` ✓
- Any future version automatically ✓

#### In Python Code (on_start_functions.py)
```python
import glob

share_dir = os.path.join(conda_prefix, "share")
for earlgrey_dir in glob.glob(os.path.join(share_dir, "earlgrey-*")):
    scripts_path = os.path.join(earlgrey_dir, "scripts")
    if os.path.isdir(scripts_path):
        possible_locations.append(scripts_path)
```

### Robustness Features

1. **Multiple Search Paths**: Checks both versioned and generic directories
2. **Validation**: Verifies presence of expected files (e.g., TEstrainer)
3. **Fallback**: Works in development mode (relative paths)
4. **Clear Errors**: Reports all searched locations if not found

## Installation Workflow

### For Users

```bash
# Install both packages (earlgrey-partea depends on earlgrey)
mamba install -c conda-forge -c bioconda earlgrey-partea

# Verify installation
earlGreyParTEA --help
```

### For Developers

```bash
# Clone the EarlGreyParTEA repository
git clone https://github.com/TobyBaril/EarlGreyParTEA.git
cd EarlGreyParTEA

# Make wrapper scripts executable
chmod +x earlGreyParTEA*

# Add to PATH (or use absolute path)
export PATH="$PWD:$PATH"

# Scripts auto-detect relative paths in development
earlGreyParTEA --generate-config example_config.yaml
```

## Version Compatibility Testing

### Test Across EarlGrey Versions

```bash
# Test with EarlGrey 7.0.3
mamba create -n test-7.0.3 earlgrey=7.0.3 earlgrey-partea
conda activate test-7.0.3
earlGreyParTEA --generate-config test.yaml

# Test with EarlGrey 7.1.0
mamba create -n test-7.1.0 earlgrey=7.1.0 earlgrey-partea
conda activate test-7.1.0
earlGreyParTEA --generate-config test.yaml
```

### Expected Behavior

- ✅ Auto-detects correct script directory for each version
- ✅ Finds Snakefile regardless of version suffix
- ✅ Validates presence of required EarlGrey components
- ✅ Reports clear errors if dependencies missing

## Upgrade Path

When EarlGrey updates from 7.0.3 to 7.1.0:

1. **No code changes needed** - glob patterns find new version
2. **No config changes needed** - paths auto-detected
3. **No user action needed** - just update packages

```bash
# Update both packages
mamba update earlgrey earlgrey-partea

# Continues working automatically
earlGreyParTEA -c config.yaml -g genomes.txt -s species.txt -t 48
```

## Troubleshooting

### If Scripts Can't Find EarlGrey

```bash
# Check installed packages
conda list | grep earlgrey

# Check directory structure
ls -la $CONDA_PREFIX/share/earlgrey*

# Manually specify script_dir in config
echo "script_dir: $CONDA_PREFIX/share/earlgrey-7.0.3-0/scripts" >> config.yaml
```

### If DAG Generation Fails

```bash
# Install graphviz
mamba install graphviz

# Or disable DAG generation
echo "generate_dag: False" >> config.yaml
```

## Benefits of Separate Packaging

1. **Independent Development**
   - Update parTEA without releasing new EarlGrey
   - Fix bugs in pangenome pipeline independently
   
2. **Version Flexibility**
   - Users can pin EarlGrey version
   - ParTEA compatible across EarlGrey versions
   
3. **Clear Dependencies**
   - Conda ensures EarlGrey installed first
   - Explicit version requirements in meta.yaml
   
4. **Smaller Footprint**
   - Users not needing pangenome don't install it
   - Separate package has fewer dependencies

5. **Clearer Maintenance**
   - Two GitHub repos or clear separation
   - Different release cycles
   - Separate issue tracking

## Recommended Repository Structure

### EarlGrey Repository (Dependency)
```
https://github.com/TobyBaril/EarlGrey

EarlGrey/
├── earlGrey                    # Main EarlGrey script
├── scripts/                    # All EarlGrey scripts
│   ├── TEstrainer/
│   ├── repeatCraft/
│   └── ...
└── conda/                      # Conda packaging for EarlGrey
    └── meta.yaml
```

### EarlGreyParTEA Repository (Separate Package)
```
https://github.com/TobyBaril/EarlGreyParTEA

EarlGreyParTEA/
├── earlGreyParTEA              # Main wrapper script
├── earlGreyParTEA_LibConstruct # Library construction wrapper
├── earlGreyParTEA_AnnotationOnly # Annotation-only wrapper
├── Snakefile                   # Workflow definition
├── rules/                      # Snakemake rule files
│   ├── annotate.smk
│   ├── annotate_simple.smk
│   ├── clustering.smk
│   └── lib_construct.smk
├── scripts/                    # Python helper scripts
│   ├── on_start_functions.py
│   ├── generate_dag.py
│   └── ...
├── config/                     # Config templates
│   └── config.yaml
├── conda/                      # Conda packaging
│   ├── meta.yaml
│   └── build.sh
└── README.md                   # Documentation
```

**Installation Paths:**
- When installed via conda, files are placed directly in `$PREFIX/share/earlgrey-partea-X.Y.Z-N/` (no subdirectories)

## Summary

The current implementation is **fully compatible** with separate conda packaging:

✅ Version-agnostic path detection with glob patterns  
✅ Works across EarlGrey 7.0.3, 7.1.0, 8.x.x automatically  
✅ Clear dependency structure (partea depends on earlgrey)  
✅ Auto-detection with multiple fallback paths  
✅ Development mode support for testing  
✅ Comprehensive error reporting  

No further code changes needed for packaging - just create the conda recipes!
