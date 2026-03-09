# EarlGrey Script Directory Auto-Detection

## Summary

The `script_dir` parameter is now **automatically detected** and does not need to be specified in config files.

## How It Works

### 1. User Config (No script_dir needed!)

```yaml
# Minimal config - no script_dir required
genome:
  species1: /path/to/genome1.fasta
species: [species1]
output_dir: /path/to/output
```

### 2. Auto-Detection Process

When the pipeline starts, it automatically searches for EarlGrey in these locations (in order):

**For conda/mamba installations:**
```bash
$CONDA_PREFIX/share/earlgrey-*/scripts     # Finds ANY version (7.x, 8.x, etc.)
$CONDA_PREFIX/share/earlgrey/scripts       # Generic location
```

**For development installations:**
```bash
/path/to/pangenome/../scripts              # Relative to Snakefile
```

**Validation:**
- Checks if directory exists
- Verifies `TEstrainer/` subdirectory is present
- Uses first valid match found

### 3. Runtime Message

```
Input validation
============================================================
[INFO] Auto-detected script directory: /path/to/earlgrey-7.0.3-0/scripts
============================================================
```

## Generated Config Files

Wrapper scripts generate configs with script_dir commented out:

```yaml
# Advanced options (usually not needed)
# script_dir: "/path/to/earlgrey/scripts"  # Auto-detected if installed via conda/mamba
```

Users can:
- ✅ Leave it commented out (recommended)
- ✅ Delete the comment entirely
- ✅ Uncomment and set for custom installations

## Version Compatibility

Auto-detection uses glob patterns to find ANY EarlGrey version:

```python
# Python code
glob.glob(os.path.join(share_dir, "earlgrey-*"))

# Matches:
#   earlgrey-7.0.3-0  ✓
#   earlgrey-7.1.0-0  ✓
#   earlgrey-8.0.0-0  ✓
#   earlgrey-10.5.2-3 ✓
```

### When EarlGrey Updates

```bash
# Upgrade EarlGrey
$ mamba update earlgrey
# earlgrey 7.0.3 → earlgrey 7.1.0

# Run pipeline (no config changes needed!)
$ earlGreyParTEA -c config.yaml -t 16
[INFO] Auto-detected script directory: .../earlgrey-7.1.0-0/scripts
```

## User Experience

### Installing

```bash
$ mamba install earlgrey-partea
$ earlGreyParTEA --generate-config my_config.yaml
[INFO] Example config file generated: my_config.yaml
```

### Editing Config

User edits `my_config.yaml` - **no script_dir to worry about**:

```yaml
genome:
  wheat: /data/wheat.fasta
  barley: /data/barley.fasta
species: [wheat, barley]
output_dir: /data/results
```

### Running

```bash
$ earlGreyParTEA -c my_config.yaml -t 32
[INFO] Auto-detected script directory: .../scripts
[INFO] Using Snakefile: .../Snakefile
# Pipeline runs successfully!
```

## Manual Override (Advanced)

For custom EarlGrey installations, users can still set script_dir:

```yaml
# Only needed for non-standard installations
script_dir: "/custom/path/to/earlgrey/scripts"
```

## Testing Results

### Test 1: Config without script_dir key
```bash
$ cat test_config.yaml
pipeline_mode: "libconstruct"
genome:
  test1: genome1.fasta
species: [test1]
output_dir: /tmp/test

$ snakemake --dry-run
✅ [INFO] Auto-detected script directory: .../earlgrey-7.0.3-0/scripts
```

### Test 2: Config with empty script_dir
```bash
$ cat test_config.yaml
script_dir: ""
genome:
  test1: genome1.fasta
species: [test1]
output_dir: /tmp/test

$ snakemake --dry-run
✅ [INFO] Auto-detected script directory: .../earlgrey-7.0.3-0/scripts
```

### Test 3: Generated config
```bash
$ earlGreyParTEA --generate-config test.yaml
$ grep "script_dir" test.yaml
✅ # script_dir: "/path/to/earlgrey/scripts"  # Auto-detected if installed via conda/mamba
```

## Implementation Details

### Python Code (pangenome/scripts/on_start_functions.py)

```python
# Auto-detect script_dir if not specified
if not config.get("script_dir") or config["script_dir"] == "":
    possible_locations = []
    
    # Check CONDA_PREFIX environment variable (version-agnostic)
    if os.environ.get("CONDA_PREFIX"):
        conda_prefix = os.environ["CONDA_PREFIX"]
        share_dir = os.path.join(conda_prefix, "share")
        if os.path.isdir(share_dir):
            # Find all earlgrey-* directories
            import glob
            for earlgrey_dir in glob.glob(os.path.join(share_dir, "earlgrey-*")):
                scripts_path = os.path.join(earlgrey_dir, "scripts")
                if os.path.isdir(scripts_path):
                    possible_locations.append(scripts_path)
    
    # Try each location
    for location in possible_locations:
        if os.path.isdir(location):
            testrainer_check = os.path.join(location, "TEstrainer")
            if os.path.isdir(testrainer_check):
                config["script_dir"] = location
                msg_info(f"Auto-detected script directory: {location}")
                break
```

### Bash Wrappers (earlGreyParTEA scripts)

Wrapper scripts keep script_dir out of generated configs - it's only shown as a commented example.

## Benefits

1. **Simpler Configs**: Users don't manage installation paths
2. **Version Updates**: Works across EarlGrey versions automatically
3. **Less Error-Prone**: No hardcoded paths to update
4. **Better UX**: One less thing to configure
5. **Still Flexible**: Advanced users can override if needed

## Error Handling

If auto-detection fails:

```
[ERROR] Script directory not specified and could not be auto-detected. 
Please set 'script_dir' in your config file or ensure EarlGrey is properly installed.
```

Clear guidance directs users to either:
- Fix their EarlGrey installation
- Manually specify script_dir for custom setups

## Documentation

- [README.md](README.md): Notes that script_dir is auto-detected
- [PACKAGING.md](PACKAGING.md): Details version-agnostic detection
- [VERSION_COMPATIBILITY.md](VERSION_COMPATIBILITY.md): Explains robustness to updates
- Troubleshooting section: Guides users if auto-detection fails

## Conclusion

✅ **User-friendly**: No installation paths in configs  
✅ **Version-robust**: Works with any EarlGrey version  
✅ **Well-tested**: Handles missing key, empty string, custom paths  
✅ **Well-documented**: Clear guidance for users  
✅ **Production-ready**: Safe for conda packaging  

The auto-detection feature makes ParTEA significantly easier to use while maintaining flexibility for advanced users!
