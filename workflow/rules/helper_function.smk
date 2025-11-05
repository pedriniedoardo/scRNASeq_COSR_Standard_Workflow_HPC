# --- Helper function ---
# This function takes wildcards, gets the fastq path string, splits it,
# strips whitespace, and marks each path as a directory dependency.
def get_input_fastq_dirs(wildcards):
    """
    Retrieves FASTQ directory paths for a sample and marks them
    as directory dependencies for Snakemake input.
    """
    # Get the comma-separated string of directory paths
    path_string = SAMPLES[wildcards.sample_name]["sample_path"]
    # Split the string into a list of paths
    dir_paths = path_string.split(',')
    # Strip leading/trailing whitespace and apply the directory() marker
    # This tells Snakemake to check if these directories exist
    return [directory(p.strip()) for p in dir_paths if p.strip()] # Added check for empty strings
