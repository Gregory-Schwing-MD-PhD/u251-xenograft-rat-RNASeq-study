import os
import re
import csv

# CONFIGURATION
# ---------------------------------------------------------
# 1. Get the directory where this script resides (ANALYSIS)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# 2. Define RAW_DIR relative to the script
#    Going up one level (..) to base, then down into DATA/METHYLARRAY/RAW
RAW_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "../DATA/METHYLARRAY/RAW"))

# Output filename (Saved in the same dir as this script)
OUTPUT_FILE = "samplesheet_methylarray.csv"

# Sample Mapping
SAMPLE_MAP = {
    "205648300021_R01C01": "IL66B",
    "205648300021_R02C01": "IL67B",
    "205648300021_R03C01": "IL68B",
    "205648300021_R04C01": "IL69B",
    "205648300021_R05C01": "NL70B",
    "205648300021_R06C01": "NL71B",
    "205648300021_R07C01": "N269B",
    "205648300021_R08C01": "C2B"
}
# ---------------------------------------------------------

def generate_samplesheet():
    samples = []
    pattern = re.compile(r"(\d+)_(\w+)_Grn.idat")

    print(f"Script Directory: {SCRIPT_DIR}")
    print(f"Target RAW Directory: {RAW_DIR}")

    if not os.path.exists(RAW_DIR):
        print(f"Error: Directory {RAW_DIR} does not exist.")
        return

    # Scan files
    for filename in os.listdir(RAW_DIR):
        match = pattern.match(filename)
        if match:
            sentrix_id = match.group(1)
            sentrix_position = match.group(2)
            map_key = f"{sentrix_id}_{sentrix_position}"
            
            # Construct absolute paths for the CSV 
            # (Best for Nextflow, even if input logic was relative)
            red_path = os.path.join(RAW_DIR, f"{sentrix_id}_{sentrix_position}_Red.idat")
            grn_path = os.path.join(RAW_DIR, filename)
            
            if map_key in SAMPLE_MAP:
                sample_name = SAMPLE_MAP[map_key]
                
                if os.path.exists(red_path):
                    samples.append({
                        "sample": sample_name,
                        "array": "EPIC",
                        "red_channel": red_path,
                        "green_channel": grn_path
                    })
                else:
                    print(f"Warning: Red channel missing for {sample_name}")

    if not samples:
        print("No valid samples found.")
        return

    # Determine Output Path
    output_path = os.path.join(SCRIPT_DIR, OUTPUT_FILE)

    # Write to CSV
    header = ["sample", "array", "red_channel", "green_channel"]
    
    try:
        with open(output_path, mode='w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            writer.writerows(samples)
            
        print(f"Successfully created: {output_path}")
        print(f"Total samples: {len(samples)}")
        print("-" * 40)
            
    except IOError as e:
        print(f"Error writing file: {e}")

if __name__ == "__main__":
    generate_samplesheet()
