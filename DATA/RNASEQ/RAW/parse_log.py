import sys
import re

def parse_log_file(log_path):
    # regex patterns to catch the specific lines
    # matches "📄 File: SomeFile.txt"
    file_pattern = re.compile(r"📄 File:\s*(.+)")
    # matches ending percentage like "100.0%"
    download_pattern = re.compile(r"Downloading:.*\|\s*([0-9.]+%)\s*$")
    # matches "✅ Passed" or "❌ Failed"
    integrity_pattern = re.compile(r"Integrity Check:\s*(.+)")

    current_file = "Unknown_File"
    current_pct = "Unknown"

    print(f"{'FILENAME':<50} | {'STATUS':<8} | {'INTEGRITY Result'}")
    print("-" * 85)

    try:
        with open(log_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                
                # 1. Check for File Name
                file_match = file_pattern.search(line)
                if file_match:
                    current_file = file_match.group(1).strip()
                    current_pct = "0.0%" # Reset percentage for new file
                    continue

                # 2. Check for Download Status
                dl_match = download_pattern.search(line)
                if dl_match:
                    current_pct = dl_match.group(1).strip()
                    continue

                # 3. Check for Integrity (Trigger Output)
                integrity_match = integrity_pattern.search(line)
                if integrity_match:
                    result = integrity_match.group(1).strip()
                    print(f"{current_file:<50} | {current_pct:<8} | {result}")

    except FileNotFoundError:
        print(f"Error: Could not find file '{log_path}'")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 parse_log.py <logfile.txt>")
    else:
        parse_log_file(sys.argv[1])
