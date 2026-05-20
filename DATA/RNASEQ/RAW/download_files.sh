#!/bin/bash

# ================= CONFIGURATION =================
# Source the token file if it exists (works in non-interactive shells like SLURM)
if [ -f "$HOME/.gdrive_token" ]; then
    source "$HOME/.gdrive_token"
fi

# Verify token is set
if [ -z "$GDRIVE_TOKEN" ]; then
    echo "❌ ERROR: GDRIVE_TOKEN is not set."
    echo "   Create ~/.gdrive_token containing:"
    echo "     export GDRIVE_TOKEN=\"ya29.a0Af...\""
    echo "   Then chmod 600 ~/.gdrive_token"
    exit 1
fi

# The Folder ID
FOLDER_ID="1bOFogbLTm_i-JQRfNKiLcOCNHOyBLUT7"
# =================================================

mkdir -p downloaded_files
cd downloaded_files

echo "Authenticating and fetching file list..."
echo "----------------------------------------------------"

# Export so python3 -c can pick them up from os.environ (safer than shell substitution into the heredoc)
export GDRIVE_TOKEN
export FOLDER_ID

python3 -c "
import urllib.request, json, sys, os, hashlib

token = os.environ['GDRIVE_TOKEN']
folder_id = os.environ['FOLDER_ID']
api_url = 'https://www.googleapis.com/drive/v3/files'

def get_all_files():
    files_list = []
    page_token = None

    while True:
        # Request ID, Name, Size, MimeType, AND md5Checksum
        params = f\"?q='{folder_id}'+in+parents+and+trashed=false&fields=nextPageToken,files(id,name,size,mimeType,md5Checksum)\"
        if page_token:
            params += f\"&pageToken={page_token}\"

        req = urllib.request.Request(api_url + params)
        req.add_header('Authorization', f'Bearer {token}')

        try:
            with urllib.request.urlopen(req) as response:
                data = json.load(response)
                files_list.extend(data.get('files', []))
                page_token = data.get('nextPageToken')
                if not page_token:
                    break
        except urllib.error.HTTPError as e:
            print(f'\n❌ Error listing files: {e}')
            if e.code == 401:
                print('   (Your Token may have expired. Please get a new one.)')
            sys.exit(1)
    return files_list

def calculate_local_md5(fname):
    # Calculates MD5 hash of local file to compare with Google's record
    hash_md5 = hashlib.md5()
    try:
        with open(fname, \"rb\") as f:
            for chunk in iter(lambda: f.read(4096), b\"\"):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    except FileNotFoundError:
        return None

def print_progress_bar(iteration, total, prefix='', suffix='', length=30, fill='█'):
    if total == 0: return
    percent = (\"{0:.1f}\".format(100 * (iteration / float(total))))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    sys.stdout.write(f'\r{prefix} |{bar}| {percent}% {suffix}')
    sys.stdout.flush()

def process_file(file_id, file_name, file_size, remote_md5):
    print(f'\n📄 File: {file_name}')

    # --- CHECK 1: Does file exist and match size? ---
    if os.path.exists(file_name):
        local_size = os.path.getsize(file_name)
        if file_size is not None and local_size == file_size:
            # --- CHECK 2: Verify Checksum (MD5) ---
            print('   • File exists. Verifying checksum...', end='')
            sys.stdout.flush()
            local_md5 = calculate_local_md5(file_name)

            if local_md5 == remote_md5:
                print(' ✅ Valid (Skipping download)')
                return
            else:
                print(' ❌ Checksum mismatch. Redownloading.')
        else:
            print(f'   • Size mismatch (Local: {local_size} vs Remote: {file_size}). Redownloading.')

    # --- DOWNLOAD START ---
    download_url = f'https://www.googleapis.com/drive/v3/files/{file_id}?alt=media'
    req = urllib.request.Request(download_url)
    req.add_header('Authorization', f'Bearer {token}')

    try:
        with urllib.request.urlopen(req) as source, open(file_name, 'wb') as out_file:
            if file_size is None:
                meta = source.info()
                file_size = int(meta.get('Content-Length', 0))

            downloaded = 0
            chunk_size = 1024 * 64 # 64KB chunks

            while True:
                chunk = source.read(chunk_size)
                if not chunk:
                    break
                out_file.write(chunk)
                downloaded += len(chunk)

                if file_size > 0:
                    print_progress_bar(downloaded, file_size, prefix='   • Downloading:', suffix='', length=30)

        print() # Newline

        # --- POST-DOWNLOAD VERIFICATION ---
        if remote_md5:
            final_md5 = calculate_local_md5(file_name)
            if final_md5 == remote_md5:
                print('   • Integrity Check: ✅ Passed')
            else:
                print('   • Integrity Check: ❌ FAILED (File corrupted)')

    except Exception as e:
        print(f'   ❌ Download Failed: {e}')

# --- Main Execution ---

all_files = get_all_files()

if not all_files:
    print('No files found.')
    sys.exit()

print(f'✅ Found {len(all_files)} files. Checking status...')
print('----------------------------------------------------')

for f in all_files:
    # Skip Google Docs
    if 'application/vnd.google-apps' in f.get('mimeType', ''):
        continue

    size = int(f.get('size', 0)) if 'size' in f else None
    md5 = f.get('md5Checksum')

    process_file(f['id'], f['name'], size, md5)

print('\n----------------------------------------------------')
print('🎉 All operations complete.')
"
