import os
import csv
import hashlib
import requests
import tarfile
import zipfile
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlparse

# Number of parallel downloads
MAX_WORKERS = 4

# Path to CSV file
CSV_FILE = "urls.csv"

# Path where files will be stored
DOWNLOAD_DIR = "../../data"

# Ensure download directory exists
os.makedirs(DOWNLOAD_DIR, exist_ok=True)

def calculate_checksum(file_path : str, hash_type : str):
    """Calculates the checksum of a file."""
    hash_func = {
        "SHA256": hashlib.sha256,
        "SHA1": hashlib.sha1,
        "MD5": hashlib.md5
    }.get(hash_type.upper())

    if not hash_func:
        print(f"Unsupported hash type: {hash_type}")
        return None

    hasher = hash_func()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hasher.update(chunk)
    return hasher.hexdigest()

def verify_checksum(file_path : str, expected_checksum : str, hash_type : str):
    """Compares the actual checksum with the expected checksum."""
    actual_checksum = calculate_checksum(file_path, hash_type)
    if actual_checksum == expected_checksum:
        print(f"{hash_type} checksum verified for {os.path.basename(file_path)}")
        return True
    else:
        print(f"{hash_type} checksum failed for {os.path.basename(file_path)}")
        print(f"Expected: {expected_checksum}")
        print(f"Got:      {actual_checksum}")
        os.remove(file_path) 
        print(f"Deleted corrupt file: {os.path.basename(file_path)}")
        return False

def extract_file(file_path : str, dir_path : str):
    """Extracts compressed files and deletes them after extraction."""
    try:
        if file_path.endswith((".tar.gz", ".tgz")):
            with tarfile.open(file_path, "r:gz") as tar:
                tar.extractall(dir_path)
        elif file_path.endswith(".tar.bz2"):
            with tarfile.open(file_path, "r:bz2") as tar:
                tar.extractall(dir_path)
        elif file_path.endswith(".tar.xz"):
            with tarfile.open(file_path, "r:xz") as tar:
                tar.extractall(dir_path)
        elif file_path.endswith(".tar"):
            with tarfile.open(file_path, "r") as tar:
                tar.extractall(dir_path)
        elif file_path.endswith(".zip"):
            with zipfile.ZipFile(file_path, "r") as zip_ref:
                zip_ref.extractall(dir_path)
        else:
            print(f"Unknown file type: {file_path}")
            return

        os.remove(file_path) 
        print(f"Extracted and deleted: {file_path}")

    except Exception as e:
        print(f"Extraction failed for {file_path}: {e}")

def download_and_verify(url : str, filename : str, expected_checksum : str, hash_type : str):
    """Downloads a file, verifies its checksum, and extracts it if valid."""
    
    dir_path = os.path.join(DOWNLOAD_DIR, os.path.splitext(filename)[0])
    file_path = os.path.join(dir_path, filename)
    os.makedirs(dir_path, exist_ok=True) 

    print(f"Downloading {filename} from {url}...")
    try:
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()

        with open(file_path, "wb") as f:
            for chunk in response.iter_content(8192):
                f.write(chunk)

        print(f"Downloaded {filename}")

    except requests.RequestException as e:
        print(f"Failed to download {filename}: {e}")
        return

    if expected_checksum and hash_type:
        if not verify_checksum(file_path, expected_checksum, hash_type):
            return  
    else:
        print(f"No checksum provided for {filename}, skipping verification.")

    # Extract file if valid
    extract_file(file_path, dir_path)

def main():
    """Reads CSV file and runs downloads in parallel."""
    with open(CSV_FILE, newline='', encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        tasks = []

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            for row in reader:
                url = row["url"].strip()
                checksum = row["checksum"].strip() if "checksum" in row else ""
                hash_type = row["hash_type"].strip().upper() if "hash_type" in row else ""
                filename = os.path.basename(urlparse(url).path)
                if not url:
                    print(f"Skipping invalid row: {row}")
                    continue

                # Submit tasks for parallel execution
                tasks.append(executor.submit(download_and_verify, url, filename, checksum, hash_type))

        # Wait for all downloads to complete
        for task in tasks:
            task.result()

    print("All downloads and extractions complete!")

if __name__ == "__main__":
    main()
