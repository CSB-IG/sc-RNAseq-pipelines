#!/bin/bash

# Function to verify checksum with different hash types (MD5, SHA1, SHA256)
verify_checksum() {
    local filename="$1"
    local expected_checksum="$2"
    local hash_type="$3"
    echo "filename: $filename"
    echo "checksum: $expected_checksum"
    echo "hash: $hash_type"

    if [[ ! -f "$filename" ]]; then
        echo "File $filename not found, skipping checksum verification."
        return 1
    fi

    echo "Verifying $hash_type checksum for $filename..."

    case "$hash_type" in
        SHA256) actual_checksum=$(sha256sum "$filename" | awk '{print $1}') ;;
        SHA1) actual_checksum=$(sha1sum "$filename" | awk '{print $1}') ;;
        MD5) actual_checksum=$(md5sum "$filename" | awk '{print $1}') ;;
        *)
            echo "Unsupported hash type: $hash_type for $filename, skipping verification."
            return 1
            ;;
    esac

    if [[ "$actual_checksum" == "$expected_checksum" ]]; then
        echo "$hash_type checksum verified successfully for $filename!"
        return 0
    else
        echo "$hash_type checksum verification failed for $filename!"
        echo "Expected: $expected_checksum"
        echo "Got:      $actual_checksum"
        echo "Deleting corrupt file: $filename"
        rm -f "$filename"
        return 1
    fi
}

# Function to download and extract
download_and_extract() {
    local url="$1"
    local hash_type="$2"
    local expected_checksum="$3"
    local filename="$(basename "$url")"
    local dirname="${filename%.*}"
    mkdir $dirname
    cd $dirname
    echo "Downloading $filename from $url..."
    wget -qc "$url"
    
    echo "hash"
    echo "$hash_type"
    echo "checksum"
    echo "$expected_checksum"
    if [[ $? -ne 0 ]]; then
        echo "Failed to download $filename" >> error.log
        return 1
    fi

    if [[ -n "$expected_checksum" && -n "$hash_type" ]]; then
        verify_checksum "$filename" "$expected_checksum" "$hash_type" || return 1
    else
        echo "No checksum provided for $filename, skipping verification." >> error.log
    fi

    echo "Extracting $filename..."
    
    case "$filename" in
        *.tar.gz|*.tgz) tar -xvzf "$filename" && rm -f "$filename" ;;
        *.tar.bz2) tar -xvjf "$filename" && rm -f "$filename" ;;
        *.tar.xz) tar -xvJf "$filename" && rm -f "$filename" ;;
        *.tar) tar -xvf "$filename" && rm -f "$filename" ;;
        *.zip) unzip -o "$filename" && rm -f "$filename" ;;
        *) echo "Unknown file type: $filename" >> error.log ;;
    esac

    echo "Completed: $filename"
}

export -f download_and_extract verify_checksum

cd ../../data || { echo "Failed to change directory"; exit 1; }

# Read CSV to obtain url, hash_type (MD5, sha256, sha1), and checksum
# download, verify checksum and extract in the background
echo "Processing downloads from CSV..."

while IFS=, read -r url hash_type checksum; do
    if [[ "$url" == "url" ]]; then
        continue
    fi

    download_and_extract "$url" "$hash_type" "$checksum" &

done < ../scripts/Hybrid/urls.csv

wait  

echo "All downloads and extractions complete!"

