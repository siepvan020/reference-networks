#!/bin/bash

# Directory to store the RDS files
TARGET_DIR="data/rds"
URL_LIST="data/urls.txt"

# Download each file
while read -r url filename; do
  echo "Downloading $filename..."
  wget -q --show-progress "$url" -O "$TARGET_DIR/$filename"
done < "$URL_LIST"

echo "All RDS files downloaded to $TARGET_DIR"