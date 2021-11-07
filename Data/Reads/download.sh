#!/usr/bin/env bash

while read line; do
  echo "Downloading file ${line}..."
  fasterq-dump "${line}"
done <ids.txt
echo "Done"
