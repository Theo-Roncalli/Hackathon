#!/usr/bin/env bash

while read line; do
  fasterq-dump "$line"
done <ids.txt