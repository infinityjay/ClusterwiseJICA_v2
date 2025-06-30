#!/bin/bash

cd /home/s4162315/cjica/ClusterwiseJICA_v2/
log_dir="./logs"
pattern='halted'

echo "Files containing 'Error' or 'halted':"

find "$log_dir" -type f -name "*.err" | while read -r file; do
    if grep -qE "$pattern" "$file"; then
        echo "$file"
    fi
done
