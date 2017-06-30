#! /bin/bash

total=$(ls ../data/calculations/*.tar.gz | wc -l)
error=$(ls ../data/calculations/*_error.tar.gz | wc -l)
fin=$(echo "${total} - ${error}" | bc)
unfin=$(ls ../data/calculations/*/ -d | wc -l)
perc=$(echo "scale=0; ${fin} * 100 / 10157." | bc)

date
echo "${error} calculations terminated incorrectly"
echo "${unfin} have not finished running"
echo "${fin} (${perc} %) have completed successfully"
