while IFS= read -r line; do
    sed -i 's/mem=10000/mem=6000/g;s/time=10:00:00/time=5:00:00/g' $line
done < ls.analyze.script.txt
