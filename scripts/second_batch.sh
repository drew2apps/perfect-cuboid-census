#!/bin/bash
# second_step.py on the open bases of the census, in two halves; output second_step_open_100.txt (in the order of the list)
P=${PYTHON:-python3}
L=($(cat open_bases_100.txt)); n=${#L[@]}; h=$((n/2))
$P second_step.py "${L[@]:0:$h}" 2>/dev/null | grep --line-buffered -v "PARI stack" > /tmp/second_a.txt &
$P second_step.py "${L[@]:$h}" 2>/dev/null | grep --line-buffered -v "PARI stack" > /tmp/second_b.txt &
wait; cat /tmp/second_a.txt /tmp/second_b.txt > second_step_open_100.txt
