#!/bin/bash

####################################################################

# compute the Linkage Disequilibrium between  SNPs
#
#####################################################################
##<inputs>
population="Line_2"

### LD between adjucent SNPs
plink1 --file ${population} --ld-window 2 --ld-window-kb 10000 -ld-window-r2 0  --chr-set 40 --make-founders --r2 --out ${population}_adj

### LD between SNPs within 0-10kb distance
plink1 --file ${population} --ld-window 999999 --ld-window-kb 10000 -ld-window-r2 0  --chr-set 30 --make-founders --r2 --out ${population}_10kb

### compute the distance in kb

awk -F" " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS=" "}{print $0,abs($5-$2)}' ${population}_10kb.ld > dist_${population}.ld

