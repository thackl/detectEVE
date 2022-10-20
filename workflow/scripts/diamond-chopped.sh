#!/usr/bin/env bash
set -e

args=("$@")
window=60000
step=50000

for ((i = 1; i <= $#; i++ )); do
    j=$(( $i + 1 ))
    h=$(( $i - 1 ))
    case "${!i}" in
        -o | --out)
            out="${!j}"
            o6_chopped="$out.chopped"
            args[$i]="$o6_chopped"
            ;;
        -q | --query)
            in="${!j}"
            fa_chopped=$( basename $in ).chopped
            args[$i]="$fa_chopped"
            ;;
        -W | --window)
            window="${!j}"
            args[$h]=""
            args[$i]=""
            ;;
        -S | --step)
            step="${!j}"
            args[$h]=""
            args[$i]=""
            ;;
    esac
done;

#[ ! -n "$step" ] && { echo "-S/--step required" >&2; exit 1;}
#[ ! -n "$window" ] && { echo "-W/--window required" >&2; exit 1;}

echo "Chopping sequences"
(set -x;
 seqkit sliding --window $window --step $step -o $fa_chopped --greedy $in;
)


echo ""
echo "Running diamond"
(set -x;
 diamond ${args[@]}
)

echo ""
echo "Unchopping diamond results"
(set -x;
 perl -ane '($start) = $F[0] =~ /_sliding:(\d+)-(\d+)$/; $F[0] =~ s/_sliding:\d+-\d+$//; $F[6]+=$start-1, $F[7]+=$start-1; print join("\t", @F), "\n";' $o6_chopped > $out
)
echo ""
echo "Removing temporary files"
(set -x;
 rm $fa_chopped $o6_chopped
)
