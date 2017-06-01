#!/bin/sh

if [[ "$1" == "--help" ]]; then
	echo -e "USAGE: ./bsub_gen.sh --outdir [OUTDIR] --fixfile [FIXFILE]"

elif [[ "$1" == "--outdir" && "$3" == "--fixfile" ]]; then
	export OUTDIR="$2"
	export FIXFILE="$4"
	envsubst < "$GERDACPTDIR/misc/bsub/template.bsub" > "$GERDACPTDIR/misc/bsub/$OUTDIR.bsub"

else
	echo -e "Bad input!"
fi

unset OUTDIR
unset FIXFILE
