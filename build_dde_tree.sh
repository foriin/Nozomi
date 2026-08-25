#!/usr/bin/env bash
set -euo pipefail

IN="${1:-sd01.txt}"
OUT="${2:-outputs/dde_tree}"
THREADS="${3:-16}"

SEED_DIR="${OUT}/seeds"
META="${OUT}/dde_metadata.tsv"
ALN="${OUT}/dde_all.mafft.fa"
TRIMMED="${OUT}/dde_all.clipkit.fa"
PREFIX="${OUT}/dde_all"

rm -rf "${SEED_DIR}"
mkdir -p "${SEED_DIR}"

printf "label\toriginal_label\tsuperfamily\n" > "${META}"

awk -v outdir="${SEED_DIR}" -v meta="${META}" '
function sanitize(x) {
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", x)
    gsub(/[^A-Za-z0-9_.-]+/, "_", x)
    return x
}

$0 ~ /^#+[^#].*#+$/ {
    family = $0
    gsub(/^#+/, "", family)
    gsub(/#+$/, "", family)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", family)

    if (family == "The superfamilies are arranged alphabetically") {
        family = ""
        next
    }

    family = sanitize(family)
    outfile = outdir "/" family ".fa"
    current_id = ""
    next
}

/^>/ {
    if (family == "") {
        print "ERROR: sequence found before a superfamily header: " $0 > "/dev/stderr"
        exit 1
    }

    original = substr($0, 2)
    id = original
    sub(/^\*/, "", id)
    id = sanitize(id)

    print ">" id >> outfile
    print id "\t" original "\t" family >> meta

    current_id = id
    next
}

{
    sequence = $0
    gsub(/[[:space:]]/, "", sequence)

    if (current_id != "" && sequence ~ /^[A-Za-z*-]+$/) {
        print sequence >> outfile
    }
}
' "${IN}"

DUPLICATES=$(
  grep -h '^>' "${SEED_DIR}"/*.fa |
    sed 's/^>//' |
    sort |
    uniq -d
)

if [[ -n "${DUPLICATES}" ]]; then
echo "ERROR: duplicated sequence identifiers:"
echo "${DUPLICATES}"
exit 1
fi

echo "Superfamily alignments:"
for f in "${SEED_DIR}"/*.fa; do
printf "%-25s %5d sequences\n" \
"$(basename "${f}" .fa)" \
"$(grep -c '^>' "${f}")"
done

seed_args=()

while IFS= read -r f; do
seed_args+=(--seed "${f}")
done < <(find "${SEED_DIR}" -name '*.fa' -type f | sort)

mafft \
--amino \
--maxiterate 1000 \
--thread "${THREADS}" \
"${seed_args[@]}" \
/dev/null \
> "${ALN}"

if command -v clipkit >/dev/null 2>&1; then
clipkit \
"${ALN}" \
-m smart-gap \
-o "${TRIMMED}"
else
  echo "ClipKIT not found; using the untrimmed MAFFT alignment."
cp "${ALN}" "${TRIMMED}"
fi

iqtree3 \
-s "${TRIMMED}" \
-st AA \
-m LG+G4 \
-B 1000 \
-alrt 1000 \
-bnni \
-nt AUTO \
-seed 1 \
-pre "${PREFIX}"

echo
echo "Alignment: ${TRIMMED}"
echo "Tree:      ${PREFIX}.contree"
echo "Metadata:  ${META}"