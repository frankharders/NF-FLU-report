process PARSE_FLUMUT {

    publishDir params.outdir, mode: 'copy'

    input:
    path flumut

    output:
    path "flumut_parsed.tsv", emit: flumut_parsed

    script:
    """
    python3 << 'EOF'
import re
import csv
from collections import defaultdict

grouped = defaultdict(set)

with open("${flumut}", "r") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        sample     = row['Sample'].strip()
        marker     = row['Marker'].strip()
        effect     = row['Effect'].strip()
        subtype    = row['Subtype'].strip()
        literature = row['Literature'].strip()

        for single_marker in marker.split(','):
            single_marker = single_marker.strip()
            if ':' not in single_marker:
                continue
            segment  = single_marker.split(':')[0]
            mutation = single_marker.split(':')[1]

            m_wt  = re.match(r'^([A-Za-z]+)', mutation)
            m_pos = re.search(r'([0-9]+)', mutation)
            m_mut = re.search(r'[0-9]+([A-Za-z]+)', mutation)

            if not (m_wt and m_pos and m_mut):
                continue

            wt_aa    = m_wt.group(1)
            position = m_pos.group(1)
            mut_aa   = m_mut.group(1)

            key = (sample, segment, wt_aa, position, mut_aa, subtype, literature)
            grouped[key].add(effect)

with open("flumut_parsed.tsv", "w", newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(["Sample", "Segment", "WT_aa", "Position", "Mut_aa", "Effects", "Subtype", "Literature"])
    for key in sorted(grouped.keys()):
        sample, segment, wt_aa, position, mut_aa, subtype, literature = key
        effects_joined = "; ".join(sorted(grouped[key]))
        writer.writerow([sample, segment, wt_aa, position, mut_aa, effects_joined, subtype, literature])

EOF
    """
}
