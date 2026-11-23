import gzip
import re
import sys

import numpy as np

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
CG = re.compile(b"CG")


def main():
    fasta, out = sys.argv[1], sys.argv[2]
    seqs, name, chunks = {}, None, []
    with gzip.open(fasta, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name in CHROMS:
                    seqs[name] = b"".join(chunks).upper()
                name, chunks = line[1:].split()[0].decode(), []
            elif name in CHROMS:
                chunks.append(line.rstrip())
    if name in CHROMS:
        seqs[name] = b"".join(chunks).upper()

    arrays, offset = {}, 0
    for c in CHROMS:
        pos = np.fromiter((m.start() for m in CG.finditer(seqs[c])), dtype=np.int32)
        arrays[f"pos_{c}"] = pos
        arrays[f"off_{c}"] = np.int64(offset)
        offset += len(pos)
        print(c, len(pos), flush=True)
    print("total CpGs", offset)
    np.savez(out, **arrays)


if __name__ == "__main__":
    main()
