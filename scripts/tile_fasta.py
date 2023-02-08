#!/usr/bin/env python3

import sys
from math import ceil


if __name__ == "__main__":
    tile_len = int(sys.argv[1])
    overlap_len = int(sys.argv[2])
    inpath = sys.argv[3]
    with open(inpath) as f:
        raw_seqs = []
        acc = []
        for line in f:
            if line.startswith(">"):
                if acc:
                    raw_seqs.append("".join(acc))
                acc = []
            else:
                acc += line.strip()
        if acc:
            raw_seqs.append("".join(acc))

    tile_offset = tile_len - overlap_len
    for j, seq in enumerate(raw_seqs):
        seq_len = len(seq)
        n_tiles = ceil((seq_len - tile_len) / tile_offset)
        for i in range(n_tiles):
            start = i * tile_offset
            stop = i * tile_offset + tile_len
            print(f">{j}_{start}_{stop}")
            print(seq[start:stop])
        if stop != seq_len:
            stop = seq_len
            start = stop - tile_len
            print(f">{j}_{start}_{stop}")
            print(seq[start:stop])
        print(seq_len, n_tiles, file=sys.stderr)
