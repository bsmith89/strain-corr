#!/usr/bin/env python3

import sfacts as sf
import sys

if __name__ == "__main__":
    mgen_path = sys.argv[1]
    ref_path = sys.argv[2]
    multiplier = int(sys.argv[3])
    outpath = sys.argv[4]

    mgen = sf.Metagenotype.load(mgen_path)
    ref = sf.Metagenotype.load(ref_path).lift(lambda x: x * multiplier)
    out = sf.data.Metagenotype.concat(
        dict(mgen=mgen, ref=ref), dim="sample", rename=False
    ).mlift("fillna", 0)
    out.dump(outpath)
