#!/usr/bin/env python3

import sys
import sfacts as sf

if __name__ == "__main__":
    sf.World.load(sys.argv[1]).reassign_plurality_strain().dump(sys.argv[2])
