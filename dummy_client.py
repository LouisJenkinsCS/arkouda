#!/usr/bin/env python3

import arkouda as ak
import math
import time

if __name__ == "__main__":
    from scipy.spatial import distance
    import numpy as np
    import argparse, sys, gc, time

    parser = argparse.ArgumentParser(description="Example of cosine distance/similarity in arkouda")
    parser.add_argument('--server', default="localhost", help='server/Hostname of arkouda server')
    parser.add_argument('--port', type=int, default=5555, help='Port of arkouda server')
    args = parser.parse_args()
   
    print("Attempting to connect to server!")
    ak.enableVerbose()
    ak.connect(server=args.server, port=args.port, timeout=1)
    x = ak.zeros(1024 * 1024 * 1024, dtype=ak.int64)
    
    while True:
        try:
            x += 1
            s = x.sum()
            assert s % (1024 * 1024 * 1024) == 0, f"Bad sum: {s}"
            print("On ", s)
            time.sleep(1)
        except KeyboardInterrupt:
            pass         
    ak.disconnect()
    
    
