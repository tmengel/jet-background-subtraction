#!/usr/bin/env python3

import sys
import argparse

def MLsubtract(config):
    
    import jetflow
    jf = jetflow.JF(config)
    jf.Train()
    jf.Eval()
    print('Done')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--configpath", required=True, help="Path to config file for dataset", default=None)
    args = parser.parse_args()
    
    if args.configpath is None:
        print('Please provide a config file')
        sys.exit(1)
        
    print('Config path: {}'.format(args.configpath))
    MLsubtract(args.configpath)
    