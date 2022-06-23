#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Check JSON file
"""

import os
import sys
import json
import numpy as np

def check_attribute(fn, **opts):
    verbose = opts.get('verbose', 0)

    # read json and get attribute
    with open(fn, 'r') as fp:
        json_obj = json.load(fp)
        attribute = json_obj.get('attribute')

    # check consistency with raw data
    status = True
    raw_fn = os.path.splitext(fn)[0] + '.raw'
    with open(raw_fn, 'r') as fp:
        for key in attribute.keys():
            item   = attribute[key]
            dtype  = np.dtype(item['datatype'])
            offset = item['offset']
            size   = item['size']
            ndim   = item['ndim']
            shape  = np.array(item['shape'])
            data   = np.atleast_1d(item['data'])
            # read
            fp.seek(offset)
            raw = np.fromfile(fp, dtype, size//dtype.itemsize)
            # check
            strrep1 = '{}'.format(data)
            strrep2 = '{}'.format(raw)
            msg = '    {:20s} : {:5s} --- {} <=> {}'
            status = status and (strrep1 == strrep2)
            if verbose > 0:
                if strrep1 == strrep2:
                    print(msg.format(key, 'Ok', strrep1, strrep2))
                else:
                    print(msg.format(key, 'Error', strrep1, strrep2))

    return status


if __name__ == '__main__':
    import argparse
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='JSON file to check',
                        nargs='*')
    parser.add_argument('-v', '--verbose', help='Verbosity level',
                        type=int, default=1)
    args = parser.parse_args()

    # numpy pretty print options
    np_print_formatter = {
        'int'   : '{:12d}'.format,
        'float' : '{:12.5e}'.format,
    }
    np.set_printoptions(formatter=np_print_formatter)

    # check attribute
    for fn in args.filename:
        if os.path.exists(fn):
            print('*** Checking {}'.format(fn))
            status = check_attribute(fn, verbose=args.verbose)
            if status:
                print('*** Success !')
            else:
                print('*** Error !')
