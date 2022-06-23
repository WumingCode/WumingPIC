#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import io
import subprocess
from subprocess import PIPE
import contextlib
import difflib
import numpy as np
import h5py


def unittest():
    # jsonio
    jsonio_args = dict(
        test_cmd = './jsonio_test',
        output   = 'jsonio_test.json',
        expected = 'jsonio_expected.json')
    jsonio_test(**jsonio_args)
    # mpiio
    mpiio_args = dict(
        test_cmd1 = 'mpiexec -n 4 ./mpiio_test1',
        test_cmd2 = 'mpiexec -n 4 ./mpiio_test2',
        hdffile   = 'mpiio_test.h5',
        rawfile   = 'mpiio_test.raw')
    mpiio_test(**mpiio_args)


def check_stderr(stderr, cmd):
    if isinstance(stderr, list):
        stderr = ''.join(stderr)

    if not stderr == '':
        print('*** Some error detected in executing: {} ***'.format(cmd))
        print(stderr)
        raise RuntimeError


def check_diff(text1, text2, operation):
    if isinstance(text1, str):
        text1 = text1.split('\n')
    if isinstance(text2, str):
        text2 = text2.split('\n')
    diff = [d for d in difflib.unified_diff(text1, text2)]
    if len(diff) != 0:
        print('*** Some error detected in {} ***'.format(operation))
        print('*** begin diff ***')
        sys.stdout.writelines(diff)
        print('*** end diff ***')
        raise RuntimeError


def exec_process(cmd):
    if isinstance(cmd, str):
        cmd = cmd.split()
    proc = subprocess.Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout = [l.decode('utf-8') for l in proc.stdout.readlines()]
    stderr = [l.decode('utf-8') for l in proc.stderr.readlines()]
    return stdout, stderr


def jsonio_test(test_cmd, output, expected):
    stdout, stderr = exec_process(test_cmd)

    # check stderr
    check_stderr(stderr, test_cmd)

    # compare output and expected result
    json1 = open(output, 'r').readlines()
    json2 = open(expected, 'r').readlines()
    check_diff(json1, json2, 'jsonio_test comparison')

    print('***')
    print('jsonio test seems to be succesful !')
    print('***')


def mpiio_test(test_cmd1, test_cmd2, rawfile, hdffile):
    stdout1, stderr1 = exec_process(test_cmd1)
    stdout2, stderr2 = exec_process(test_cmd2)

    # check stderr
    check_stderr(stderr1, test_cmd1)
    check_stderr(stderr2, test_cmd2)

    # compare outputs
    check_diff(stdout1, stdout2, 'mpiio_test1 and mpiio_test2 comparison')

    # read in python
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        create_hdf5(hdffile, rawfile)
    buf.seek(0)
    stdout0 = buf.readlines()

    # compare outputs
    check_diff(stdout0, stdout1, 'comparison between python and test1')
    check_diff(stdout0, stdout2, 'comparison between python and test2')

    print('***')
    print('mpiio test seems to be succesful !')
    print('HDF5 file {} has been generated'.format(hdffile))
    print('***')


def create_hdf5(hdffile, rawfile):
    # open and read small data
    with open(rawfile, 'rb') as fp:
        endian = np.fromfile(fp, 'i4', 1)[0]
        if endian == 1:
            bo = '<'
        elif endian == 16777216:
            bo = '>'
        else:
            errmsg = 'Error: unrecognized endian flag: {}'.format(endian)
            raise RuntimeError(errmsg)
        nproc  = np.fromfile(fp, bo+'i4', 1)[0]
        dims   = np.fromfile(fp, bo+'i4', 2)
        nx     = np.fromfile(fp, bo+'i4', 1)[0]
        ny     = np.fromfile(fp, bo+'i4', 1)[0]
        var_i4 = np.fromfile(fp, bo+'i4', 1)[0]
        var_i8 = np.fromfile(fp, bo+'i8', 1)[0]
        var_r4 = np.fromfile(fp, bo+'f4', 1)[0]
        var_r8 = np.fromfile(fp, bo+'f8', 1)[0]
        ldesc  = np.fromfile(fp, bo+'i4', 1)[0]
        desc   = fp.read(ldesc).decode('utf-8')
        offset = fp.tell()
        print('{:30s} : {:8d}'.format('# endian flag', endian))
        print('{:30s} : {:8d}'.format('# MPI process', nproc))
        print('{:30s} : {:8d}'.format('# MPI process in x', dims[0]))
        print('{:30s} : {:8d}'.format('# MPI process in y', dims[1]))
        print('{:30s} : {:8d}'.format('# local grid in x', nx))
        print('{:30s} : {:8d}'.format('# local grid in y', ny))
        print('{:30s} : {:8d}'.format('# global grid in x', nx*dims[0]))
        print('{:30s} : {:8d}'.format('# global grid in y', ny*dims[1]))
        print('{:30s} : {:8d}'.format('# integer(4)', var_i4))
        print('{:30s} : {:8d}'.format('# integer(8)', var_i8))
        print('{:30s} : {:10.4f}'.format('# real(4)', var_r4))
        print('{:30s} : {:10.4f}'.format('# real(8)', var_r8))
        print('{:30s} : {:8d}'.format('# character(*) len', ldesc))
        print('{:30s} : {:s}'.format('# character(*)', desc))

    # create hdf5 and write dataset for large data
    with h5py.File(hdffile, 'w') as fp:
        n = nx*ny * nproc
        s = (ny*dims[1], nx*dims[0],)

        # int32
        print('{:30s} : {:16d}'.format('# offset for integer(4)', offset))
        size = 4 * n
        ext  = ((rawfile, offset, size),)
        fp.create_dataset('dat_i4', shape=s, dtype=bo+'i4', external=ext)
        offset = offset + size

        # int64
        print('{:30s} : {:16d}'.format('# offset for integer(8)', offset))
        size = 8 * n
        ext  = ((rawfile, offset, size),)
        fp.create_dataset('dat_i8', shape=s, dtype=bo+'i8', external=ext)
        offset = offset + size

        # float32
        print('{:30s} : {:16d}'.format('# offset for real(4)', offset))
        size = 4 * n
        ext  = ((rawfile, offset, size),)
        fp.create_dataset('dat_r4', shape=s, dtype=bo+'f4', external=ext)
        offset = offset + size

        # int64
        print('{:30s} : {:16d}'.format('# offset for real(8)', offset))
        size = 8 * n
        ext  = ((rawfile, offset, size),)
        fp.create_dataset('dat_r8', shape=s, dtype=bo+'f8', external=ext)
        offset = offset + size


if __name__ == '__main__':
    unittest()
