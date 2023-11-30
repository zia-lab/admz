#!/usr/bin/env python3

import os
import gzip
import shutil
import subprocess
import numpy as np

def list_to_symmetric_matrix(lst, dtype = np.float_):
    # Calculate the size of the matrix
    N = int((np.sqrt(8*len(lst) + 1) - 1) / 2)

    # Initialize an NxN matrix with zeros
    matrix = np.zeros((N, N), dtype=dtype)

    # Fill in the upper triangular part of the matrix
    idx = 0
    for i in range(N):
        for j in range(0, i+1):
            matrix[i, j] = lst[idx]
            idx += 1

    # Since the matrix is symmetric, copy the upper part to the lower part
    matrix = matrix + matrix.T - np.diag(matrix.diagonal())

    return matrix

def load_gz_matrix(fname, upper_triangular=False, assume_real = True):
    '''
    Loads a gzip compressed matrix from a file in which the matrix elements are
    stored in a single column. With their real and imaginary parts given.
    '''
    with gzip.open(fname, 'rt') as f:
        # Assuming the file contains one number per line
        data = np.loadtxt(f)
    # assuming that the matrix is real
    if assume_real:
        the_loaded_matrix = data[:,0]
        dtype = np.float_
    else:
        the_loaded_matrix = data[:,0] + 1j*data[:,1]
        dtype = np.complex_
    if upper_triangular:
        the_full_matrix = list_to_symmetric_matrix(the_loaded_matrix, dtype=dtype)
    else:
        dim = int(np.sqrt(len(the_loaded_matrix)))
        the_full_matrix = np.reshape(the_loaded_matrix, 
                                     (dim, dim))
    return the_full_matrix

def extract_variable_names_from_format(template_string):
    # Create a regular expression pattern for identifying format variables
    # This pattern looks for {variable}
    pattern = r'\{([^{}]+)\}'

    # Find all matches in the template string
    matches = re.findall(pattern, template_string)

    # Filtering out any numerical indices or formatting specifications
    variable_names = [match for match in matches if not match.isdigit() and ':' not in match]

    return variable_names

def execute_and_print_realtime(command, workdir='./', verbose=True):
    # Start the subprocess and get its output in real-time
    command = f"module load lapack; {command}"
    print(command)
    process = subprocess.Popen(command, cwd= workdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # Read the output line by line as it is generated
    while True:
        output = process.stdout.readline()
        if output == b'' and process.poll() is not None:
            break
        if output:
            if verbose:
                print(output.strip().decode())

    rc = process.poll()
    return rc

param_template = '''1.0 1.0 1.0 0.0 0.0 0.0
{F2} {F4} {F6} {zeta}
{alpha} {beta} {gamma}
{T2} {T3} {T4} {T6} {T7} {T8}
{Bx} {By} {Bz}
'''

cfp_template = '''0 0 {B00}
1 -1 {B1m1}
1 0 {B10}
1 1 {B11}
2 -2 {B2m2}
2 -1 {B2m1}
2 0 {B20}
2 1 {B21}
2 2 {B22}
3 -3 {B3m3}
3 -2 {B3m2}
3 -1 {B3m1}
3 0 {B30}
3 1 {B31}
3 2 {B32}
3 3 {B33}
4 -4 {B4m4}
4 -3 {B4m3}
4 -2 {B4m2}
4 -1 {B4m1}
4 0 {B40}
4 1 {B41}
4 2 {B42}
4 3 {B43}
4 4 {B44}
5 -5 {B5m5}
5 -4 {B5m4}
5 -3 {B5m3}
5 -2 {B5m2}
5 -1 {B5m1}
5 0 {B50}
5 1 {B51}
5 2 {B52}
5 3 {B53}
5 4 {B54}
5 5 {B55}
6 -6 {B6m6}
6 -5 {B6m5}
6 -4 {B6m4}
6 -3 {B6m3}
6 -2 {B6m2}
6 -1 {B6m1}
6 0 {B60}
6 1 {B61}
6 2 {B62}
6 3 {B63}
6 4 {B64}
6 5 {B65}
6 6 {B66}
7 -7 {B7m7}
7 -6 {B7m6}
7 -5 {B7m5}
7 -4 {B7m4}
7 -3 {B7m3}
7 -2 {B7m2}
7 -1 {B7m1}
7 0 {B70}
7 1 {B71}
7 2 {B72}
7 3 {B73}
7 4 {B74}
7 5 {B75}
7 6 {B76}
7 7 {B77}
'''

def param_inp_file_maker(numE, F2=0., F4=0., F6=0., zeta=0.,
                         alpha=0., beta=0., gamma=0.,
                         T2=0., T3=0., T4=0., T6=0., T7=0., T8=0.,
                         Bx=0., By=0., Bz=0.):
    file_content = param_template.format(F2=F2, F4=F4, F6=F6, zeta=zeta, alpha=alpha, beta=beta, gamma=gamma, T2=T2, T3=T3, T4=T4, T6=T6, T7=T7, T8=T8, Bx=Bx, By=By, Bz=Bz)
    with open('./f%d/param.inp' % numE, 'w') as f:
        f.write(file_content)
    return file_content

def cfp_inp_file_maker(numE, Bkq):
    default_Bkq = {'B00': 0.0,'B1m1': 0.0,'B10': 0.0,'B11': 0.0,'B2m2': 0.0,'B2m1': 0.0,'B20': 0.0,'B21': 0.0,'B22': 0.0,'B3m3': 0.0,'B3m2': 0.0,'B3m1': 0.0,'B30': 0.0,'B31': 0.0,'B32': 0.0,'B33': 0.0,'B4m4': 0.0,'B4m3': 0.0,'B4m2': 0.0,'B4m1': 0.0,'B40': 0.0,'B41': 0.0,'B42': 0.0,'B43': 0.0,'B44': 0.0,'B5m5': 0.0,'B5m4': 0.0,'B5m3': 0.0,'B5m2': 0.0,'B5m1': 0.0,'B50': 0.0,'B51': 0.0,'B52': 0.0,'B53': 0.0,'B54': 0.0,'B55': 0.0,'B6m6': 0.0,'B6m5': 0.0,'B6m4': 0.0,'B6m3': 0.0,'B6m2': 0.0,'B6m1': 0.0,'B60': 0.0,'B61': 0.0,'B62': 0.0,'B63': 0.0,'B64': 0.0,'B65': 0.0,'B66': 0.0,'B7m7': 0.0,'B7m6': 0.0,'B7m5': 0.0,'B7m4': 0.0,'B7m3': 0.0,'B7m2': 0.0,'B7m1': 0.0,'B70': 0.0,'B71': 0.0,'B72': 0.0,'B73': 0.0,'B74': 0.0,'B75': 0.0,'B76': 0.0,'B77': 0.0}
    for key, val in Bkq.items():
        reg_key = 'B%d%d' % (key[0], key[1])
        req_key = reg_key.replace('-', 'm')
        default_Bkq[req_key] = val
    Bkq = default_Bkq
    re_and_im ={}
    for key, val in Bkq.items():
        re_and_im[key] = '%.1f %.1f' % (val.real, val.imag)
    file_content = cfp_template.format(**re_and_im)
    with open('./f%d/cfp.inp' % numE, 'w') as f:
        f.write(file_content)
    return file_content

def lanthanide(numE, F2=0., F4=0., F6=0., zeta=0.,
               alpha=0., beta=0., gamma=0.,
               T2=0., T3=0., T4=0., T6=0., T7=0., T8=0.,
               Bx=0., By=0., Bz=0.,
               Bkq={}, verbose=True):
    numStates = {1:14,
                 2:91,
                 3:364,
                 4:1001,
                 5:2002,
                 6:3003,
                 7:3432,
                 8:3003,
                 9:2002,
                 10:1001,
                 11:364,
                 12:91,
                 13:14}[numE]
    param_inp_file_maker(numE, F2=F2, F4=F4, F6=F6, zeta=zeta, alpha=alpha, beta=beta, gamma=gamma, T2=T2, T3=T3, T4=T4, T6=T6, T7=T7, T8=T8, Bx=Bx, By=By, Bz=Bz)
    cfp_inp_file_maker(numE, Bkq)
    if not os.path.exists('./lanthanide'):
        print("lanthanide binary not found, compiling ...")
        command = 'gcc -O3 -o ./lanthanide ./lanthanide.c  -llapack -lblas -lm -lgfortran -lz'
        execute_and_print_realtime(command)
    if not os.path.exists('./f%d/lanthanide' % numE):
        print("lanthanide binary not found in execution directory, copying ...")
        shutil.copy('./lanthanide', './f%d/lanthanide' % numE)
    if not os.path.exists('./f%d/vk1k2k3.gz' % numE):
        print("vk1k2k3.gz not found, running lanthanide 0 %d %d..." % (numStates, numE))
        command = ['./lanthanide','0','%d' % numStates, '%d' % numE]
        command = ' '.join(command)
        execute_and_print_realtime(command, workdir='./f%d' % numE, verbose=verbose)
    print("Running lanthanide 1 %d %d..." % (numStates, numE))
    command = ['./lanthanide','1', '%d' % numStates, '%d' % numE]
    command = ' '.join(command)
    execute_and_print_realtime(command, workdir='./f%d' % numE, verbose=verbose)
    energies = np.genfromtxt('./f%d/energies' % numE)
    return energies