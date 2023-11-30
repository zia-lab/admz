#!/usr/bin/env python3

import numpy as np
import h5py
import os

def save_to_h5(data, filename, group=None, comments='', overwrite=False):
    '''
    Save numeric values, numpy arrays, strings, and nested dictionaries
    from the provided dictionary to an H5 file.
    If a value of the provided dictionary is a list, then
    the list is converted to a dictionary with keys equal
    to the indices of the lists (as strings).    

    Parameters
    ----------
    data: dict
        The input dictionary.
    filename: str
        The path to the H5 file to save.
    group: (h5py.Group, optional)
        Current group for recursion. Defaults to None.
    
    Returns
    -------
    None

    Example
    -------
    data_dict = {
        'a': [1,2,3],
        'b': {
            'c': "hello",
            'd': np.array([1, 2, 3]),
            'e': {
                'f': 3.14,
                'g': "world"
            }
        },
        'h': np.array([[1, 2], [3, 4]])
    }
    save_to_h5(data_dict, 'output.h5')
    '''
    if group == None and os.path.exists(filename):
        if overwrite:
            print("File already exists, overwriting ...")
        else:
            print("File already exists, doing nothing ...")
            return None
    
    # Initialize the file/group on the first call
    if group is None:
        with h5py.File(filename, 'w') as h5f:
            if comments:
                h5f.attrs['comments'] = comments.encode('utf-8')
            save_to_h5(data, filename, group=h5f)
    else:
        for key, value in data.items():
            # If the value is a dictionary, recurse into it
            if isinstance(value, dict):
                if isinstance(key, int):
                    key = str(key)
                subgroup = group.create_group(key)
                save_to_h5(value, filename, group=subgroup)
                
            # Convert lists to dictionaries
            elif isinstance(value, list):
                list_dict = {str(i): v for i, v in enumerate(value)}
                subgroup = group.create_group(key)
                save_to_h5(list_dict, filename, group=subgroup)
            
            else:
                # If the value is a string, encode it to UTF-8
                if isinstance(value, str):
                    value = value.encode('utf-8')
                group.create_dataset(key, data=value)
