import numpy as np
import math, os
import pandas as pd

def get_acf(x_list, y_list, mean_shift=False):
    '''
    x_list must be evenly spaced
    '''
    n = len(x_list)
    if len(y_list) != n:
        return 'x_list, y_list have different length'
    y_list = np.array(y_list)
    if mean_shift:
        y_list -= y_list.mean()
    dx = x_list[1] - x_list[0]
    _x_list = []
    _acf_list = []
    for i in range(int(n/2)):
        Dx = i * dx
        acf = 0.
        for j in range(n-i):
            acf += y_list[j] * y_list[j+i]
        acf /= (n-i)
        _x_list.append(Dx)
        _acf_list.append(acf)
    return np.array(_x_list), np.array(_acf_list)

def get_integral(x_list, y_list):
    dx = x_list[1] - x_list[0]
    acf_integral = y_list[0] * dx * 0.5
    acf_integral_list = [acf_integral]
    for i in range(1, len(x_list)):
        acf_integral += y_list[i] * dx
        acf_integral_list.append(acf_integral)
    return (np.array(x_list) + 0.5 * dx), np.array(acf_integral_list)

def get_integral_acf(x_list, y_list):
    _x_list, _acf_list = get_acf(x_list, y_list)
    return get_integral(_x_list, _acf_list)

def get_block_average(list, n_block=None):
    if n_block is None:
        return list
    list = np.array(list)
    list_out = []
    for i in range(math.floor(list.size / n_block)):
        list_out.append(list[i * n_block:(i + 1) * n_block].mean())
    return np.array(list_out)

def get_t_property_list(property, dir=None, name=None, weight=0.00):
    Property_dict = {
        'viscosity': {'abbr': 'vis', 'property_unit': 'viscosity(mPa·s)'},
        'electrical conductivity': {'abbr': 'econ', 'property_unit': 'electrical_conductivity(S/m)'},
        'diffusion constant': {'abbr': 'diff', 'property_unit': 'diffusion_constant(cm^2/s)'},
    }
    if Property_dict.get(property) is None:
        raise Exception('invalid property: %s', property)
    file_name = Property_dict.get(property).get('abbr')
    if name is not None:
        file_name += '-%s' % (name)
    if weight != 0.00:
        file_name += '-%.2f' % (weight)
    file_name += '.txt'
    if dir is not None:
        file_name = os.path.join(dir, file_name)
    info = pd.read_csv(file_name, sep='\s+', header=0)
    t_list = np.array(info['#time(ps)'])
    property_list = np.array(info[Property_dict.get(property).get('property_unit')])
    return t_list, property_list

def get_std_out(input):
    if type(input)==float:
        return float('%#.5e' % (input))
    if type(input)==list:
        output = []
        for i in input:
            output.append(float('%#.5e' % (i)))
        return output