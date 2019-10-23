import numpy as np
from mstools.panedr import edr_to_df

def get_acf(x_list, y_list):
    '''
    x_list must be evenly spaced
    '''
    n = len(x_list)
    if len(y_list) != n:
        return 'x_list, y_list have different length'
    y_list = np.array(y_list)
    y_list -= y_list.mean()
    dx = x_list[1] - x_list[0]
    _x_list = []
    _acf_list = []
    for i in range(int(n/2)):
        Dx = i * dx
        acf = 0
        for j in range(n-i):
            acf += y_list[j] * y_list[j+i]
        acf /= (n-i)
        _x_list.append(Dx)
        _acf_list.append(acf)
    return _x_list, _acf_list

def get_integral_acf(x_list, y_list):
    _x_list, _acf_list = get_acf(x_list, y_list)
    dx = _x_list[1] - _x_list[0]
    acf_integral = _acf_list[0] * dx * 0.5
    acf_integral_list = [acf_integral]
    for i in range(1, len(_x_list)):
        acf_integral += _acf_list[i] * dx
        acf_integral_list.append(acf_integral)
    return (np.array(_x_list) + 0.5 * dx), np.array(acf_integral_list)


