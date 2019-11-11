import matplotlib.pyplot as plt
import subprocess
from subprocess import Popen, PIPE

def plot(x_list, y_list):
    plt.plot(x_list, y_list)
    plt.show()

def gnuplot(output, xlabel, ylabel, title, txt_list=[], type_list=[], title_list=[], x_min=None, x_max=None, y_min=None, y_max=None, x_log_scale=False, y_log_scale=False):

    f = open('%s.gpi' % (output), 'w')
    info = 'set terminal pngcairo size 1000,800 enhanced font \'Times New Roman,25\'\n'
    info += 'set output "%s.png"\n' % (output)
    info += 'set title "%s"\n' % (title)
    info += 'set border lw 1.5\n'

    info += 'set xlabel "%s"\n' % (xlabel)
    info += 'set ylabel "%s"\n' % (ylabel)

    if x_min is not None:
        info += 'set xr [%f:%f]\n' % (x_min, x_max)
    if y_min is not None:
        info += 'set yr [%f:%f]\n' % (y_min, y_max)

    info += 'set style line 1 lt 1 lw 1 lc rgb "#000000" # black\n'
    info += 'set style line 2 lt 1 lw 1 lc rgb "#ff0000" # red\n'
    info += 'set style line 3 lt 1 lw 1 lc rgb "#0000ff" # blue\n'
    info += 'set style line 4 lt 1 lw 1 lc rgb "#00eeee" # dark-cyan\n'

    if x_log_scale:
        info += 'set logscale x\n'
    if y_log_scale:
        info += 'set logscale y\n'

    info += 'pl '
    for i, txt in enumerate(txt_list):
        if type_list[i] == 'errorbars':
            info += '"%s" u 1:2:3 with errorbars ls %i title "%s"' % (txt, i + 1, title_list[i])
        elif type_list[i] == 'errorlines':
            info += '"%s" u 1:2:3 with errorlines ls %i title "%s"' % (txt, i + 1, title_list[i])
        elif type_list[i] == 'lines':
            info += '"%s" u 1:2 with lines ls %i title "%s"' % (txt, i + 1, title_list[i])
        if i + 1 != len(txt_list):
            info += ', \\\n'
    f.write(info)

    cmd = 'gnuplot < %s.gpi' % (output)
    sp = Popen(cmd.split(), stdin=PIPE)
    sp.communicate()