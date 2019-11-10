import matplotlib.pyplot as plt

def plot(x_list, y_list):
    plt.plot(x_list, y_list)
    plt.show()

def gnuplot(output, x_min, x_max, y_min, y_max, xlabel, ylabel, title, txt_list=[], type_list=[], title_list=[], x_log_scale=False, y_log_scale=False):
    f = open('1.gpi', 'w')
    info = 'set terminal pngcairo size 1000,800 enhanced font \'Times New Roman,25\'\n'
    info += 'set output "%s"\n' % (output)
    info += 'set title "%s"\n' % (title)
    info += 'set border lw 1.5\n'

    info += 'set xlabel "%s"\n' % (xlabel)
    info += 'set ylabel "%s"\n' % (ylabel)

    info += 'set xr [%f, %f]\n' % (x_min, x_max)
    info += 'set yr [%f, %f]\n' % (y_min, y_max)

    info += 'set style line 1 lt 1 lw 3 lc rgb "000000" # black\n'
    info += 'set style line 2 lt 1 lw 3 lc rgb "ff0000" # red\n'
    info += 'set style line 3 lt 1 lw 3 lc rgb "0000ff" # blue\n'
    info += 'set style line 4 lt 1 lw 3 lc rgb "00eeee" # dark-cyan\n'

    if x_log_scale:
        info += 'set logscale x\n'
    if y_log_scale:
        info += 'set logscale y\n'

    info += 'pl '
    for i, txt in enumerate(txt_list):
        if type_list[i] == 'errorbars':
            info += '"%s" u 1:2:3 with errorbars ls %i title "%s"' % (txt, i+1, title_list[i])
        if i + 1 != len(txt_list):
            info += ', \\\n'

