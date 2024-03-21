#!/usr/bin/env python3

import sys
import os
from math import ceil

keys = ['n', 'c', 'g', 'f', 'eps', 's']

num_args: int = len(sys.argv)

if num_args <= 7:
    instance_name: str = ''
    for i in range(num_args - 1):
        instance_name += ('' if i == 0 else '_') + keys[i] + '_' + str(sys.argv[i + 1])
    for i in range(num_args - 1, 6):
        value: str = str(input('Enter the value for ' + keys[i] + ': '))
        instance_name += ('' if i == 0 else '_') + keys[i] + '_' + value

    num_bins: int = int(input('Enter the number of bins: '))

elif num_args == 8:
    instance_name: str = ''
    for i in range(6):
        instance_name += ('' if i == 0 else '_') + keys[i] + '_' + str(sys.argv[i + 1])

    num_bins: int = int(sys.argv[7])
    if num_bins <= 0:
        print('Number of bins must be positive.')
        sys.exit(1)

else:
    print('Too many arguments passed => aborting.')
    sys.exit(1)

script_dir: str = os.path.dirname(os.path.realpath(sys.argv[0]))
instance_dir: str = os.path.join(os.path.dirname(script_dir), 'instances', instance_name)
path_to_instance_data: str = os.path.join(instance_dir, 'results.txt')
if not os.path.isfile(path_to_instance_data):
    print('For the given instance no data is available.')
    sys.exit(1)

if num_bins <= 0:
    print('Number of bins must be positive.')
    sys.exit(1)

os.chdir(instance_dir)

with open('results.txt', 'r') as file:
    num_states: int = int(file.readline()) # read in number of recorded states
    next(file) # ignore optimal solution value
    next(file) # ignore final approximation ratio
    data: list[tuple[float, float]] = []
    for i in range(num_states):
        line: list[str] = next(file).split()
        data.append((float(line[0]), float(line[1])))

# create list of bins with their total probability
plot_data: list[float] = [0.0 for _ in range(num_bins)]
for s in range(num_states):
    for b in range(num_bins):
        if b * 1. / num_bins <= data[s][0] <= (b + 1) * 1. / num_bins:
            plot_data[b] += data[s][1]

# create LaTeX file for plotting
latex_header: str = ('\\documentclass{article}\n'
                     + '\\usepackage{pgfplots}\n'
                     + '\\pgfplotsset{compat = newest}\n'
                     + '\\usepackage{tikz}\n'
                     + '\\usepackage{color}\n'
                     + '\\definecolor{QCPbeige}{HTML}{E7D5AB} %231 213 171\n'
                     + '\\definecolor{QCPblue}{HTML}{479093} %71 144 147\n'
                     + '\\definecolor{QCPdarkblue}{HTML}{2C6071} %44  96 113\n'
                     + '\\definecolor{QCPgray}{HTML}{333743} %51  55  67\n'
                     + '\\pgfrealjobname{MakeProbabilityPlot}\n')

latex_main_body: str = ('\\begin{document}\n'
                        + '\\beginpgfgraphicnamed{'
                        + instance_name
                        + '_plot}\n'
                        + '\\begin{tikzpicture}\n'
                        + '\t\\begin{axis}[\n'
                        + '\t\tname=theaxis,\n'
                        + '\t\tscale only axis,\n'
                        + '\t\tybar interval,\n'
                        + '\t\tenlargelimits=0.0,\n'
                        + '\t\tminor y tick num = 4,\n'
                        + '\t\tymin=0,\n'
                        + f'\t\tymax={min(0.05 * (ceil(max(plot_data) / 0.05) + 1), 1.0)},\n'
                        + '\t\theight=0.4\\textwidth,\n'
                        + '\t\tgrid,\n'
                        + '\t\tgrid style=dashed,\n'
                        + '\t\txtick style={draw=none},\n'
                        + '\t\txminorgrids=false,\n'
                        + '\t\txmajorgrids=false,\n'
                        + '\t\tyminorgrids=true,\n'
                        + '\t\tymajorgrids=true,\n'
                        + '\t\txlabel=Approximation ratio (lower bound),\n'
                        + '\t\tylabel=Sample probability,\n'
                        + '\t\txtick=data,\n'
                        + '\t\txticklabels = {0.0')
for b in range(1, num_bins):
    latex_main_body += f', {b / num_bins}'
latex_main_body += ('},\n'
                    + '\t\txticklabel style={\n'
                    + '\t\t\tyshift=0.25cm,\n'
                    + '\t\t\trotate=75,\n'
                    + '\t\t\tfont=\\scriptsize,\n'
                    + '\t\t}\n'
                    + '\t]\n'
                    + '\t\\addplot [draw=QCPdarkblue, fill=QCPblue] coordinates {')
for b in range(num_bins):
    latex_main_body += f' ({b / num_bins}, {plot_data[b]})'
latex_main_body += (' (1.0, 0.0)};\n'
                    + '\t\\end{axis}\n'
                    + '\\end{tikzpicture}\n'
                    + '\\endpgfgraphicnamed\n'
                    + '\\end{document}\n')

with open('tmp.tex', 'w') as latex_file:
    latex_file.write(latex_header)
    latex_file.write(latex_main_body)

os.system(f'pdflatex -shell-escape -jobname={instance_name}_plot tmp.tex')
os.remove('tmp.tex')
os.remove(instance_name + '_plot.log')
