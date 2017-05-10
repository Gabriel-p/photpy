
from os import getcwd
from os.path import join, realpath, dirname, isfile
import sys


def main():
    """
    Read parameter values from 'params_input.dat' file.
    """
    mypath = realpath(join(getcwd(), dirname(__file__)))
    pars_f = join(mypath.replace('tasks', ''), 'params_input.dat')
    if not isfile(pars_f):
        print("Parameters file is missing. Exit.")
        sys.exit()

    pars = {'mypath': mypath}
    with open(pars_f, 'r') as f:
        for line in f:
            if not line.startswith('#') and line != '\n':
                lin = line.replace('\n', '').split()
                key = lin[0]
                value = [_.replace(',', '') for _ in lin[1:]]
                if len(value) > 1:
                    # If the key exists, append instead of re-writing.
                    if key in pars:
                        pars[key].append(value)
                    else:
                        pars[key] = [value]
                else:
                    pars[key] = value[0]

    return pars


if __name__ == '__main__':
    main()
