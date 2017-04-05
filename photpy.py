
import argparse

from tasks import getdata
from tasks import id_standard


def main():
    """
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('task', type=str, help='Name of task to execute')

    args = parser.parse_args()
    task = args.task
    print("\nCalling task: {}".format(task))

    if task == 'getdata':
        getdata.main()
    elif task == 'id_standard':
        id_standard.main()
    else:
        print("Unrecognized task.")


if __name__ == '__main__':
    main()
