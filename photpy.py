
import argparse

from tasks import getdata


def main():
    """
    """

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('task', type=str, help='task to execute')

    args = parser.parse_args()
    task = args.task
    print("\nCalling task: {}".format(task))

    if task == 'getdata':
        getdata.main()


if __name__ == '__main__':
    main()
