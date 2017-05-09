
import argparse


def main():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('task', type=str, help='Name of task to execute')

    args = parser.parse_args()
    task = args.task
    print("\nCalling task: {}".format(task))

    if task == 'fitstats':
        from tasks import fitstats
        fitstats.main()
    elif task == 'align_crop':
        from tasks import align_crop
        align_crop.main()
    elif task == 'id_standard':
        from tasks import id_standard
        id_standard.main()
    else:
        print("Unrecognized task.")


if __name__ == '__main__':
    main()
