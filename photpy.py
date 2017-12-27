
import argparse


def main():
    """
    Call as:

    $ python photpy.py name_of_task

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
    elif task == 'aperphot_standards':
        from tasks import aperphot_standards
        aperphot_standards.main()
    elif task == 'fit_standard':
        from tasks import fit_standard
        fit_standard.main()
    elif task == 'find_stars':
        from tasks import find_stars
        find_stars.main()
    elif task == 'psf_phot':
        from tasks import psf_phot
        psf_phot.main()
    elif task == 'match':
        from tasks import match
        match.main()
    elif task == 'transf':
        from tasks import transf
        transf.main()
    else:
        print("Unrecognized task.")


if __name__ == '__main__':
    main()
