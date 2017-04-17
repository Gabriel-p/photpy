

def main():
    """
    """

    # Read .coo file with standard stars coordinates in the system of the
    # observed frame.
    read_standard_coo()

    # For each .fits standard file.
    for imname in fits_list:

        # Load .fits file.
        hdulist = fits.open(imname)
        # Extract header and data.
        hdr, hdu_data = hdulist[0].header, hdulist[0].data

        # Background subtraction.
        # aperture_photometry() assumes that the data have been
        # background-subtracted




if __name__ == '__main__':
    main()
