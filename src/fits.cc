#include "thresher/thresher.h"
#include <fitsio.h>


void printerror(int status)
{
    if (status) {
        fits_report_error(stderr, status);
        throw Thresher::FITSError;
    }
}


MatrixXd Thresher::readImage(const char *fn)
{
    fitsfile *f;
    int status = 0;

    // Open the file.
    if (fits_open_file(&f, fn, READONLY, &status)) printerror(status);

    // Get the dimensions.
    int nfound;
    long naxes[2];
    if (fits_read_keys_lng(f, "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);

    // Allocate the image.
    MatrixXd img(naxes[0], naxes[1]);

    // Read in the data.
    int anynull;
    long fpixel = 1;
    double nullval = 0;
    if (fits_read_img(f, TDOUBLE, fpixel, naxes[0] * naxes[1], &nullval,
                      &img(0), &anynull, &status))
        printerror(status);

    if (fits_close_file(f, &status))
        printerror(status);

    return img;
}
