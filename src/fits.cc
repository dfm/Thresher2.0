#include "thresher/thresher.h"
#include <fitsio.h>

using std::vector;
using std::string;

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


void Thresher::writeImage(const char *fn, MatrixXd data)
{
    // The file pointer.
    fitsfile *fptr;

    // Parameters of the image.
    int bitpix = DOUBLE_IMG;
    long naxis = 2;
    vector<long> naxes;
    naxes.push_back(data.rows());
    naxes.push_back(data.cols());

    // Clobber any existing file.
    remove(fn);

    // Create the FITS file.
    int status = 0;
    if (fits_create_file(&fptr, fn, &status))
        printerror(status);

    // Create the image header.
    if (fits_create_img(fptr, bitpix, naxis, &naxes[0], &status))
         printerror(status);

    // Write the image to the primary header.
    long fpixel = 1,
         nelements = naxes[0] * naxes[1];
    if (fits_write_img(fptr, TDOUBLE, fpixel, nelements, &data(0), &status))
        printerror(status);

    if (fits_close_file(fptr, &status))
        printerror(status);
}

void Thresher::appendLuckyMetadata(const char *fn,
                                   MatrixXi dim_min, MatrixXi dim_max,
                                   VectorXd ranks, vector<string> fns)
{
    // The file pointer.
    fitsfile *fptr;

    // Open the FITS file.
    int status = 0;
    if (fits_open_file(&fptr, fn, READWRITE, &status)) printerror(status);

    // Move to the second HDU to write the table of metadata.
    int hdutype;
    if (fits_movabs_hdu(fptr, 1, &hdutype, &status)) printerror(status);

    // Set up the columns.
    char *extname = (char*)"METATDATA_Binary";
    char *colname[] = {(char*)"xmin", (char*)"ymin", (char*)"xmax",
                       (char*)"ymax", (char*)"rank", (char*)"filename"};
    char *colfmt[] = {(char*)"1J", (char*)"1J", (char*)"1J", (char*)"1J",
                      (char*)"1E", (char*)"1024A"};
    char *colunit[] = {(char*)"px", (char*)"px", (char*)"px", (char*)"px",
                       (char*)"\0", (char*)"\0"};

    // Append an empty binary table.
    int nrows = dim_min.rows();
    if (fits_create_tbl(fptr, BINARY_TBL, nrows, 6,
                        colname, colfmt, colunit, extname, &status))
         printerror(status);

    // Write the coordinate columns.
    for (int i = 0; i < 2; ++i)
        fits_write_col(fptr, TINT, i + 1, 1, 1, nrows,
                       &dim_min.col(i)(0), &status);
    for (int i = 0; i < 2; ++i)
        fits_write_col(fptr, TINT, i + 3, 1, 1, nrows,
                       &dim_max.col(i)(0), &status);

    // Write the rank and filename columns.
    fits_write_col(fptr, TDOUBLE, 5, 1, 1, nrows, &ranks(0), &status);
    fits_write_col(fptr, TSTRING, 6, 1, 1, nrows, &fns[0], &status);

    // Commit the changes.
    if (fits_close_file(fptr, &status))
        printerror(status);
}
