#include <iostream>
#include <gflags/gflags.h>
#include "thresher/thresher.h"


using Thresher::readImage;
using Thresher::writeImage;
using Thresher::appendLuckyMetadata;

using std::cout;
using std::cin;
using std::getline;
using std::endl;
using std::vector;
using std::string;


int main(int argc, char **argv)
{
    // Parse the command line arguments.
    google::ParseCommandLineFlags(&argc, &argv, true);

    // Find input files.
    vector<string> fns;
    string fn;

    cout << "Enter data filenames:" << endl;
    while (getline(cin, fn))
        fns.push_back(fn);

    // The number of images.
    int N = fns.size();

    // Load the images into RAM and compute the centers.
    vector <MatrixXd> imgs;
    VectorXd values(N);
    MatrixXi dim_min(N, 2),
             dim_max(N, 2);
    for (int i = 0; i < N; ++i) {
        // Load the image.
        imgs.push_back(readImage(fns[i].c_str()));

        // Find the value and coordinates of the maximum pixel.
        int mx, my;
        values(i) = imgs[i].maxCoeff(&mx, &my);

        // Compute the bounds of the image.
        float h = 0.5 * imgs[i].rows(),
              w = 0.5 * imgs[i].cols();
        dim_min(i, 0) = mx - floor(h);
        dim_min(i, 1) = my - floor(w);
        dim_max(i, 0) = mx + ceil(h);
        dim_max(i, 1) = my + ceil(w);
    }

    // Compute the maximum and minimum offsets.
    VectorXi mn = dim_min.colwise().minCoeff(),
             mx = dim_max.colwise().maxCoeff(),
             s = mx - mn;

    // Offset the dimension vector.
    dim_min.rowwise() -= mn.transpose();
    dim_max.rowwise() -= mn.transpose();

    // Compute the co-add.
    MatrixXd result = MatrixXd::Zero(s(0), s(1));
    MatrixXi counts = MatrixXi::Zero(s(0), s(1));
    for (int i = 0; i < N; ++i) {
        counts.block(dim_min(i, 0), dim_min(i, 1),
                     imgs[i].rows(), imgs[i].cols())
                            += MatrixXi::Ones(imgs[i].rows(), imgs[i].cols());
        result.block(dim_min(i, 0), dim_min(i, 1),
                     imgs[i].rows(), imgs[i].cols()) += imgs[i];
    }

    for (int i = 0, l1 = result.rows(); i < l1; ++i)
        for (int j = 0, l2 = result.cols(); j < l2; ++j)
            if (counts(i, j) > 0)
                result(i, j) /= counts(i, j);

    writeImage("lucky.fits", result);
    appendLuckyMetadata("lucky.fits", dim_min, dim_max, values, fns);
    return 0;
}
