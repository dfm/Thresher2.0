#include <iostream>
#include <gflags/gflags.h>
#include "thresher/thresher.h"


using Thresher::readImage;

using std::cout;
using std::cin;
using std::getline;
using std::endl;
using std::vector;
using std::string;

using ceres::Problem;
using ceres::AutoDiffCostFunction;
using ceres::Solver;
using ceres::Solve;


// Command line options.
DEFINE_string(i, "",
              "The path to the file containing an initial guess at the scene.");


int main(int argc, char **argv)
{
    int psf_hw = 3, img_hw = 10, dim = 2 * img_hw + 1 + 2 * psf_hw;

    // Parse the command line arguments.
    google::ParseCommandLineFlags(&argc, &argv, true);
    if (FLAGS_i == "") {
        cout << "You must provide an initial file for now." << endl;
        return -1;
    }

    // Find input files.
    vector<string> fns;
    string fn;

    cout << "Enter data filenames:" << endl;
    while (getline(cin, fn))
        fns.push_back(fn);

    for (int i = 0, l = fns.size(); i < l; ++i) {
        MatrixXd img = readImage(fns[i].c_str());
        cout << fns[i] << ": " << img.rows() << " " << img.cols() << endl;
    }

    Problem p;

    return 0;
}
