#include <iostream>
#include <gflags/gflags.h>
#include "thresher/thresher.h"


using Thresher::ThreshCostFunction;
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

    /* MatrixXd scene = MatrixXd::Zero(dim, dim); */
    /* std::vector<MatrixXd> psfs; */
    /* std::vector<double> skies; */

    /* for (int ind = 0; ind < 2; ++ind) { */
    /*     MatrixXd img = readImage(fns[ind].c_str()); */

    /*     int h = img.rows(), w = img.cols(); */
    /*     img = img.block(int(h / 2 - img_hw), int(w / 2 - img_hw), 2 * img_hw + 1, 2 * img_hw + 1); */

    /*     h = img.rows(); */
    /*     w = img.cols(); */
    /*     double mu = img.sum() / h / w; */

    /*     for (int i = psf_hw; i < h + psf_hw; ++i) */
    /*         for (int j = psf_hw; j < w + psf_hw; ++j) */
    /*             scene(i, j) += img(i - psf_hw, j - psf_hw) - mu; */

    /*     MatrixXd psf = MatrixXd::Zero(2 * psf_hw + 1, 2 * psf_hw + 1); */
    /*     psf(psf_hw, psf_hw) = 1.0; */
    /*     psfs.push_back(psf); */

    /*     skies.push_back(0.0); */

    /*     p.AddResidualBlock(new ThreshCostFunction(img, psf_hw), NULL, &scene(0), &(psfs[ind](0)), &(skies[ind])); */
    /* } */

    /* scene /= 2; */

    /* Solver::Options options; */
    /* /1* options.check_gradients = true; *1/ */
    /* options.max_num_iterations = 100; */
    /* options.minimizer_progress_to_stdout = true; */

    /* Solver::Summary summary; */
    /* Solve(options, &p, &summary); */

    /* cout << summary.BriefReport() << endl; */
    /* /1* cout << scene << endl; *1/ */

    return 0;
}
