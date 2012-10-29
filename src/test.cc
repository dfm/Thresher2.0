#include <iostream>
#include "image.h"


using namespace Thresher;


using std::cout;
using std::endl;

using ceres::Problem;
using ceres::AutoDiffCostFunction;
using ceres::Solver;
using ceres::Solve;


int main()
{
    const int width = 11, height = 11;
    const int dim = width * height;

    PSF psf = gaussianPSF(6, 2);

    Scene<double> scene(width, height);
    scene(2, 5, 1.0);
    scene(6, 6, 2.0);
    MatrixXd img = scene.render(psf);

    MatrixXd p0 = MatrixXd::Zero(img.rows(), img.cols());

    Problem p;
    p.AddResidualBlock(new ThreshCostFunction(img, psf.values()), NULL,
                        &p0(0));

    Solver::Options options;
    options.max_num_iterations = 25;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = false;  // true;

    Solver::Summary summary;
    Solve(options, &p, &summary);

    /* cout << summary.BriefReport() << endl; */
    cout << p0 << endl;

    return 0;
}
