#include <iostream>
#include "image.h"
#include "ceres/ceres.h"


using namespace Thresher;


int main()
{
    PSF psf = gaussianPSF(6, 2);

    Scene<double> scene(11, 11);
    scene(2, 5, 1.0);
    scene(6, 6, 2.0);
    MatrixXd img = scene.render(psf);

    MatrixXd p0 = MatrixXd::Zero(img.rows(), img.cols());

    int dim = img.rows() * img.cols();

    ceres::Problem p;
    p.AddResidualBlock(
            new ceres::AutoDiffCostFunction<ThreshCostFunction, 121, 121>(
                                new ThreshCostFunction(img, psf.values())), NULL,
                                &p0(0));

    ceres::Solver::Options options;
    options.max_num_iterations = 25;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = false;  // true;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &p, &summary);
    /* std::cout << summary.BriefReport() << "\n"; */
    std::cout << p0 << std::endl;

    return 0;
}
