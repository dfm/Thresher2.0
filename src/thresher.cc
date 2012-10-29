#include "thresher/thresher.h"


using namespace Thresher;


PSF Thresher::gaussianPSF(int hw, double var) {
    int P = 2 * hw + 1;
    MatrixXd p = MatrixXd::Zero(P, P);
    for (int xi = 0; xi < P; ++xi)
        for (int yi = 0; yi < P; ++yi)
            p(xi, yi) = exp(-0.5 * ((xi - hw) * (xi - hw) + (yi - hw) * (yi - hw)) / var) / sqrt(2 * M_PI * var);
    PSF psf(p);
    return psf;
};


bool ThreshCostFunction::Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const
{

    const double *scene = parameters[0];

    int hw = (psf_.rows() - 1) / 2, dim = width_ * height_;
    int pd[] = {psf_.rows(), psf_.cols()},
        sd[] = {width_, height_};

    for (int i = 0, n = 0; i < width_; ++i)
        for (int j = 0; j < height_; ++j, ++n)
            residuals[n] = data_(i, j);

    for (int i = 0, n = 0; i < width_; ++i)
        for (int j = 0; j < height_; ++j, ++n)
            for (int x = fmax(0, i - hw), xi = x - i + hw; x < fmin(width_, i + hw + 1); ++x, ++xi)
                for (int y = fmax(0, j - hw), yi = y - j + hw; y < fmin(height_, j + hw + 1); ++y, ++yi)
                    residuals[n] -= scene[x * height_ + y] * psf_(xi, yi);

    if (jacobians != NULL && jacobians[0] != NULL) {
        for (int i = 0; i < width_ * height_; ++i)
            for (int j = 0; j < width_ * height_; ++j)
                jacobians[0][i * width_ * height_ + j] = 0.0;

        for (int i = 0, n = 0; i < width_; ++i)
            for (int j = 0; j < height_; ++j, ++n)
                for (int x = fmax(0, i - hw), xi = x - i + hw; x < fmin(width_, i + hw + 1); ++x, ++xi)
                    for (int y = fmax(0, j - hw), yi = y - j + hw; y < fmin(height_, j + hw + 1); ++y, ++yi)
                        jacobians[0][n * width_ * height_ + x * height_ + y] = -psf_(xi, yi);
    }

    return true;

};
