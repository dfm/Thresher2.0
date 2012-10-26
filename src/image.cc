#include <iostream>
#include "image.h"


using namespace Thresher;


PSF gaussianPSF(int hw, double var) {
    int P = 2 * hw + 1;
    MatrixXd p = MatrixXd::Zero(P, P);
    for (int xi = 0; xi < P; ++xi)
        for (int yi = 0; yi < P; ++yi)
            p(xi, yi) = exp(-0.5 * ((xi - hw) * (xi - hw) + (yi - hw) * (yi - hw)) / var) / sqrt(2 * M_PI * var);
    PSF psf(p);
    return psf;
};


template <class T> Matrix<T, Dynamic, Dynamic> Scene::render(PSF psf)
{
    Matrix<T, Dynamic, Dynamic> img = Matrix<T, Dynamic, Dynamic>::Zero(width_, height_);
    int hw = psf.halfwidth();

    for (int i = 0; i < width_; ++i)
        for (int j = 0; j < height_; ++j)
            for (int x = fmax(0, i - hw), xi = x - i + hw; x < fmin(width_, i + hw + 1); ++x, ++xi)
                for (int y = fmax(0, j - hw), yi = y - j + hw; y < fmin(height_, j + hw + 1); ++y, ++yi)
                    img(i, j) += T(values_(x, y) * psf(xi, yi));

    return img;
};
