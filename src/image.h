#ifndef THRESHER_IMAGE_H
#define THRESHER_IMAGE_H


#include <iostream>

#include <Eigen/Dense>
#include <complex>
#include <fftw3.h>

#include "ceres/ceres.h"


using namespace Eigen;


namespace Thresher {

    class PSF {

        public:

            PSF(MatrixXd values)
                : halfwidth_((values.rows() - 1) / 2),
                  values_(values) {};

            int halfwidth() const { return halfwidth_; };
            MatrixXd values() const { return values_; };

            double operator () (int i, int j) const { return values_(i, j); };

        private:

            const int halfwidth_;
            const MatrixXd values_;

    };

    PSF gaussianPSF(int hw, double var) {
        int P = 2 * hw + 1;
        MatrixXd p = MatrixXd::Zero(P, P);
        for (int xi = 0; xi < P; ++xi)
            for (int yi = 0; yi < P; ++yi)
                p(xi, yi) = exp(-0.5 * ((xi - hw) * (xi - hw) + (yi - hw) * (yi - hw)) / var) / sqrt(2 * M_PI * var);
        PSF psf(p);
        return psf;
    };

    template <typename T>
    class Scene {

        public:

            Scene(int width, int height) {
                width_ = width;
                height_ = height;
                values_ = MatrixXd::Zero(width, height);
            };
            Scene(MatrixXd img) {
                width_ = img.rows();
                height_ = img.cols();
                values_ = img;
            };

            MatrixXd values() { return values_; };
            void values(MatrixXd newValues) { values_ = newValues; };

            double operator () (int i, int j) const { return values_(i, j); };
            void operator () (int i, int j, double val) { values_(i, j) = val; };

            Matrix<T, Dynamic, Dynamic> render(PSF psf) {
                Matrix<T, Dynamic, Dynamic> img = Matrix<T, Dynamic, Dynamic>::Zero(width_, height_);
                int hw = psf.halfwidth();
                for (int i = 0; i < width_; ++i)
                    for (int j = 0; j < height_; ++j)
                        for (int x = fmax(0, i - hw), xi = x - i + hw; x < fmin(width_, i + hw + 1); ++x, ++xi)
                            for (int y = fmax(0, j - hw), yi = y - j + hw; y < fmin(height_, j + hw + 1); ++y, ++yi)
                                img(i, j) += values_(x, y) * T(psf(xi, yi));
                return img;
            };

        private:

            int width_, height_;
            Matrix<double, Dynamic, Dynamic> values_;

    };

    /* template <typename T> */
    /* int fftconvolve(const double *a, int *da, const T *b, int *db, T *out) { */
    /*     fftw_complex *a2 = (fftw_complex*) fftw_malloc(da[0] * da[1] * sizeof(fftw_complex)); */
    /*     T *a0 = (T*) fftw_malloc(da[0] * da[1] * sizeof(T)); */
    /*     for (int i = 0; i < da[0] * da[1]; ++i) */
    /*         a0[i] = a[i]; */

    /*     fftw_plan p = fftw_plan_dft_r2c_2d(da[0], da[1], a0, a2, FFTW_ESTIMATE); */
    /*     /1* fftw_plan p2 = fftw_plan_dft_c2r_2d(da[0], da[1], a2, a, FFTW_ESTIMATE); *1/ */

    /*     fftw_execute(p); */

    /*     fftw_destroy_plan(p); */
    /*     fftw_free(a2); */

    /*     return 0; */
    /* }; */

    class ThreshCostFunction : public ceres::CostFunction {

        public:

            ThreshCostFunction(MatrixXd data, MatrixXd givenpsf)
                : data_(data),
                  psf_(givenpsf),
                  width_(data.rows()),
                  height_(data.cols())
            {
                set_num_residuals(width_ * height_);
                mutable_parameter_block_sizes()->push_back(width_ * height_);
            };

            bool Evaluate(double const* const* parameters,
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

        private:

            const int width_, height_;
            const MatrixXd data_;
            const MatrixXd psf_;

    };

};

#endif
