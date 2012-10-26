#ifndef THRESHER_IMAGE_H
#define THRESHER_IMAGE_H


#include <Eigen/Dense>


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

    class ThreshCostFunction {

        public:

            ThreshCostFunction(MatrixXd data, MatrixXd givenpsf)
                : data_(data),
                  psf_(givenpsf),
                  width_(data.rows()),
                  height_(data.cols()) {};

            template <typename T> bool operator () (const T* const scene, T* e) const {

                int hw = (psf_.rows() - 1) / 2, dim = width_ * height_;

                for (int i = 0, n = 0; i < width_; ++i)
                    for (int j = 0; j < height_; ++j, ++n)
                        e[n] = T(data_(i, j));

                for (int i = 0, n = 0; i < width_; ++i)
                    for (int j = 0; j < height_; ++j, ++n)
                        for (int x = fmax(0, i - hw), xi = x - i + hw; x < fmin(width_, i + hw + 1); ++x, ++xi)
                            for (int y = fmax(0, j - hw), yi = y - j + hw; y < fmin(height_, j + hw + 1); ++y, ++yi)
                                e[n] -= scene[x * height_ + y] * T(psf_(xi, yi));

                return true;

            };

        private:

            const int width_, height_;
            const MatrixXd data_;
            const MatrixXd psf_;

    };

};

#endif
