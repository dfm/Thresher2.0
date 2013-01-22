#ifndef THRESHER_H
#define THRESHER_H


#include <iostream>
#include <cmath>

#include <Eigen/Dense>
#include <complex>
#include <fftw3.h>

#include "ceres/ceres.h"


using namespace Eigen;
using std::cout;


namespace Thresher {

    // =====================================================================
    //                                                                ERRORS
    // =====================================================================

    const int DimensionMismatch = -1;
    const int FITSError = -2;
    const int ConvolutionError = -3;


    // =====================================================================
    //                                             DATA MANIPULATION HELPERS
    // =====================================================================

    MatrixXd optimalPSF(MatrixXd scene, MatrixXd data, int psf_hw);
    MatrixXd fftconvolve(MatrixXd data, MatrixXd kernel);


    // =====================================================================
    //                                             CERES-STYLE COST FUNCTION
    // =====================================================================

    class ThreshCostFunction : public ceres::CostFunction {

        public:

            ThreshCostFunction(MatrixXd data, int psf_hw)
                : data_(data),
                  psf_hw_(psf_hw),
                  width_(data.rows()),
                  height_(data.cols())
            {
                set_num_residuals(width_ * height_);

                // Scene parameter block.
                mutable_parameter_block_sizes()->push_back((width_ + 2 * psf_hw) * (height_ + 2 * psf_hw));

                // PSF parameter block.
                int p = 2 * psf_hw + 1;
                mutable_parameter_block_sizes()->push_back(p * p);

                // Sky parameter block.
                mutable_parameter_block_sizes()->push_back(1);
            };

            bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const;

        private:

            const int width_, height_, psf_hw_;
            const MatrixXd data_;

    };


    // =====================================================================
    //                                                        FITSIO HELPERS
    // =====================================================================

    MatrixXd readImage(const char *fn);

};

#endif // THRESHER_H
