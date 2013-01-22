#include "thresher/thresher.h"


using std::cout;
using std::endl;


MatrixXd Thresher::optimalPSF(MatrixXd scene, MatrixXd data, int psf_hw)
{
    // Compute and check the dimensions.
    int dh = data.rows(), dw = data.cols(), D = dh * dw,
        sh = scene.rows(), sw = scene.cols(), S = sh * sw,
        ph = 2 * psf_hw + 1, P = ph * ph;

    if (dh + 2 * psf_hw != sh || dw + 2 * psf_hw != sw)
        throw DimensionMismatch;

    // First, construct the scene matrix.
    MatrixXd sm(D, P + 1);

    // Loop over the scene and populate the matrix.
    for (int sx = psf_hw, d = 0; sx < sh - psf_hw; ++sx)
        for (int sy = psf_hw; sy < sw - psf_hw; ++sy, ++d) {

            for (int px = -psf_hw, p = 0; px <= psf_hw; ++px)
                for (int py = -psf_hw; py <= psf_hw; ++py, ++p)
                    sm(d, p) = scene(sx + px, sy + py);

            // Add in the sky component.
            sm(d, P) = 1.0;
        }

    // Then build the flattened data vector.
    VectorXd dv(D);
    for (int i = 0, n = 0; i < dh; ++i)
        for (int j = 0; j < dw; ++j, ++n)
            dv(n) = data(i, j);

    // Solve the system.
    VectorXd psf = sm.jacobiSvd(ComputeThinU | ComputeThinV).solve(dv);

    cout << "The sky value is: " << psf(P) << std::endl;

    // Reshape the PSF for output.
    MatrixXd result(ph, ph);
    for (int i = 0, n = 0; i < ph; ++i)
        for (int j = 0; j < ph; ++j, ++n)
            result(i, j) = psf(n);

    // Compute the residual.
    VectorXd diff = dv - sm * psf;
    cout << "The error is: " <<  diff.norm() / dv.norm() << "\n";

    return result;
}

bool Thresher::ThreshCostFunction::Evaluate(double const* const* parameters,
                                  double* residuals,
                                  double** jacobians) const
{

    int hw = psf_hw_, p = 2 * hw + 1,
        w = width_ + 2 * hw, h = height_ + 2 * hw,
        dim = w * h;

    MatrixXd scene(h, w);
    for (int i = 0, n = 0; i < h; ++i)
        for (int j = 0; j < w; ++j, ++n)
            scene(i, j) = parameters[0][n];

    MatrixXd psf(p, p);
    for (int i = 0, n = 0; i < p; ++i)
        for (int j = 0; j < p; ++j, ++n)
            psf(i, j) = parameters[1][n];

    double sky = parameters[2][0];

    // Compute the residuals.
    for (int i = hw, n = 0; i < height_ + hw; ++i)
        for (int j = hw; j < width_ + hw; ++j, ++n) {
            residuals[n] = data_(i - hw, j - hw);
            for (int k = -hw; k <= hw; ++k)
                for (int l = -hw; l <= hw; ++l)
                    residuals[n] -= (psf(k + hw, l + hw) + sky) * scene(i + k, j + l);
        }

    // Compute the Jacobian if requested.
    if (jacobians != NULL && jacobians[0] != NULL) {
        // Scene gradients.
        for (int i = hw, ind = 0; i < height_ + hw; ++i)
            for (int j = hw; j < width_ + hw; ++j)
                for (int n = 0; n < h; ++n)
                    for (int m = 0; m < w; ++m, ++ind)
                        if (n - i >= -hw && n - i <= hw && m - j >= -hw && m - j <= hw)
                            jacobians[0][ind] = -psf(n - i + hw, m - j + hw) - sky;
                        else jacobians[0][ind] = 0.0;

        // PSF gradients.
        for (int i = hw, ind = 0; i < height_ + hw; ++i)
            for (int j = hw; j < width_ + hw; ++j)
                for (int n = -hw; n <= hw; ++n)
                    for (int m = -hw; m <= hw; ++m, ++ind)
                        jacobians[1][ind] = -scene(i + n, j + m);

        // Sky gradient.
        for (int i = hw, ind = 0; i < height_ + hw; ++i)
            for (int j = hw; j < width_ + hw; ++j, ++ind) {
                jacobians[2][ind] = 0.0;
                for (int n = -hw; n <= hw; ++n)
                    for (int m = -hw; m <= hw; ++m, ++ind)
                        jacobians[2][ind] -= scene(i + n, j + m);
            }

    }

    return true;

};


int imax(int a, int b) {
    if (a >= b) return a;
    return b;
}

MatrixXd Thresher::fftconvolve(MatrixXd data, MatrixXd kernel)
{
    int dh = data.rows(), dw = data.cols(),
        kh = kernel.rows(), kw = kernel.cols(),
        ky = floor(0.5 * (kh - 1)), kx = floor(0.5 * (kw - 1));

    if (kh >= dh || kw >= dw) {
        throw ConvolutionError;
    }

    // Pad out the kernel to match the data.
    MatrixXd padded = MatrixXd::Zero(dh, dw);

    padded.block(dh - kx - 1, dw - ky - 1, kx, ky) = kernel.block(0, 0, kx, ky);
    padded.block(0, 0, kh - kx, kw - ky) = kernel.block(kx, ky, kh - kx, kw - ky);
    padded.block(dh - kx - 1, 0, kx, kw - ky) = kernel.block(0, ky, kx, kw - ky);
    padded.block(0, dw - ky - 1, kh - kx, ky) = kernel.block(kx, 0, kh - kx, ky);

    /* int na = a.rows() * floor(0.5 * a.cols() + 1), */
    /*     nb = b.rows() * floor(0.5 * b.cols() + 1), */
    /*     nx = imax(a.rows(), b.rows()), */
    /*     ny = imax(a.cols(), b.cols()), */
    /*     nmx = nx * floor(0.5 * ny + 1); */

    /* MatrixXd result( */

    /* fftw_complex *fa = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * na), */
    /*              *fb = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nb), */
    /*              *mid = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nmx); */

    /* fftw_plan pa = fftw_plan_dft_r2c_2d(a.rows(), a.cols(), */
    /*                                     &(a(0)), fa, FFTW_ESTIMATE), */
    /*           pb = fftw_plan_dft_r2c_2d(b.rows(), b.cols(), */
    /*                                     &(b(0)), fb, FFTW_ESTIMATE), */
    /*           pb = fftw_plan_dft_r2c_2d(b.rows(), b.cols(), */
    /*                                     &(b(0)), fb, FFTW_ESTIMATE); */

    /* fftw_execute(pa); */
    /* fftw_execute(pb); */



    /* fftw_destroy_plan(pa); */
    /* fftw_destroy_plan(pb); */
    /* fftw_destroy_plan(pinv); */

    /* fftw_free(fa); */
    /* fftw_free(fb); */
    /* fftw_free(mid); */
}
