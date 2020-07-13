/**
 * Copyright by Ruman Gerst
 * Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
 * https://www.leibniz-hki.de/en/applied-systems-biology.html
 * HKI-Center for Systems Biology of Infection
 * Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
 * Adolf-Reichwein-Straße 23, 07745 Jena, Germany
 *
 * This code is licensed under BSD 2-Clause
 * See the LICENSE file provided with this code for the full license.
 */

#include "misaxx_convolve.h"
#include <cmath>


namespace cv::images {
    using grayscale8u = cv::Mat1b;
    using grayscale32f = cv::Mat1f;
    using mask = cv::Mat1b;
    using labels = cv::Mat1i;
    using complex = cv::Mat2f;
}

namespace {

    cv::Size get_fft_size(const cv::Mat &img, const cv::Mat &kernel) {
        return cv::Size(img.size().width + kernel.size().width - 1,
                        img.size().height + kernel.size().height - 1);
    }

    cv::images::grayscale32f
    fftunpad(const cv::images::grayscale32f &deconvolved, const cv::Size &target_size, const cv::Size &source_size) {
        cv::Size ap{};
        if (source_size.width % 2 == 0)
            ap.width = 1;
        if (source_size.height % 2 == 0)
            ap.height = 1;

        int padded_width = deconvolved.size().width;
        int padded_height = deconvolved.size().height;

        int bleft = (padded_width - source_size.width - ap.width) / 2;
        int btop = (padded_height - source_size.height - ap.height) / 2;
        cv::Rect roi{bleft, btop, source_size.width, source_size.height};
        cv::images::grayscale32f result{source_size, 0};
        deconvolved(roi).copyTo(result);
        return result;
    }



    cv::images::grayscale32f fftpad(const cv::images::grayscale32f &img, const cv::Size &target_size, bool shift = false) {
        cv::Size ap{};
        if (img.size().width % 2 == 0)
            ap.width = 1;
        if (img.size().height % 2 == 0)
            ap.height = 1;

        cv::Size c{};
        c.width = (target_size.width - img.size().width - ap.width) / 2;
        c.height = (target_size.height - img.size().height - ap.height) / 2;

        int bleft = c.width;
        int btop = c.height;
        int bright = c.width + ap.width;
        int bbottom = c.height + ap.height;

        // Further pad to optimal FFT size
        {
            int currentWidth = bleft + bright + img.size().width;
            int currentHeight = btop + bbottom + img.size().height;
            int optimalWidth = cv::getOptimalDFTSize(currentWidth);
            int optimalHeight = cv::getOptimalDFTSize(currentHeight);

            int dow = optimalWidth - currentWidth;
            int doh = optimalHeight - currentHeight;
            int ow0 = dow / 2 + 1;
            int ow1 = dow - ow0;
            int oh0 = doh / 2 + 1;
            int oh1 = doh - oh0;

            // Add to padding
            bleft += ow0;
            btop += oh0;
            bright += ow1;
            bbottom += oh1;
        }

        cv::images::grayscale32f padded = img;
        cv::copyMakeBorder(img, padded, btop, bbottom, bleft, bright, cv::BORDER_CONSTANT, cv::Scalar::all(0));


        if (shift) {
            cv::images::grayscale32f tmp {};
            cv::repeat(padded, 2, 2, tmp);
            int sx = padded.cols;
            int sy = padded.rows;
            padded = tmp(cv::Rect(sx / 2, sy / 2, sx, sy));
        }

        return padded;
    }

    cv::images::complex fft(const cv::images::grayscale32f &padded) {
        cv::images::complex result{padded.size(), cv::Vec2f{0, 0}};
        cv::dft(padded, result, cv::DFT_COMPLEX_OUTPUT);
        return result;
    }

    cv::images::grayscale32f ifft(const cv::images::complex &fft) {
        cv::images::grayscale32f result{fft.size(), 0};
        cv::idft(fft, result, cv::DFT_REAL_OUTPUT | cv::DFT_SCALE);
        return result;
    }

    inline cv::Vec2f complex_add(cv::Vec2f a, cv::Vec2f b) {
        return cv::Vec2f{a[0] + b[0], a[1] + b[1]};
    }

    inline cv::Vec2f complex_mul(cv::Vec2f a, cv::Vec2f b) {
        return cv::Vec2f{a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]};
    }

    inline cv::Vec2f scalar_mul(cv::Vec2f a, float b) {
        return cv::Vec2f{a[0] * b, a[1] * b};
    }

    inline cv::Vec2f complex_div(cv::Vec2f a, cv::Vec2f b) {
        const float A0 = a[0] * b[0] + a[1] * b[1];
        const float B0 = a[1] * b[0] - a[0] * b[1];
        const float D = b[0] * b[0] + b[1] * b[1];
        return cv::Vec2f{A0 / D, B0 / D};
    }

    inline float complex_abs(cv::Vec2f a) {
        return std::sqrt(a[0] * a[0] + a[1] * a[1]);
    }

}

/**
 * Perform FFT on input image and convole/multiply image in Fourier space with following IFFT.
 *
 * @param input Image on which the FFT filter is applied.
 * @param mask Which is used as filter.
 * @param dst Resulting image after IFFT.
 */
void misaxx_convolve::work(const cv::Mat &input, const cv::Mat1f &mask, cv::Mat &dst) {

    cv::Size target_size = get_fft_size(input, mask);

    cv::images::complex Y = fft(fftpad(input, target_size));
    cv::images::grayscale32f H = fftpad(mask, target_size, true);
    cv::images::complex X{Y.size(), cv::Vec2f(0, 0)};

    for(int y = 0; y < Y.rows; ++y) {

        const float *h = H[y];
        const cv::Vec2f *rY = Y[y];
        cv::Vec2f *rX = X[y];

        for(int x = 0; x < Y.cols; ++x) {
            rX[x] = scalar_mul(rY[x], h[x]);
        }
    }

    cv::images::grayscale32f multiplied = fftunpad(ifft(X), target_size, input.size());

    dst = multiplied;

}
