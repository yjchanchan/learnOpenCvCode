#include "features2d.hpp"
#include <algorithm>
#include "../core/define.hpp"

/** @brief Applies a separable linear filter to an image.

The function applies a separable linear filter to the image. That is, first, every row of src is
filtered with the 1D kernel kernelX. Then, every column of the result is filtered with the 1D
kernel kernelY. The final result shifted by delta is stored in dst .
*/
void sepFilter2D(uint8_t *src, int16_t *dst, int w, int h, double *kernelX, double *kernelY, int ksize)
{
}

static void getSobelKernels(double *kx, double *ky, int dx, int dy, int ker_size, bool normalize)
{
    int i, j;
    int ksizeX = ker_size; // default 3
    int ksizeY = ker_size;
    if (ksizeX == 1 && dx > 0)
        ksizeX = 3;
    if (ksizeY == 1 && dy > 0)
        ksizeY = 3;

    int tmp_size = std::max(ksizeX, ksizeY) + 1;
    int *kerI = (int *)calloc(tmp_size, sizeof(int));

    for (int k = 0; k < 2; k++)
    {
        // first x, then y
        double *kernel = (k == 0) ? kx : ky;
        int order = (k == 0) ? dx : dy;
        int ksize = (k == 0) ? ksizeX : ksizeY;

        assert(ksize > order);

        if (ksize == 1)
            kerI[0] = 1;
        else if (ksize == 3)
        {
            if (order == 0)
                kerI[0] = 1, kerI[1] = 2, kerI[2] = 1;
            else if (order == 1)
                kerI[0] = -1, kerI[1] = 0, kerI[2] = 1;
            else
                kerI[0] = 1, kerI[1] = -2, kerI[2] = 1;
        }
        else
        {
            int oldval, newval;
            kerI[0] = 1;
            for (i = 0; i < ksize; i++)
                kerI[i + 1] = 0;

            for (i = 0; i < ksize - order - 1; i++)
            {
                oldval = kerI[0];
                for (j = 1; j <= ksize; j++)
                {
                    newval = kerI[j] + kerI[j - 1];
                    kerI[j - 1] = oldval;
                    oldval = newval;
                }
            }

            for (i = 0; i < order; i++)
            {
                oldval = -kerI[0];
                for (j = 1; j <= ksize; j++)
                {
                    newval = kerI[j - 1] - kerI[j];
                    kerI[j - 1] = oldval;
                    oldval = newval;
                }
            }
        }

        double scale = !normalize ? 1. : 1. / (1 << (ksize - order - 1));
        for (int i = 0; i < ksize; i++)
        {
            kernel[i] = kerI[i] * scale;
        }
    }

    safeFree(kerI);
}

void sobel(uint8_t *src, int16_t *dst, int w, int h, int dx, int dy, int ksize, double scale)
{
    double *kx = (double *)calloc(ksize, sizeof(double));
    double *ky = (double *)calloc(ksize, sizeof(double));
    getSobelKernels(kx, ky, dx, dy, ksize, false);
    if (scale != 1)
    {
        // usually the smoothing part is the slowest to compute,
        // so try to scale it instead of the faster differentiating part
        if (dx == 0)
        {
            for (int i = 0; i < ksize; i++)
            {
                kx[i] *= scale;
            }
        }
        else
        {
            for (int i = 0; i < ksize; i++)
            {
                ky[i] *= scale;
            }
        }
    }

    sepFilter2D(src, dst, w, h, kx, ky, ksize);

    safeFree(kx);
    safeFree(ky);
}

void cornerHarris(uint8_t *src, uint8_t *dst, int w, int h, int block_size, int grad_size, double k)
{
    // calc scale coeff
    double scale = (double)(1 << ((grad_size > 0 ? grad_size : 3) - 1)) * block_size;
    if (grad_size < 0)
        scale *= 2.0;

    scale *= 255.0;
    scale = 1.0 / scale;

    // calc gradient
    int16_t *Ix = (int16_t *)calloc(w * h, sizeof(int16_t));
    int16_t *Iy = (int16_t *)calloc(w * h, sizeof(int16_t));

    sobel(src, Ix, w, h, 1, 0, grad_size, scale);
    sobel(src, Iy, w, h, 0, 1, grad_size, scale);

    safeFree(Ix);
    safeFree(Iy);
}