#include "calib.hpp"
#include <math.h>
#include <opencv.hpp>

// todo, undistortPoints() learn

inline float linearInterp(float v0, float v1, float w)
{
    return v0 + (v1 - v0) * w;
}

void distortMappingLut(float *mapx, float *mapy, const LDC_LUT &lut, int direction)
{
    // r_undistort * (1 + distort) = r_distort
    int width = lut.w;
    int height = lut.h;
    int w_half = width >> 1;
    int h_half = height >> 1;
    float max_r = std::sqrt(w_half * w_half + h_half * h_half); // relatived to lut[-1] param
    float step = max_r / (lut.len - 1);                         // lut[0] is always 0, for r=0

    for (int i = 0; i < height; i++)
    {
        float y = i - lut.cy;
        for (int j = 0; j < width; j++)
        {
            float x = j - lut.cx;
            float r = std::sqrt(x * x + y * y);

            bool keep_region = std::abs(r) < 1e-5;
            if (keep_region)
            {
                // mapx[i * width + j] = j;
                // mapy[i * width + j] = i;
                mapx[i * width + j] = 0;
                mapy[i * width + j] = 0;
                continue;
            }

            float map_u = 0;
            float map_v = 0;
            if (LDC_MAP_BACKWARD == direction)
            {
                // (i, j) is in undistort img
                float ratio = r / step;
                int ratio_i = std::floor(ratio);
                ratio_i = std::min(ratio_i, lut.len - 2); // r may > max_r
                float wei = ratio - ratio_i;              // can larger than 1 for extra interp

                float d0 = lut.distort[ratio_i];
                float d1 = lut.distort[ratio_i + 1];

                float dist = linearInterp(d0, d1, wei);
                float scale = 1 + dist;
                float xd = x * scale;
                float yd = y * scale;
                map_u = xd + lut.cx;
                map_v = yd + lut.cy;
            }
            else if (LDC_MAP_FORWARD == direction)
            {
            }

            mapx[i * width + j] = map_u - j;
            mapy[i * width + j] = map_v - i;
        }
    }
}

void distortMappingFishEyes(float x_in, float y_in, float *u, float *v, double *K, double *Kinv, double *dist_coeffs)
{
    const double *k = dist_coeffs; // dist_coeffs[4]

    double cx = K[0 * 3 + 2], cy = K[1 * 3 + 2];
    double fx = K[0 * 3 + 0], fy = K[1 * 3 + 1];

    // [_x, _y, _w] = Kinv \cdot [x, y, 1]
    double _x = x_in * Kinv[0] + y_in * Kinv[1] + Kinv[2];
    double _y = x_in * Kinv[3] + y_in * Kinv[4] + Kinv[5];
    double _w = x_in * Kinv[6] + y_in * Kinv[7] + Kinv[8];
    // (x, y) is 3d point without dist, at z=1 plane
    double w = 1. / _w;
    double x = _x * w;
    double y = _y * w;

    double r = std::sqrt(x * x + y * y);
    double theta = std::atan(r);
    double theta2 = theta * theta, theta4 = theta2 * theta2, theta6 = theta4 * theta2, theta8 = theta4 * theta4;
    double theta_d = theta * (1 + k[0] * theta2 + k[1] * theta4 + k[2] * theta6 + k[3] * theta8);

    double scale = (r == 0) ? 1.0 : theta_d / r;
    double xd = x * scale;
    double yd = y * scale;
    *u = fx * xd + cx;
    *v = fy * yd + cy;
}

void distortMapping(float x_in, float y_in, float *u, float *v, double *K, double *Kinv, double *dist_coeffs)
{
    double k1 = dist_coeffs[0]; // dist_coeffs[5]
    double k2 = dist_coeffs[1];
    double p1 = dist_coeffs[2];
    double p2 = dist_coeffs[3];
    double k3 = dist_coeffs[4];

    double k4 = 0;
    double k5 = 0;
    double k6 = 0;
    double s1 = 0;
    double s2 = 0;
    double s3 = 0;
    double s4 = 0;
    double tauX = 0;
    double tauY = 0;

    double cx = K[0 * 3 + 2], cy = K[1 * 3 + 2];
    double fx = K[0 * 3 + 0], fy = K[1 * 3 + 1];

    // [_x, _y, _w] = Kinv \cdot [x, y, 1]
    double _x = x_in * Kinv[0] + y_in * Kinv[1] + Kinv[2];
    double _y = x_in * Kinv[3] + y_in * Kinv[4] + Kinv[5];
    double _w = x_in * Kinv[6] + y_in * Kinv[7] + Kinv[8];
    // (x, y) is 3d point without dist, at z=1 plane
    double w = 1. / _w;
    double x = _x * w;
    double y = _y * w;
    double x2 = x * x;
    double y2 = y * y;
    double r2 = x2 + y2;
    double _2xy = 2 * x * y;
    double kr = (1 + ((k3 * r2 + k2) * r2 + k1) * r2) / (1 + ((k6 * r2 + k5) * r2 + k4) * r2);
    // (xd, yd) is 3d point after dist, at z=1 plane
    double xd = (x * kr + p1 * _2xy + p2 * (r2 + 2 * x2) + s1 * r2 + s2 * r2 * r2);
    double yd = (y * kr + p1 * (r2 + 2 * y2) + p2 * _2xy + s3 * r2 + s4 * r2 * r2);
    // no considering for trapezoidal distortion of tilted image sensor
    // (u, v) is img point of (xd, yd) with K inner matrix
    *u = fx * xd + cx;
    *v = fy * yd + cy;
}

void initUndistortRectifyMap(
    double K[9],    // cameraMatrix
    double Knew[9], // newCameraMatrix
    float *map1,    // output dist map x
    float *map2,    // output dist map y
    int width, int height, double *dist_coeffs, int mode)
{
    cv::Mat mat_k(3, 3, CV_64FC1, Knew);
    cv::Mat mat_inv;
    cv::invert(mat_k, mat_inv, cv::DECOMP_LU); // todo, use inv method by self

    double *Kinv = (double *)mat_inv.data;

    // (j, i) is 2d img point in new undist dst img
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            float u = 0.0f, v = 0.0f;
            if (mode == LDC_NORM)
            {
                distortMapping(j, i, &u, &v, K, Kinv, dist_coeffs);
            }
            else if (mode == LDC_FISHEYE)
            {
                distortMappingFishEyes(j, i, &u, &v, K, Kinv, dist_coeffs);
            }
            map1[i * width + j] = u;
            map2[i * width + j] = v;
        }
    }
}
