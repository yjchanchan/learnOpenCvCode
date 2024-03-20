#include "calib.hpp"
#include "opencv2/opencv.hpp"
#include <string>
#include "../core/define.hpp"
#include "../core/utils.hpp"

int main()
{
    std::string path_in = "E:/data/404.jpg";
    // cv::Mat img_in = cv::imread(path_in, cv::IMREAD_GRAYSCALE);
    int w = 800;
    int h = 600;

    double K[9] = {573.8534, -0.6974, 406.0101,
                   0, 575.0448, 309.0112,
                   0, 0, 1};

    double dist_coeff[5] = {3, 2, 0.1, 0.1, 1};
    float *mapx = (float *)calloc(w * h, sizeof(float));
    float *mapy = (float *)calloc(w * h, sizeof(float));

    // initUndistortRectifyMap(K, K, mapx, mapy, w, h, dist_coeff, LDC_NORM);

    // ============
    LDC_LUT lut;
    float distort[11] = {
        0,
        0.001,
        0.002,
        0.003,
        0.004,
        0.005,
        0.006,
        0.008,
        0.010,
        0.012,
        0.015};
    lut.len = 11;
    lut.distort = distort;
    lut.w = w;
    lut.h = h;
    lut.cx = w >> 1;
    lut.cy = h >> 1;

    distortMappingLut(mapx, mapy, lut, LDC_MAP_BACKWARD);

    cv::Mat mapx_mat(h, w, CV_32FC1, mapx);
    cv::Mat mapy_mat(h, w, CV_32FC1, mapy);
    cv::Mat flow_map = drawMotionMap(mapx_mat, mapy_mat);

    mapx_mat.convertTo(mapx_mat, CV_8UC1);
    mapy_mat.convertTo(mapy_mat, CV_8UC1);
    cv::imwrite("./mapx.bmp", mapx_mat);
    cv::imwrite("./mapy.bmp", mapy_mat);
    cv::imwrite("./flow.bmp", flow_map);

    printf("end\n");

    safeFree(mapx);
    safeFree(mapy);
    return 0;
}