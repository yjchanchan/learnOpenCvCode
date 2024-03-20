#include "features2d.hpp"
#include "opencv2/opencv.hpp"
#include <string>
#include "../core/define.hpp"

int main()
{
    std::string path_in = "E:/data/404.jpg";
    cv::Mat img_in = cv::imread(path_in, cv::IMREAD_GRAYSCALE);
    // cv::imshow("hh", img_in);
    // cv::waitKey(0);
    int w = img_in.cols;
    int h = img_in.rows;
    uint8_t *harris_res = (uint8_t *)calloc(w * h, sizeof(uint8_t));
    cornerHarris(img_in.data, harris_res, w, h);

    printf("end\n");

    safeFree(harris_res);
    return 0;
}