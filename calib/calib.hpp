
enum
{
    LDC_NORM = 0,
    LDC_FISHEYE
};

enum
{
    LDC_MAP_BACKWARD = 0,
    LDC_MAP_FORWARD
};

struct LDC_LUT
{
    float *distort; // r=0, to r=max
    int len;
    int w;
    int h;
    float cx; // center point pos without distort
    float cy;
};

void distortMappingLut(float *mapx, float *mapy, const LDC_LUT &lut, int direction);

void initUndistortRectifyMap(
    double K[9],    // cameraMatrix
    double Knew[9], // newCameraMatrix
    float *map1,    // output dist map x
    float *map2,    // output dist map y
    int width, int height, double *distCoeffs, int mode);