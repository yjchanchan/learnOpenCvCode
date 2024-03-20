#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

void cornerHarris(uint8_t *src, uint8_t *dst, int w, int h, int block_size = 3, int grad_size = 3, double k = 0.04);