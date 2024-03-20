#pragma once
#include <stdlib.h>

#define safeFree(ptr)  \
    if (ptr)           \
    {                  \
        free(ptr);     \
        ptr = nullptr; \
    }
