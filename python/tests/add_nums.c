#include <stdio.h>
#include <stdlib.h>

#include "add_nums.h"

float add_float(float a, float b)
{
    return a + b;
}

int add_int(int a, int b)
{
    return a + b;
}

int add_float_ref(float *a, float *b, float *c)
{
    *c = *a + *b;
    return 0;
}
