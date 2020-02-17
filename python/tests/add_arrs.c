#include <stdio.h>
#include <stdlib.h>

#include "add_arrs.h"

int add_int_array(int *a, int *b, int *c, int n)
{
    int i;
    for(i = 0; i < n; ++i)
        c[i] = a[i] + b[i];

    return 0;
}

float dot_product(float *a, float *b, int n)
{
    float res = 0.0;
    int i;
    for(i = 0; i < n; ++i)
        res += a[i] * b[i];

    return res;
}
