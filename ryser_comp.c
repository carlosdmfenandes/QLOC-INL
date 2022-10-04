#include <math.h>
#include <complex.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>

#define COMP complex

unsigned long long  bits2grayll(unsigned long long bits){
    unsigned int gray = bits ^ (bits >> 1);
    return gray;
};

unsigned int bits2gray(unsigned int bits){
    unsigned int gray = bits ^ (bits >> 1);
    return gray;
};

/*
int graybitshift(unsigned int *pbits){
    int i = 0;
    while(*pbits%2){
        i++:
        *pbits >> 1;
    };
    return res;
};


int graybitshiftll2(unsigned long long bits){
    unsigned long long gray1 = bits2grayll(bits);
    unsigned long long gray2 = bits2grayll(bits - 1);
    unsigned long long d = gray1 ^ gray2;
    while(bits%2){
        i++;
        bits >> 1;
    };
    return res;
};
*/
int graybitshiftll(unsigned long long * restrict bits){
    int res = 0;
    while(*bits%2){
        *bits >>= 1;
        res++;
    };
    *bits >>= 1;
    return res;
};

double COMP prodd(const double COMP * const vec, size_t len){
    double COMP prod = 1.0;
    for(const double COMP *ptr = vec; ptr - vec < len; prod *= *ptr++);
    return prod;
};

double COMP *addvecd(double COMP *dest, const double COMP *term, size_t len){
    for(double COMP *ptr = dest; ptr - dest < len; *ptr++ += *term++);
    return dest;
};

double COMP *subvecd(double COMP *dest, const double COMP *term, size_t len){
    for(double COMP *ptr = dest; ptr - dest < len; *ptr++ -= *term++);
    return dest;
};

double COMP ryserperm(const double COMP *mat, size_t dim){
    double COMP retval = 0.0;
    double COMP * vec = calloc(dim, sizeof(double COMP));
    if(vec == NULL){
        puts(strerror(ENOMEM));
        exit(ENOMEM);
    };
    for(unsigned long long i = 0; i < (1llu << dim) - 1; i++){
        unsigned long long bits = i;
        unsigned int pos = graybitshiftll(&bits);
        if(bits%2){
            addvecd(vec, mat+pos*dim, dim);
        } else {
            subvecd(vec, mat+pos*dim, dim);
        };
//        printf("%d\n", pos);
//        for(double *p = vec; p-vec < dim; printf("%f\n", *p++));
        double COMP prod = prodd(vec, dim);
        if(i%2){
            retval -= prod;
        } else {
            retval += prod;
        };
    };
    retval = dim%2 ? -retval : retval;
    return retval;
};
