#include <math.h>

#include "fourier.h"

void normalize(complex s[], int n) {
    for (int k = 0; k < n; k++) {
        s[k].a /= n;
        s[k].b /= n;
    }
}

void nft(complex s[], complex t[], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k].a = 0;
        t[k].b = 0;

        for (int j = 0; j < n; j++) {
            double x = sign * 2 * PI * k * j / n;

            double cosx = cos(x);
            double sinx = sin(x);

            t[k].a += s[j].a * cosx - s[j].b * sinx;
            t[k].b += s[j].a * sinx + s[j].b * cosx;
        }
    }
}

void nft_forward(complex s[], complex t[], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(complex t[], complex s[], int n) {
    nft(t, s, n, 1);
    normalize(s, n);
}

void fft(complex s[], complex t[], int n, int sign) {
    if (n == 1) {
        t[0] = s[0];
        return;
    }

    int m = n / 2;
    complex even[m], odd[m];
    complex even_t[m], odd_t[m];

    for (int i = 0; i < m; i++) {
        even[i] = s[2 * i];
        odd[i]  = s[2 * i + 1];
    }

    fft(even, even_t, m, sign);
    fft(odd, odd_t, m, sign);

    for (int k = 0; k < m; k++) {
        double angle = sign * -2.0 * PI * k / n;
        double cos_a = cos(angle);
        double sin_a = sin(angle);

        complex w;
        w.a = cos_a;
        w.b = sin_a;

        complex wodd;
        wodd.a = w.a * odd_t[k].a - w.b * odd_t[k].b;
        wodd.b = w.a * odd_t[k].b + w.b * odd_t[k].a;

        t[k].a     = even_t[k].a + wodd.a;
        t[k].b     = even_t[k].b + wodd.b;

        t[k+m].a   = even_t[k].a - wodd.a;
        t[k+m].b   = even_t[k].b - wodd.b;
    }
}

void fft_forward(complex s[], complex t[], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(complex t[], complex s[], int n) {
    fft(t, s, n, 1);
    normalize(s, n);
}

void fft_forward_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
}

void fft_inverse_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
}

void filter(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;
            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x].a = g * input[y][x].a;
            output[y][x].b = g * input[y][x].b;
        }
    }
}

void filter_lp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}
