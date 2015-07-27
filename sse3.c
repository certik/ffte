/*
      FFTE: A FAST FOURIER TRANSFORM PACKAGE

      (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED
                 BY
          DAISUKE TAKAHASHI
          GRADUATE SCHOOL OF SYSTEMS AND INFORMATION ENGINEERING
          UNIVERSITY OF TSUKUBA
          1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
          E-MAIL: daisuke@cs.tsukuba.ac.jp


      RADIX-2, 3, 4, 5 AND 8 FFT KERNEL ROUTINE (FOR SSE3)

      C SOURCE PROGRAM

      WRITTEN BY DAISUKE TAKAHASHI
*/

#include <pmmintrin.h>

__m128d ZMUL(__m128d a, __m128d b);

static __inline __m128d ZMUL(__m128d a, __m128d b)
{
    __m128d ar, ai;

    ar = _mm_movedup_pd(a);       /* ar = [a.r a.r] */
    ar = _mm_mul_pd(ar, b);       /* ar = [a.r*b.r a.r*b.i] */
    ai = _mm_unpackhi_pd(a, a);   /* ai = [a.i a.i] */
    b = _mm_shuffle_pd(b, b, 1);  /* b = [b.i b.r] */
    ai = _mm_mul_pd(ai, b);       /* ai = [a.i*b.i a.i*b.r] */

    return _mm_addsub_pd(ar, ai); /* [a.r*b.r-a.i*b.i a.r*b.i+a.i*b.r] */
}

int fft2_(double *a, double *b, int *m)
{
    int i, i0, i1;
    /* double x0, y0, x1, y1; */
    __m128d t0, t1;

    for (i = 0; i < *m ; i++) {
	i0 = i << 1;
	i1 = i0 + (*m << 1);
	/* x0 = a[i0];
	y0 = a[i0 + 1];
	x1 = a[i1];
	y1 = a[i1 + 1]; */
	t0 = _mm_load_pd(&a[i0]);
	t1 = _mm_load_pd(&a[i1]);
	/* b[i0] = x0 + x1;
	b[i0 + 1] = y0 + y1;
	b[i1] = x0 - x1;
	b[i1 + 1] = y0 - y1; */
	_mm_store_pd(&b[i0], _mm_add_pd(t0, t1));
	_mm_store_pd(&b[i1], _mm_sub_pd(t0, t1));
    }
    return 0;
}

int fft3a_(double *a, double *b, double *w, int *l)
{
    /* static double c31 = .86602540378443865;
    static double c32 = .5; */
    static __m128d c31, c32;

    int j, j0, j1, j2, j3, j4, j5;
    /* double x0, y0, x1, y1, x2, y2, wi1, wi2, wr1, wr2; */
    __m128d t0, t1, t2, t3, w1, w2;

    c31 = _mm_set1_pd(0.86602540378443865);
    c32 = _mm_set1_pd(0.5);

    for (j = 0; j < *l; j++) {
        j0 = j << 1;
	j1 = j0 + (*l << 1);
	j2 = j1 + (*l << 1);
	j3 = j * 6;
	j4 = j3 + 2;
	j5 = j4 + 2;
	/* wr1 = w[j0];
	wi1 = w[j0 + 1];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1; */
	w1 = _mm_load_pd(&w[j0]);
	w2 = ZMUL(w1, w1);
	/* x0 = a[j1] + a[j2];
	y0 = a[j1 + 1] + a[j2 + 1];
	x1 = a[j0] - c32 * x0;
	y1 = a[j0 + 1] - c32 * y0;
	x2 = c31 * (a[j1 + 1] - a[j2 + 1]);
	y2 = c31 * (a[j2] - a[j1]); */
	t1 = _mm_load_pd(&a[j1]);
	t2 = _mm_load_pd(&a[j2]);
	t0 = _mm_add_pd(t1, t2);
	t2 = _mm_xor_pd(_mm_sub_pd(t1, t2), _mm_set_sd(-0.0));
	t2 = _mm_mul_pd(c31, _mm_shuffle_pd(t2, t2, 1));
	t3 = _mm_load_pd(&a[j0]);
	t1 = _mm_sub_pd(t3, _mm_mul_pd(c32, t0));
	/* b[j3] = a[j0] + x0;
	b[j3 + 1] = a[j0 + 1] + y0;
	b[j4] = wr1 * (x1 + x2) - wi1 * (y1 + y2);
	b[j4 + 1] = wr1 * (y1 + y2) + wi1 * (x1 + x2);
	b[j5] = wr2 * (x1 - x2) - wi2 * (y1 - y2);
	b[j5 + 1] = wr2 * (y1 - y2) + wi2 * (x1 - x2); */
	_mm_store_pd(&b[j3], _mm_add_pd(t3, t0));
	_mm_store_pd(&b[j4], ZMUL(w1, _mm_add_pd(t1, t2)));
	_mm_store_pd(&b[j5], ZMUL(w2, _mm_sub_pd(t1, t2)));
    }
    return 0;
}

int fft3b_(double *a, double *b, double *w, int *m, int *l)
{
    /* static double c31 = .86602540378443865;
    static double c32 = .5; */
    static __m128d c31, c32;

    int i, i0, i1, i2, i3, i4, i5, j, j0;
    /* double x0, y0, x1, y1, x2, y2, wi1, wi2, wr1, wr2; */
    __m128d t0, t1, t2, t3, w1, w2;

    c31 = _mm_set1_pd(0.86602540378443865);
    c32 = _mm_set1_pd(0.5);

    for (i = 0; i < *m; i++) {
        i0 = i << 1;
	i1 = i0 + (*m * *l << 1);
	i2 = i1 + (*m * *l << 1);
	i3 = i << 1;
	i4 = i3 + (*m << 1);
	i5 = i4 + (*m << 1);
	/* x0 = a[i1] + a[i2];
	y0 = a[i1 + 1] + a[i2 + 1];
	x1 = a[i0] - c32 * x0;
	y1 = a[i0 + 1] - c32 * y0;
	x2 = c31 * (a[i1 + 1] - a[i2 + 1]);
	y2 = c31 * (a[i2] - a[i1]); */
	t1 = _mm_load_pd(&a[i1]);
	t2 = _mm_load_pd(&a[i2]);
	t0 = _mm_add_pd(t1, t2);
	t2 = _mm_xor_pd(_mm_sub_pd(t1, t2), _mm_set_sd(-0.0));
	t2 = _mm_mul_pd(c31, _mm_shuffle_pd(t2, t2, 1));
	t3 = _mm_load_pd(&a[i0]);
	t1 = _mm_sub_pd(t3, _mm_mul_pd(c32, t0));
	/* b[i3] = a[i0] + x0;
	b[i3 + 1] = a[i0 + 1] + y0;
	b[i4] = x1 + x2;
	b[i4 + 1] = y1 + y2;
	b[i5] = x1 - x2;
	b[i5 + 1] = y1 - y2; */
	_mm_store_pd(&b[i3], _mm_add_pd(t3, t0));
	_mm_store_pd(&b[i4], _mm_add_pd(t1, t2));
	_mm_store_pd(&b[i5], _mm_sub_pd(t1, t2));
    }
    for (j = 1; j < *l; j++) {
        j0 = j << 1;
	/* wr1 = w[j0];
	wi1 = w[j0 + 1];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1; */
	w1 = _mm_load_pd(&w[j0]);
	w2 = ZMUL(w1, w1);
	for (i = 0; i < *m; i++) {
	    i0 = (i << 1) + (j * *m << 1);
	    i1 = i0 + (*m * *l << 1);
	    i2 = i1 + (*m * *l << 1);
	    i3 = (i << 1) + (j * *m * 6);
	    i4 = i3 + (*m << 1);
	    i5 = i4 + (*m << 1);
	    /* x0 = a[i1] + a[i2];
	    y0 = a[i1 + 1] + a[i2 + 1];
	    x1 = a[i0] - x0 * .5;
	    y1 = a[i0 + 1] - y0 * .5;
	    x2 = c31 * (a[i1 + 1] - a[i2 + 1]);
	    y2 = c31 * (a[i2] - a[i1]); */
	    t1 = _mm_load_pd(&a[i1]);
	    t2 = _mm_load_pd(&a[i2]);
	    t0 = _mm_add_pd(t1, t2);
	    t2 = _mm_xor_pd(_mm_sub_pd(t1, t2), _mm_set_sd(-0.0));
	    t2 = _mm_mul_pd(c31, _mm_shuffle_pd(t2, t2, 1));
	    t3 = _mm_load_pd(&a[i0]);
	    t1 = _mm_sub_pd(t3, _mm_mul_pd(c32, t0));
	    /* b[i3] = a[i0] + x0;
	    b[i3 + 1] = a[i0 + 1] + y0;
	    b[i4] = wr1 * (x1 + x2) - wi1 * (y1 + y2);
	    b[i4 + 1] = wr1 * (y1 + y2) + wi1 * (x1 + x2);
	    b[i5] = wr2 * (x1 - x2) - wi2 * (y1 - y2);
	    b[i5 + 1] = wr2 * (y1 - y2) + wi2 * (x1 - x2); */
	    _mm_store_pd(&b[i3], _mm_add_pd(t3, t0));
	    _mm_store_pd(&b[i4], ZMUL(w1, _mm_add_pd(t1, t2)));
	    _mm_store_pd(&b[i5], ZMUL(w2, _mm_sub_pd(t1, t2)));
	}
    }
    return 0;
}

int fft4a_(double *a, double *b, double *w, int *l)
{
    int j, j0, j1, j2, j3, j4, j5, j6, j7;
    /* double x0, y0, x1, y1, x2, y2, x3, y3, wi1, wi2, wi3, wr1, wr2, wr3; */
    __m128d t0, t1, t2, t3, t4, w1, w2, w3;

    for (j = 0; j < *l; j++) {
	j0 = j << 1;
	j1 = j0 + (*l << 1);
	j2 = j1 + (*l << 1);
	j3 = j2 + (*l << 1);
	j4 = j << 3;
	j5 = j4 + 2;
	j6 = j5 + 2;
	j7 = j6 + 2;
	/* wr1 = w[j0];
	wi1 = w[j0 + 1];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	wr3 = wr1 * wr2 - wi1 * wi2;
	wi3 = wr1 * wi2 + wi1 * wr2; */
	w1 = _mm_load_pd(&w[j0]);
	w2 = ZMUL(w1, w1);
	w3 = ZMUL(w1, w2);
	/* x0 = a[j0] + a[j2];
	y0 = a[j0 + 1] + a[j2 + 1];
	x1 = a[j0] - a[j2];
	y1 = a[j0 + 1] - a[j2 + 1];
	x2 = a[j1] + a[j3];
	y2 = a[j1 + 1] + a[j3 + 1];
	x3 = a[j1 + 1] - a[j3 + 1];
	y3 = a[j3] - a[j1]; */
	t0 = _mm_load_pd(&a[j0]);
	t2 = _mm_load_pd(&a[j2]);
	t1 = _mm_sub_pd(t0, t2);
	t0 = _mm_add_pd(t0, t2);
	t3 = _mm_load_pd(&a[j1]);
	t4 = _mm_load_pd(&a[j3]);
	t2 = _mm_add_pd(t3, t4);
	t3 = _mm_xor_pd(_mm_sub_pd(t3, t4), _mm_set_sd(-0.0));
	t3 = _mm_shuffle_pd(t3, t3, 1);
	/* b[j4] = x0 + x2;
	b[j4 + 1] = y0 + y2;
	b[j6] = wr2 * (x0 - x2) - wi2 * (y0 - y2);
	b[j6 + 1] = wr2 * (y0 - y2) + wi2 * (x0 - x2);
	b[j5] = wr1 * (x1 + x3) - wi1 * (y1 + y3);
	b[j5 + 1] = wr1 * (y1 + y3) + wi1 * (x1 + x3);
	b[j7] = wr3 * (x1 - x3) - wi3 * (y1 - y3);
	b[j7 + 1] = wr3 * (y1 - y3) + wi3 * (x1 - x3); */
	_mm_store_pd(&b[j4], _mm_add_pd(t0, t2));
	_mm_store_pd(&b[j6], ZMUL(w2, _mm_sub_pd(t0, t2)));
	_mm_store_pd(&b[j5], ZMUL(w1, _mm_add_pd(t1, t3)));
	_mm_store_pd(&b[j7], ZMUL(w3, _mm_sub_pd(t1, t3)));
    }
    return 0;
}

int fft4b_(double *a, double *b, double *w, int *m, int *l)
{
    int i, i0, i1, i2, i3, i4, i5, i6, i7, j, j0;
    /* double x0, y0, x1, y1, x2, y2, x3, y3, wi1, wi2, wi3, wr1, wr2, wr3; */
    __m128d t0, t1, t2, t3, t4, w1, w2, w3;

    for (i = 0; i < *m; i++) {
	i0 = i << 1;
	i1 = i0 + (*m * *l << 1);
	i2 = i1 + (*m * *l << 1);
	i3 = i2 + (*m * *l << 1);
	i4 = i << 1;
	i5 = i4 + (*m << 1);
	i6 = i5 + (*m << 1);
	i7 = i6 + (*m << 1);
	/* x0 = a[i0] + a[i2];
	y0 = a[i0 + 1] + a[i2 + 1];
	x1 = a[i0] - a[i2];
	y1 = a[i0 + 1] - a[i2 + 1];
	x2 = a[i1] + a[i3];
	y2 = a[i1 + 1] + a[i3 + 1];
	x3 = a[i1 + 1] - a[i3 + 1];
	y3 = a[i3] - a[i1]; */
	t0 = _mm_load_pd(&a[i0]);
	t2 = _mm_load_pd(&a[i2]);
	t1 = _mm_sub_pd(t0, t2);
	t0 = _mm_add_pd(t0, t2);
	t3 = _mm_load_pd(&a[i1]);
	t4 = _mm_load_pd(&a[i3]);
	t2 = _mm_add_pd(t3, t4);
	t3 = _mm_xor_pd(_mm_sub_pd(t3, t4), _mm_set_sd(-0.0));
	t3 = _mm_shuffle_pd(t3, t3, 1);
	/* b[i4] = x0 + x2;
	b[i4 + 1] = y0 + y2;
	b[i6] = x0 - x2;
	b[i6 + 1] = y0 - y2;
	b[i5] = x1 + x3;
	b[i5 + 1] = y1 + y3;
	b[i7] = x1 - x3;
	b[i7 + 1] = y1 - y3; */
	_mm_store_pd(&b[i4], _mm_add_pd(t0, t2));
	_mm_store_pd(&b[i6], _mm_sub_pd(t0, t2));
	_mm_store_pd(&b[i5], _mm_add_pd(t1, t3));
	_mm_store_pd(&b[i7], _mm_sub_pd(t1, t3));
    }
    for (j = 1; j < *l; j++) {
	j0 = j << 1;
	/* wr1 = w[j0];
	wi1 = w[j0 + 1];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	wr3 = wr1 * wr2 - wi1 * wi2;
	wi3 = wr1 * wi2 + wi1 * wr2; */
	w1 = _mm_load_pd(&w[j0]);
	w2 = ZMUL(w1, w1);
	w3 = ZMUL(w1, w2);
	for (i = 0; i < *m; i++) {
	    i0 = (i << 1) + (j * *m << 1);
	    i1 = i0 + (*m * *l << 1);
	    i2 = i1 + (*m * *l << 1);
	    i3 = i2 + (*m * *l << 1);
	    i4 = (i << 1) + (j * *m << 3);
	    i5 = i4 + (*m << 1);
	    i6 = i5 + (*m << 1);
	    i7 = i6 + (*m << 1);
	    /* x0 = a[i0] + a[i2];
	    y0 = a[i0 + 1] + a[i2 + 1];
	    x1 = a[i0] - a[i2];
	    y1 = a[i0 + 1] - a[i2 + 1];
	    x2 = a[i1] + a[i3];
	    y2 = a[i1 + 1] + a[i3 + 1];
	    x3 = a[i1 + 1] - a[i3 + 1];
	    y3 = a[i3] - a[i1]; */
	    t0 = _mm_load_pd(&a[i0]);
	    t2 = _mm_load_pd(&a[i2]);
	    t1 = _mm_sub_pd(t0, t2);
	    t0 = _mm_add_pd(t0, t2);
	    t3 = _mm_load_pd(&a[i1]);
	    t4 = _mm_load_pd(&a[i3]);
	    t2 = _mm_add_pd(t3, t4);
	    t3 = _mm_xor_pd(_mm_sub_pd(t3, t4), _mm_set_sd(-0.0));
	    t3 = _mm_shuffle_pd(t3, t3, 1);
	    /* b[i4] = x0 + x2;
	    b[i4 + 1] = y0 + y2;
	    b[i6] = wr2 * (x0 - x2) - wi2 * (y0 - y2);
	    b[i6 + 1] = wr2 * (y0 - y2) + wi2 * (x0 - x2);
	    b[i5] = wr1 * (x1 + x3) - wi1 * (y1 + y3);
	    b[i5 + 1] = wr1 * (y1 + y3) + wi1 * (x1 + x3);
	    b[i7] = wr3 * (x1 - x3) - wi3 * (y1 - y3);
	    b[i7 + 1] = wr3 * (y1 - y3) + wi3 * (x1 - x3); */
	    _mm_store_pd(&b[i4], _mm_add_pd(t0, t2));
	    _mm_store_pd(&b[i6], ZMUL(w2, _mm_sub_pd(t0, t2)));
	    _mm_store_pd(&b[i5], ZMUL(w1, _mm_add_pd(t1, t3)));
	    _mm_store_pd(&b[i7], ZMUL(w3, _mm_sub_pd(t1, t3)));
	}
    }
    return 0;
}

int fft5a_(double *a, double *b, double *w, int *l)
{
    /* static double c51 = .95105651629515357;
    static double c52 = .61803398874989485;
    static double c53 = .55901699437494742;
    static double c54 = .25; */
    static __m128d c51, c52, c53, c54;

    int j, j0, j1, j2, j3, j4, j5, j6, j7, j8, j9;
    /* double x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7,
                x8, y8, x9, y9, x10, y10, wi1, wi2, wi3, wi4, wr1, wr2, wr3, wr4; */
    __m128d t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, w1, w2, w3, w4;

    c51 = _mm_set1_pd(0.95105651629515357);
    c52 = _mm_set1_pd(0.61803398874989485);
    c53 = _mm_set1_pd(0.55901699437494742);
    c54 = _mm_set1_pd(0.25);

    for (j = 0; j < *l; j++) {
        j0 = j << 1;
	j1 = j0 + (*l << 1);
	j2 = j1 + (*l << 1);
	j3 = j2 + (*l << 1);
	j4 = j3 + (*l << 1);
	j5 = j * 10;
	j6 = j5 + 2;
	j7 = j6 + 2;
	j8 = j7 + 2;
	j9 = j8 + 2;
	/* wr1 = w[j0];
	wi1 = w[j0 + 1];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	wr3 = wr1 * wr2 - wi1 * wi2;
	wi3 = wr1 * wi2 + wi1 * wr2;
	wr4 = wr2 * wr2 - wi2 * wi2;
	wi4 = wr2 * wi2 + wr2 * wi2; */
	w1 = _mm_load_pd(&w[j0]);
	w2 = ZMUL(w1, w1);
	w3 = ZMUL(w1, w2);
	w4 = ZMUL(w2, w2);
	/* x0 = a[j1] + a[j4];
	y0 = a[j1 + 1] + a[j4 + 1];
	x1 = a[j2] + a[j3];
	y1 = a[j2 + 1] + a[j3 + 1];
	x2 = c51 * (a[j1] - a[j4]);
	y2 = c51 * (a[j1 + 1] - a[j4 + 1]);
	x3 = c51 * (a[j2] - a[j3]);
	y3 = c51 * (a[j2 + 1] - a[j3 + 1]);
	x4 = x0 + x1;
	y4 = y0 + y1;
	x5 = c53 * (x0 - x1);
	y5 = c53 * (y0 - y1);
	x6 = a[j0] - c54 * x4;
	y6 = a[j0 + 1] - c54 * y4;
	x7 = x6 + x5;
	y7 = y6 + y5;
	x8 = x6 - x5;
	y8 = y6 - y5;
	x9 = y2 + c52 * y3;
	y9 = -x2 - c52 * x3;
	x10 = c52 * y2 - y3;
	y10 = x3 - c52 * x2; */
	t1 = _mm_load_pd(&a[j1]);
	t4 = _mm_load_pd(&a[j4]);
	t0 = _mm_add_pd(t1, t4);
	t2 = _mm_mul_pd(c51, _mm_sub_pd(t1, t4));
	t1 = _mm_load_pd(&a[j2]);
	t4 = _mm_load_pd(&a[j3]);
	t3 = _mm_mul_pd(c51, _mm_sub_pd(t1, t4));
	t1 = _mm_add_pd(t1, t4);
	t4 = _mm_add_pd(t0, t1);
	t5 = _mm_mul_pd(c53, _mm_sub_pd(t0, t1));
	t0 = _mm_load_pd(&a[j0]);
	t6 = _mm_sub_pd(t0, _mm_mul_pd(c54, t4));
	t7 = _mm_add_pd(t6, t5);
	t8 = _mm_sub_pd(t6, t5);
	t9 = _mm_xor_pd(_mm_add_pd(t2, _mm_mul_pd(c52, t3)), _mm_set_sd(-0.0));
	t9 = _mm_shuffle_pd(t9, t9, 1);
	t10 = _mm_sub_pd(t3, _mm_mul_pd(c52, t2));
	t10 = _mm_xor_pd(_mm_shuffle_pd(t10, t10, 1), _mm_set_sd(-0.0));
	/* b[j5] = a[j0] + x4;
	b[j5 + 1] = a[j0 + 1] + y4;
	b[j6] = wr1 * (x7 + x9) - wi1 * (y7 + y9);
	b[j6 + 1] = wr1 * (y7 + y9) + wi1 * (x7 + x9);
	b[j7] = wr2 * (x8 + x10) - wi2 * (y8 + y10);
	b[j7 + 1] = wr2 * (y8 + y10) + wi2 * (x8 + x10);
	b[j8] = wr3 * (x8 - x10) - wi3 * (y8 - y10);
	b[j8 + 1] = wr3 * (y8 - y10) + wi3 * (x8 - x10);
	b[j9] = wr4 * (x7 - x9) - wi4 * (y7 - y9);
	b[j9 + 1] = wr4 * (y7 - y9) + wi4 * (x7 - x9); */
	_mm_store_pd(&b[j5], _mm_add_pd(t0, t4));
	_mm_store_pd(&b[j6], ZMUL(w1, _mm_add_pd(t7, t9)));
	_mm_store_pd(&b[j7], ZMUL(w2, _mm_add_pd(t8, t10)));
	_mm_store_pd(&b[j8], ZMUL(w3, _mm_sub_pd(t8, t10)));
	_mm_store_pd(&b[j9], ZMUL(w4, _mm_sub_pd(t7, t9)));
    }
    return 0;
}

int fft5b_(double *a, double *b, double *w, int *m, int *l)
{
    /* static double c51 = .95105651629515357;
    static double c52 = .61803398874989485;
    static double c53 = .55901699437494742;
    static double c54 = .25; */
    static __m128d c51, c52, c53, c54;

    int i, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, j, j0;
    /* double x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7,
                x8, y8, x9, y9, x10, y10, wi1, wi2, wi3, wi4, wr1, wr2, wr3, wr4; */
    __m128d t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, w1, w2, w3, w4;

    c51 = _mm_set1_pd(0.95105651629515357);
    c52 = _mm_set1_pd(0.61803398874989485);
    c53 = _mm_set1_pd(0.55901699437494742);
    c54 = _mm_set1_pd(0.25);

    for (i = 0; i < *m; i++) {
        i0 = i << 1;
	i1 = i0 + (*m * *l << 1);
	i2 = i1 + (*m * *l << 1);
	i3 = i2 + (*m * *l << 1);
	i4 = i3 + (*m * *l << 1);
	i5 = i << 1;
	i6 = i5 + (*m << 1);
	i7 = i6 + (*m << 1);
	i8 = i7 + (*m << 1);
	i9 = i8 + (*m << 1);
	/* x0 = a[i1] + a[i4];
	y0 = a[i1 + 1] + a[i4 + 1];
	x1 = a[i2] + a[i3];
	y1 = a[i2 + 1] + a[i3 + 1];
	x2 = c51 * (a[i1] - a[i4]);
	y2 = c51 * (a[i1 + 1] - a[i4 + 1]);
	x3 = c51 * (a[i2] - a[i3]);
	y3 = c51 * (a[i2 + 1] - a[i3 + 1]);
	x4 = x0 + x1;
	y4 = y0 + y1;
	x5 = c53 * (x0 - x1);
	y5 = c53 * (y0 - y1);
	x6 = a[i0] - c54 * x4;
	y6 = a[i0 + 1] - c54 * y4;
	x7 = x6 + x5;
	y7 = y6 + y5;
	x8 = x6 - x5;
	y8 = y6 - y5;
	x9 = y2 + c52 * y3;
	y9 = -x2 - c52 * x3;
	x10 = c52 * y2 - y3;
	y10 = x3 - c52 * x2; */
	t1 = _mm_load_pd(&a[i1]);
	t4 = _mm_load_pd(&a[i4]);
	t0 = _mm_add_pd(t1, t4);
	t2 = _mm_mul_pd(c51, _mm_sub_pd(t1, t4));
	t1 = _mm_load_pd(&a[i2]);
	t4 = _mm_load_pd(&a[i3]);
	t3 = _mm_mul_pd(c51, _mm_sub_pd(t1, t4));
	t1 = _mm_add_pd(t1, t4);
	t4 = _mm_add_pd(t0, t1);
	t5 = _mm_mul_pd(c53, _mm_sub_pd(t0, t1));
	t0 = _mm_load_pd(&a[i0]);
	t6 = _mm_sub_pd(t0, _mm_mul_pd(c54, t4));
	t7 = _mm_add_pd(t6, t5);
	t8 = _mm_sub_pd(t6, t5);
	t9 = _mm_xor_pd(_mm_add_pd(t2, _mm_mul_pd(c52, t3)), _mm_set_sd(-0.0));
	t9 = _mm_shuffle_pd(t9, t9, 1);
	t10 = _mm_sub_pd(t3, _mm_mul_pd(c52, t2));
	t10 = _mm_xor_pd(_mm_shuffle_pd(t10, t10, 1), _mm_set_sd(-0.0));
	/* b[i5] = a[i0] + x4;
	b[i5 + 1] = a[i0 + 1] + y4;
	b[i6] = x7 + x9;
	b[i6 + 1] = y7 + y9;
	b[i7] = x8 + x10;
	b[i7 + 1] = y8 + y10;
	b[i8] = x8 - x10;
	b[i8 + 1] = y8 - y10;
	b[i9] = x7 - x9;
	b[i9 + 1] = y7 - y9; */
	_mm_store_pd(&b[i5], _mm_add_pd(t0, t4));
	_mm_store_pd(&b[i6], _mm_add_pd(t7, t9));
	_mm_store_pd(&b[i7], _mm_add_pd(t8, t10));
	_mm_store_pd(&b[i8], _mm_sub_pd(t8, t10));
	_mm_store_pd(&b[i9], _mm_sub_pd(t7, t9));
    }
    for (j = 1; j < *l; j++) {
        j0 = j << 1;
	/* wr1 = w[j0];
	wi1 = w[j0 + 1];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	wr3 = wr1 * wr2 - wi1 * wi2;
	wi3 = wr1 * wi2 + wi1 * wr2;
	wr4 = wr2 * wr2 - wi2 * wi2;
	wi4 = wr2 * wi2 + wr2 * wi2; */
	w1 = _mm_load_pd(&w[j0]);
	w2 = ZMUL(w1, w1);
	w3 = ZMUL(w1, w2);
	w4 = ZMUL(w2, w2);
	for (i = 0; i < *m; i++) {
	    i0 = (i << 1) + (j * *m << 1);
	    i1 = i0 + (*m * *l << 1);
	    i2 = i1 + (*m * *l << 1);
	    i3 = i2 + (*m * *l << 1);
	    i4 = i3 + (*m * *l << 1);
	    i5 = (i << 1) + (j * *m * 10);
	    i6 = i5 + (*m << 1);
	    i7 = i6 + (*m << 1);
	    i8 = i7 + (*m << 1);
	    i9 = i8 + (*m << 1);
	    /* x0 = a[i1] + a[i4];
	    y0 = a[i1 + 1] + a[i4 + 1];
	    x1 = a[i2] + a[i3];
	    y1 = a[i2 + 1] + a[i3 + 1];
	    x2 = c51 * (a[i1] - a[i4]);
	    y2 = c51 * (a[i1 + 1] - a[i4 + 1]);
	    x3 = c51 * (a[i2] - a[i3]);
	    y3 = c51 * (a[i2 + 1] - a[i3 + 1]);
	    x4 = x0 + x1;
	    y4 = y0 + y1;
	    x5 = c53 * (x0 - x1);
	    y5 = c53 * (y0 - y1);
	    x6 = a[i0] - c54 * x4;
	    y6 = a[i0 + 1] - c54 * y4;
	    x7 = x6 + x5;
	    y7 = y6 + y5;
	    x8 = x6 - x5;
	    y8 = y6 - y5;
	    x9 = y2 + c52 * y3;
	    y9 = -x2 - c52 * x3;
	    x10 = c52 * y2 - y3;
	    y10 = x3 - c52 * x2; */
	    t1 = _mm_load_pd(&a[i1]);
	    t4 = _mm_load_pd(&a[i4]);
	    t0 = _mm_add_pd(t1, t4);
	    t2 = _mm_mul_pd(c51, _mm_sub_pd(t1, t4));
	    t1 = _mm_load_pd(&a[i2]);
	    t4 = _mm_load_pd(&a[i3]);
	    t3 = _mm_mul_pd(c51, _mm_sub_pd(t1, t4));
	    t1 = _mm_add_pd(t1, t4);
	    t4 = _mm_add_pd(t0, t1);
	    t5 = _mm_mul_pd(c53, _mm_sub_pd(t0, t1));
	    t0 = _mm_load_pd(&a[i0]);
	    t6 = _mm_sub_pd(t0, _mm_mul_pd(c54, t4));
	    t7 = _mm_add_pd(t6, t5);
	    t8 = _mm_sub_pd(t6, t5);
	    t9 = _mm_xor_pd(_mm_add_pd(t2, _mm_mul_pd(c52, t3)), _mm_set_sd(-0.0));
	    t9 = _mm_shuffle_pd(t9, t9, 1);
	    t10 = _mm_sub_pd(t3, _mm_mul_pd(c52, t2));
	    t10 = _mm_xor_pd(_mm_shuffle_pd(t10, t10, 1), _mm_set_sd(-0.0));
	    /* b[i5] = a[i0] + x4;
	    b[i5 + 1] = a[i0 + 1] + y4;
	    b[i6] = wr1 * (x7 + x9) - wi1 * (y7 + y9);
	    b[i6 + 1] = wr1 * (y7 + y9) + wi1 * (x7 + x9);
	    b[i7] = wr2 * (x8 + x10) - wi2 * (y8 + y10);
	    b[i7 + 1] = wr2 * (y8 + y10) + wi2 * (x8 + x10);
	    b[i8] = wr3 * (x8 - x10) - wi3 * (y8 - y10);
	    b[i8 + 1] = wr3 * (y8 - y10) + wi3 * (x8 - x10);
	    b[i9] = wr4 * (x7 - x9) - wi4 * (y7 - y9);
	    b[i9 + 1] = wr4 * (y7 - y9) + wi4 * (x7 - x9); */
	    _mm_store_pd(&b[i5], _mm_add_pd(t0, t4));
	    _mm_store_pd(&b[i6], ZMUL(w1, _mm_add_pd(t7, t9)));
	    _mm_store_pd(&b[i7], ZMUL(w2, _mm_add_pd(t8, t10)));
	    _mm_store_pd(&b[i8], ZMUL(w3, _mm_sub_pd(t8, t10)));
	    _mm_store_pd(&b[i9], ZMUL(w4, _mm_sub_pd(t7, t9)));
	}
    }
    return 0;
}

int fft8a_(double *a, double *b, double *w, int *l)
{
    /* static double c81 = .70710678118654752; */
    static __m128d c81;

    int j, j0, j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14, j15;
    /* double u0, v0, u1, x0, y0, x1, y1, x2, y2, x3, y3, v1, x4, y4, x5, y5,
             x6, y6, x7, y7, u2, v2, u3, v3, wi1, wi2, wi3, wi4, wi5, wi6,
             wi7, wr1, wr2, wr3, wr4, wr5, wr6, wr7; */
    __m128d t0, t1, t2, t3, t4, t5, t6, t7, t8, u0, u1, u2, u3, w1, w2, w3, w4, w5, w6, w7;

    c81 = _mm_set1_pd(0.70710678118654752);

    for (j = 0; j < *l; j++) {
        j0 = j << 1;
        j1 = j0 + (*l << 1);
        j2 = j1 + (*l << 1);
        j3 = j2 + (*l << 1);
        j4 = j3 + (*l << 1);
        j5 = j4 + (*l << 1);
        j6 = j5 + (*l << 1);
        j7 = j6 + (*l << 1);
        j8 = j << 4;
        j9 = j8 + 2;
        j10 = j9 + 2;
        j11 = j10 + 2;
        j12 = j11 + 2;
        j13 = j12 + 2;
        j14 = j13 + 2;
        j15 = j14 + 2;
	/* wr1 = w[j0];
	wi1 = w[j0 + 1];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	wr3 = wr1 * wr2 - wi1 * wi2;
	wi3 = wr1 * wi2 + wi1 * wr2;
	wr4 = wr2 * wr2 - wi2 * wi2;
	wi4 = wr2 * wi2 + wr2 * wi2;
	wr5 = wr2 * wr3 - wi2 * wi3;
	wi5 = wr2 * wi3 + wi2 * wr3;
	wr6 = wr3 * wr3 - wi3 * wi3;
	wi6 = wr3 * wi3 + wr3 * wi3;
	wr7 = wr3 * wr4 - wi3 * wi4;
	wi7 = wr3 * wi4 + wi3 * wr4; */
	w1 = _mm_load_pd(&w[j0]);
	w2 = ZMUL(w1, w1);
	w3 = ZMUL(w1, w2);
	w4 = ZMUL(w2, w2);
	w5 = ZMUL(w2, w3);
	w6 = ZMUL(w3, w3);
	w7 = ZMUL(w3, w4);
	/* x0 = a[j0] + a[j4];
	y0 = a[j0 + 1] + a[j4 + 1];
	x1 = a[j0] - a[j4];
	y1 = a[j0 + 1] - a[j4 + 1];
	x2 = a[j2] + a[j6];
	y2 = a[j2 + 1] + a[j6 + 1];
	x3 = a[j2 + 1] - a[j6 + 1];
	y3 = a[j6] - a[j2]; */
	t0 = _mm_load_pd(&a[j0]);
	t2 = _mm_load_pd(&a[j4]);
	t1 = _mm_sub_pd(t0, t2);
	t0 = _mm_add_pd(t0, t2);
	t3 = _mm_load_pd(&a[j2]);
	t4 = _mm_load_pd(&a[j6]);
	t2 = _mm_add_pd(t3, t4);
	t3 = _mm_xor_pd(_mm_sub_pd(t3, t4), _mm_set_sd(-0.0));
	t3 = _mm_shuffle_pd(t3, t3, 1);
	/* u0 = x0 + x2;
	v0 = y0 + y2;
	u1 = x0 - x2;
	v1 = y0 - y2; */
	u0 = _mm_add_pd(t0, t2);
	u1 = _mm_sub_pd(t0, t2);
	/* x4 = a[j1] + a[j5];
	y4 = a[j1 + 1] + a[j5 + 1];
	x5 = a[j1] - a[j5];
	y5 = a[j1 + 1] - a[j5 + 1];
	x6 = a[j3] + a[j7];
	y6 = a[j3 + 1] + a[j7 + 1];
	x7 = a[j3] - a[j7];
	y7 = a[j3 + 1] - a[j7 + 1]; */
	t4 = _mm_load_pd(&a[j1]);
	t6 = _mm_load_pd(&a[j5]);
	t5 = _mm_sub_pd(t4, t6);
	t4 = _mm_add_pd(t4, t6);
	t7 = _mm_load_pd(&a[j3]);
	t8 = _mm_load_pd(&a[j7]);
	t6 = _mm_add_pd(t7, t8);
	t7 = _mm_sub_pd(t7, t8);
	/* u2 = x4 + x6;
	v2 = y4 + y6;
	u3 = y4 - y6;
	v3 = x6 - x4; */
	u2 = _mm_add_pd(t4, t6);
	u3 = _mm_xor_pd(_mm_sub_pd(t4, t6), _mm_set_sd(-0.0));
	u3 = _mm_shuffle_pd(u3, u3, 1);
	/* b[j8] = u0 + u2;
	b[j8 + 1] = v0 + v2;
	b[j12] = wr4 * (u0 - u2) - wi4 * (v0 - v2);
	b[j12 + 1] = wr4 * (v0 - v2) + wi4 * (u0 - u2);
	b[j10] = wr2 * (u1 + u3) - wi2 * (v1 + v3);
	b[j10 + 1] = wr2 * (v1 + v3) + wi2 * (u1 + u3);
	b[j14] = wr6 * (u1 - u3) - wi6 * (v1 - v3);
	b[j14 + 1] = wr6 * (v1 - v3) + wi6 * (u1 - u3); */
	_mm_store_pd(&b[j8], _mm_add_pd(u0, u2));
	_mm_store_pd(&b[j12], ZMUL(w4, _mm_sub_pd(u0, u2)));
	_mm_store_pd(&b[j10], ZMUL(w2, _mm_add_pd(u1, u3)));
	_mm_store_pd(&b[j14], ZMUL(w6, _mm_sub_pd(u1, u3)));
	/* u0 = x1 + c81 * (x5 - x7);
	v0 = y1 + c81 * (y5 - y7);
	u1 = x1 - c81 * (x5 - x7);
	v1 = y1 - c81 * (y5 - y7);
	u2 = x3 + c81 * (y5 + y7);
	v2 = y3 - c81 * (x5 + x7);
	u3 = x3 - c81 * (y5 + y7);
	v3 = y3 + c81 * (x5 + x7); */
	u1 = _mm_mul_pd(c81, _mm_sub_pd(t5, t7));
	u0 = _mm_add_pd(t1, u1);
	u1 = _mm_sub_pd(t1, u1);
	u3 = _mm_xor_pd(_mm_mul_pd(c81, _mm_add_pd(t5, t7)), _mm_set_sd(-0.0));
	u3 = _mm_shuffle_pd(u3, u3, 1);
	u2 = _mm_add_pd(t3, u3);
	u3 = _mm_sub_pd(t3, u3);
	/* b[j9] = wr1 * (u0 + u2) - wi1 * (v0 + v2);
	b[j9 + 1] = wr1 * (v0 + v2) + wi1 * (u0 + u2);
	b[j13] = wr5 * (u1 + u3) - wi5 * (v1 + v3);
	b[j13 + 1] = wr5 * (v1 + v3) + wi5 * (u1 + u3);
	b[j11] = wr3 * (u1 - u3) - wi3 * (v1 - v3);
	b[j11 + 1] = wr3 * (v1 - v3) + wi3 * (u1 - u3);
	b[j15] = wr7 * (u0 - u2) - wi7 * (v0 - v2);
	b[j15 + 1] = wr7 * (v0 - v2) + wi7 * (u0 - u2); */
	_mm_store_pd(&b[j9], ZMUL(w1, _mm_add_pd(u0, u2)));
	_mm_store_pd(&b[j13], ZMUL(w5, _mm_add_pd(u1, u3)));
	_mm_store_pd(&b[j11], ZMUL(w3, _mm_sub_pd(u1, u3)));
	_mm_store_pd(&b[j15], ZMUL(w7, _mm_sub_pd(u0, u2)));
    }
    return 0;
}

int fft8b_(double *a, double *b, double *w, int *m, int *l)
{
    /* static double c81 = .70710678118654752; */
    static __m128d c81;

    int i, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, j, j0;
    /* double u0, v0, u1, x0, y0, x1, y1, x2, y2, x3, y3, v1, x4, y4, x5, y5,
             x6, y6, x7, y7, u2, v2, u3, v3, wi1, wi2, wi3, wi4, wi5, wi6,
             wi7, wr1, wr2, wr3, wr4, wr5, wr6, wr7; */
    __m128d t0, t1, t2, t3, t4, t5, t6, t7, t8, u0, u1, u2, u3, w1, w2, w3, w4, w5, w6, w7;

    c81 = _mm_set1_pd(0.70710678118654752);

    for (i = 0; i < *m; i++) {
        i0 = i << 1;
	i1 = i0 + (*m * *l << 1);
	i2 = i1 + (*m * *l << 1);
	i3 = i2 + (*m * *l << 1);
	i4 = i3 + (*m * *l << 1);
	i5 = i4 + (*m * *l << 1);
	i6 = i5 + (*m * *l << 1);
	i7 = i6 + (*m * *l << 1);
	i8 = i << 1;
	i9 = i8 + (*m << 1);
	i10 = i9 + (*m << 1);
	i11 = i10 + (*m << 1);
	i12 = i11 + (*m << 1);
	i13 = i12 + (*m << 1);
	i14 = i13 + (*m << 1);
	i15 = i14 + (*m << 1);
	/* x0 = a[i0] + a[i4];
	y0 = a[i0 + 1] + a[i4 + 1];
	x1 = a[i0] - a[i4];
	y1 = a[i0 + 1] - a[i4 + 1];
	x2 = a[i2] + a[i6];
	y2 = a[i2 + 1] + a[i6 + 1];
	x3 = a[i2 + 1] - a[i6 + 1];
	y3 = a[i6] - a[i2]; */
	t0 = _mm_load_pd(&a[i0]);
	t2 = _mm_load_pd(&a[i4]);
	t1 = _mm_sub_pd(t0, t2);
	t0 = _mm_add_pd(t0, t2);
	t3 = _mm_load_pd(&a[i2]);
	t4 = _mm_load_pd(&a[i6]);
	t2 = _mm_add_pd(t3, t4);
	t3 = _mm_xor_pd(_mm_sub_pd(t3, t4), _mm_set_sd(-0.0));
	t3 = _mm_shuffle_pd(t3, t3, 1);
	/* u0 = x0 + x2;
	v0 = y0 + y2;
	u1 = x0 - x2;
	v1 = y0 - y2; */
	u0 = _mm_add_pd(t0, t2);
	u1 = _mm_sub_pd(t0, t2);
	/* x4 = a[i1] + a[i5];
	y4 = a[i1 + 1] + a[i5 + 1];
	x5 = a[i1] - a[i5];
	y5 = a[i1 + 1] - a[i5 + 1];
	x6 = a[i3] + a[i7];
	y6 = a[i3 + 1] + a[i7 + 1];
	x7 = a[i3] - a[i7];
	y7 = a[i3 + 1] - a[i7 + 1]; */
	t4 = _mm_load_pd(&a[i1]);
	t6 = _mm_load_pd(&a[i5]);
	t5 = _mm_sub_pd(t4, t6);
	t4 = _mm_add_pd(t4, t6);
	t7 = _mm_load_pd(&a[i3]);
	t8 = _mm_load_pd(&a[i7]);
	t6 = _mm_add_pd(t7, t8);
	t7 = _mm_sub_pd(t7, t8);
	/* u2 = x4 + x6;
	v2 = y4 + y6;
	u3 = y4 - y6;
	v3 = x6 - x4; */
	u2 = _mm_add_pd(t4, t6);
	u3 = _mm_xor_pd(_mm_sub_pd(t4, t6), _mm_set_sd(-0.0));
	u3 = _mm_shuffle_pd(u3, u3, 1);
	/* b[i8] = u0 + u2;
	b[i8 + 1] = v0 + v2;
	b[i12] = u0 - u2;
	b[i12 + 1] = v0 - v2;
	b[i10] = u1 + u3;
	b[i10 + 1] = v1 + v3;
	b[i14] = u1 - u3;
	b[i14 + 1] = v1 - v3; */
	_mm_store_pd(&b[i8], _mm_add_pd(u0, u2));
	_mm_store_pd(&b[i12], _mm_sub_pd(u0, u2));
	_mm_store_pd(&b[i10], _mm_add_pd(u1, u3));
	_mm_store_pd(&b[i14], _mm_sub_pd(u1, u3));
	/* u0 = x1 + c81 * (x5 - x7);
	v0 = y1 + c81 * (y5 - y7);
	u1 = x1 - c81 * (x5 - x7);
	v1 = y1 - c81 * (y5 - y7);
	u2 = x3 + c81 * (y5 + y7);
	v2 = y3 - c81 * (x5 + x7);
	u3 = x3 - c81 * (y5 + y7);
	v3 = y3 + c81 * (x5 + x7); */
	u1 = _mm_mul_pd(c81, _mm_sub_pd(t5, t7));
	u0 = _mm_add_pd(t1, u1);
	u1 = _mm_sub_pd(t1, u1);
	u3 = _mm_xor_pd(_mm_mul_pd(c81, _mm_add_pd(t5, t7)), _mm_set_sd(-0.0));
	u3 = _mm_shuffle_pd(u3, u3, 1);
	u2 = _mm_add_pd(t3, u3);
	u3 = _mm_sub_pd(t3, u3);
	/* b[i9] = u0 + u2;
	b[i9 + 1] = v0 + v2;
	b[i13] = u1 + u3;
	b[i13 + 1] = v1 + v3;
	b[i11] = u1 - u3;
	b[i11 + 1] = v1 - v3;
	b[i15] = u0 - u2;
	b[i15 + 1] = v0 - v2; */
	_mm_store_pd(&b[i9], _mm_add_pd(u0, u2));
	_mm_store_pd(&b[i13], _mm_add_pd(u1, u3));
	_mm_store_pd(&b[i11], _mm_sub_pd(u1, u3));
	_mm_store_pd(&b[i15], _mm_sub_pd(u0, u2));
    }
    for (j = 1; j < *l; j++) {
        j0 = j << 1;
	/* wr1 = w[j0];
	wi1 = w[j0 + 1];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	wr3 = wr1 * wr2 - wi1 * wi2;
	wi3 = wr1 * wi2 + wi1 * wr2;
	wr4 = wr2 * wr2 - wi2 * wi2;
	wi4 = wr2 * wi2 + wr2 * wi2;
	wr5 = wr2 * wr3 - wi2 * wi3;
	wi5 = wr2 * wi3 + wi2 * wr3;
	wr6 = wr3 * wr3 - wi3 * wi3;
	wi6 = wr3 * wi3 + wr3 * wi3;
	wr7 = wr3 * wr4 - wi3 * wi4;
	wi7 = wr3 * wi4 + wi3 * wr4; */
	w1 = _mm_load_pd(&w[j0]);
	w2 = ZMUL(w1, w1);
	w3 = ZMUL(w1, w2);
	w4 = ZMUL(w2, w2);
	w5 = ZMUL(w2, w3);
	w6 = ZMUL(w3, w3);
	w7 = ZMUL(w3, w4);
	for (i = 0; i < *m; i++) {
	    i0 = (i << 1) + (j * *m << 1);
	    i1 = i0 + (*m * *l << 1);
	    i2 = i1 + (*m * *l << 1);
	    i3 = i2 + (*m * *l << 1);
	    i4 = i3 + (*m * *l << 1);
	    i5 = i4 + (*m * *l << 1);
	    i6 = i5 + (*m * *l << 1);
	    i7 = i6 + (*m * *l << 1);
	    i8 = (i << 1) + (j * *m << 4);
	    i9 = i8 + (*m << 1);
	    i10 = i9 + (*m << 1);
	    i11 = i10 + (*m << 1);
	    i12 = i11 + (*m << 1);
	    i13 = i12 + (*m << 1);
	    i14 = i13 + (*m << 1);
	    i15 = i14 + (*m << 1);
	    /* x0 = a[i0] + a[i4];
	    y0 = a[i0 + 1] + a[i4 + 1];
	    x1 = a[i0] - a[i4];
	    y1 = a[i0 + 1] - a[i4 + 1];
	    x2 = a[i2] + a[i6];
	    y2 = a[i2 + 1] + a[i6 + 1];
	    x3 = a[i2 + 1] - a[i6 + 1];
	    y3 = a[i6] - a[i2]; */
	    t0 = _mm_load_pd(&a[i0]);
	    t2 = _mm_load_pd(&a[i4]);
	    t1 = _mm_sub_pd(t0, t2);
	    t0 = _mm_add_pd(t0, t2);
	    t3 = _mm_load_pd(&a[i2]);
	    t4 = _mm_load_pd(&a[i6]);
	    t2 = _mm_add_pd(t3, t4);
	    t3 = _mm_xor_pd(_mm_sub_pd(t3, t4), _mm_set_sd(-0.0));
	    t3 = _mm_shuffle_pd(t3, t3, 1);
	    /* u0 = x0 + x2;
	    v0 = y0 + y2;
	    u1 = x0 - x2;
	    v1 = y0 - y2; */
	    u0 = _mm_add_pd(t0, t2);
	    u1 = _mm_sub_pd(t0, t2);
	    /* x4 = a[i1] + a[i5];
	    y4 = a[i1 + 1] + a[i5 + 1];
	    x5 = a[i1] - a[i5];
	    y5 = a[i1 + 1] - a[i5 + 1];
	    x6 = a[i3] + a[i7];
	    y6 = a[i3 + 1] + a[i7 + 1];
	    x7 = a[i3] - a[i7];
	    y7 = a[i3 + 1] - a[i7 + 1]; */
	    t4 = _mm_load_pd(&a[i1]);
	    t6 = _mm_load_pd(&a[i5]);
	    t5 = _mm_sub_pd(t4, t6);
	    t4 = _mm_add_pd(t4, t6);
	    t7 = _mm_load_pd(&a[i3]);
	    t8 = _mm_load_pd(&a[i7]);
	    t6 = _mm_add_pd(t7, t8);
	    t7 = _mm_sub_pd(t7, t8);
	    /* u2 = x4 + x6;
	    v2 = y4 + y6;
	    u3 = y4 - y6;
	    v3 = x6 - x4; */
	    u2 = _mm_add_pd(t4, t6);
	    u3 = _mm_xor_pd(_mm_sub_pd(t4, t6), _mm_set_sd(-0.0));
	    u3 = _mm_shuffle_pd(u3, u3, 1);
	    /* b[i8] = u0 + u2;
	    b[i8 + 1] = v0 + v2;
	    b[i12] = wr4 * (u0 - u2) - wi4 * (v0 - v2);
	    b[i12 + 1] = wr4 * (v0 - v2) + wi4 * (u0 - u2);
	    b[i10] = wr2 * (u1 + u3) - wi2 * (v1 + v3);
	    b[i10 + 1] = wr2 * (v1 + v3) + wi2 * (u1 + u3);
	    b[i14] = wr6 * (u1 - u3) - wi6 * (v1 - v3);
	    b[i14 + 1] = wr6 * (v1 - v3) + wi6 * (u1 - u3); */
	    _mm_store_pd(&b[i8], _mm_add_pd(u0, u2));
	    _mm_store_pd(&b[i12], ZMUL(w4, _mm_sub_pd(u0, u2)));
	    _mm_store_pd(&b[i10], ZMUL(w2, _mm_add_pd(u1, u3)));
	    _mm_store_pd(&b[i14], ZMUL(w6, _mm_sub_pd(u1, u3)));
	    /* u0 = x1 + c81 * (x5 - x7);
	    v0 = y1 + c81 * (y5 - y7);
	    u1 = x1 - c81 * (x5 - x7);
	    v1 = y1 - c81 * (y5 - y7);
	    u2 = x3 + c81 * (y5 + y7);
	    v2 = y3 - c81 * (x5 + x7);
	    u3 = x3 - c81 * (y5 + y7);
	    v3 = y3 + c81 * (x5 + x7); */
	    u1 = _mm_mul_pd(c81, _mm_sub_pd(t5, t7));
	    u0 = _mm_add_pd(t1, u1);
	    u1 = _mm_sub_pd(t1, u1);
	    u3 = _mm_xor_pd(_mm_mul_pd(c81, _mm_add_pd(t5, t7)), _mm_set_sd(-0.0));
	    u3 = _mm_shuffle_pd(u3, u3, 1);
	    u2 = _mm_add_pd(t3, u3);
	    u3 = _mm_sub_pd(t3, u3);
	    /* b[i9] = wr1 * (u0 + u2) - wi1 * (v0 + v2);
	    b[i9 + 1] = wr1 * (v0 + v2) + wi1 * (u0 + u2);
	    b[i13] = wr5 * (u1 + u3) - wi5 * (v1 + v3);
	    b[i13 + 1] = wr5 * (v1 + v3) + wi5 * (u1 + u3);
	    b[i11] = wr3 * (u1 - u3) - wi3 * (v1 - v3);
	    b[i11 + 1] = wr3 * (v1 - v3) + wi3 * (u1 - u3);
	    b[i15] = wr7 * (u0 - u2) - wi7 * (v0 - v2);
	    b[i15 + 1] = wr7 * (v0 - v2) + wi7 * (u0 - u2); */
	    _mm_store_pd(&b[i9], ZMUL(w1, _mm_add_pd(u0, u2)));
	    _mm_store_pd(&b[i13], ZMUL(w5, _mm_add_pd(u1, u3)));
	    _mm_store_pd(&b[i11], ZMUL(w3, _mm_sub_pd(u1, u3)));
	    _mm_store_pd(&b[i15], ZMUL(w7, _mm_sub_pd(u0, u2)));
	}
    }
    return 0;
}
