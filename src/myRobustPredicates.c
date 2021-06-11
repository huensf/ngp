#include "myRobustPredicates.h"

static double splitter, epsilon;
static double resulterrbound;
static double ccwerrboundA, ccwerrboundB, ccwerrboundC;

void my_exactinit() {

#ifdef CPU86
#ifdef SINGLE
    _control87(_PC_24, _MCW_PC); /* Set FPU control word for single precision. */
#else /* not SINGLE */
    _control87(_PC_53, _MCW_PC); /* Set FPU control word for double precision. */
#endif /* not SINGLE */
#endif /* CPU86 */

    double lastcheck;
    int every_other = 1;
    double half = 0.5;
    epsilon = 1.0;
    splitter = 1.0;
    double check = 1.0;

    do {
        lastcheck = check;
        epsilon *= half;
        if (every_other) 
            splitter *= 2.0;
        every_other = !every_other;
        check = 1.0 + epsilon;
    } while ((check != 1.0) && (check != lastcheck));
    splitter += 1.0;

    resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
    ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon;
    ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon;
    ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;

    return;
}

double my_estimate(int elen, double* e) {
    double Q;
    int eindex;

    Q = e[0];
    for (eindex = 1; eindex < elen; eindex++)
        Q += e[eindex];
    return Q;
}

double my_Absolute(double a) {
    return (a >= 0.0 ? a : (-a));
}

void my_Fast_Two_Sum_Tail(double a, double b, double x, double y) {
    double bvirt = x - a;
    y = b - bvirt;
}

void my_Fast_Two_Sum(double a, double b, double x, double y) {
    x = (double)(a + b);
    my_Fast_Two_Sum_Tail(a, b, x, y);
}

void my_Two_Sum_Tail(double a, double b, double x, double y) {
    double bvirt = (double)(x - a);
    double avirt = x - bvirt;
    double bround = b- bvirt;
    double around = a - avirt;
    y = around + bround;
}

void my_Two_Sum(double a, double b, double x, double y) {
    x = (double)(a + b);
    my_Two_Sum_Tail(a, b, x, y);
}

void my_Two_Diff_Tail(double a, double b, double x, double y) {
    double bvirt = (double)(a - x);
    double avirt = x + bvirt;
    double bround = bvirt - b;
    double around = a - avirt;
    y = around + bround;
}

void my_Two_Diff(double a, double b, double x, double y) {
    x = (double)(a - b);
    my_Two_Diff_Tail(a, b, x, y);
}

void my_Two_One_Sum(double a1, double a0, double b,
    double x2, double x1, double x0) {
    double _i = 0;
    my_Two_Sum(a0, b, _i, x0);
    my_Two_Sum(a1, _i, x2, x1);
}

void my_Two_One_Diff(double a1, double a0, double b,
    double x2, double x1, double x0) {
    double _i = 0;
    my_Two_Diff(a0, b, _i, x0);
    my_Two_Sum(a1, _i, x2, x1);
}

void my_Two_Two_Sum(double a1, double a0, double b1, double b0,
    double x3, double x2, double x1, double x0) {
    double _j = 0, _0 = 0;
    my_Two_One_Sum(a1, a0, b0, _j, _0, x0);
    my_Two_One_Sum(_j, _0, b1, x3, x2, x1);
}

void my_Two_Two_Diff(double a1, double a0, double b1, double b0,
    double x3, double x2, double x1, double x0) {
    double _j = 0, _0 = 0;
    my_Two_One_Diff(a1, a0, b0, _j, _0, x0);
    my_Two_One_Diff(_j, _0, b1, x3, x2, x1);
}

void my_Split(double a, double ahi, double alo) {
    double c = (double)(splitter * a);
    double abig = (double)(c - a);
    ahi = c - abig;
    alo = a - ahi;
}

void my_Two_Product_Tail(double a, double b, double x, double y) {
    double ahi = 0, alo = 0, bhi = 0, blo = 0;
    my_Split(a, ahi, alo);
    my_Split(b, bhi, blo);
    double err1 = x - (ahi * bhi);
    double err2 = err1 - (alo * bhi);
    double err3 = err2 - (ahi * blo);
    y = (alo * blo) - err3;
}

void my_Two_Product(double a, double b, double x, double y) {
    x = (double)(a * b);
    my_Two_Product_Tail(a, b, x, y);
}

int my_fast_expansion_sum_zeroelim(int elen, double* e, int flen, double* f, double* h) {
    double Q;
    double Qnew = 0;
    double hh;
    double bvirt;
    double avirt, bround, around;
    int eindex, findex, hindex;
    double enow, fnow;

    enow = e[0];
    fnow = f[0];
    eindex = findex = 0;
    if ((fnow > enow) == (fnow > -enow)) {
        Q = enow;
        enow = e[++eindex];
    }
    else {
        Q = fnow;
        fnow = f[++findex];
    }
    hindex = 0;
    if ((eindex < elen) && (findex < flen)) {
        if ((fnow > enow) == (fnow > -enow)) {
            my_Fast_Two_Sum(enow, Q, Qnew, hh);
            enow = e[++eindex];
        }
        else {
            my_Fast_Two_Sum(fnow, Q, Qnew, hh);
            fnow = f[++findex];
        }
        Q = Qnew;
        if (hh != 0.0) {
            h[hindex++] = hh;
        }
        while ((eindex < elen) && (findex < flen)) {
            if ((fnow > enow) == (fnow > -enow)) {
                my_Two_Sum(Q, enow, Qnew, hh);
                enow = e[++eindex];
            }
            else {
                my_Two_Sum(Q, fnow, Qnew, hh);
                fnow = f[++findex];
            }
            Q = Qnew;
            if (hh != 0.0) {
                h[hindex++] = hh;
            }
        }
    }
    while (eindex < elen) {
        my_Two_Sum(Q, enow, Qnew, hh);
        enow = e[++eindex];
        Q = Qnew;
        if (hh != 0.0) {
            h[hindex++] = hh;
        }
    }
    while (findex < flen) {
        my_Two_Sum(Q, fnow, Qnew, hh);
        fnow = f[++findex];
        Q = Qnew;
        if (hh != 0.0) {
            h[hindex++] = hh;
        }
    }
    if ((Q != 0.0) || (hindex == 0)) {
        h[hindex++] = Q;
    }
    return hindex;
}

double my_orient2dadapt(double* pa, double* pb, double* pc, double detsum) {
    double acx, acy, bcx, bcy;
    double acxtail = 0, acytail = 0, bcxtail = 0, bcytail = 0;
    double detleft = 0, detright = 0;
    double detlefttail = 0, detrighttail = 0;
    double det, errbound;
    double B[4], C1[8], C2[12], D[16];
    double B3 = 0;
    int C1length, C2length, Dlength;
    double u[4];
    double u3;
    double s1, t1;
    double s0, t0;

    double bvirt;
    double avirt, bround, around;
    double c;
    double abig;
    double ahi, alo, bhi, blo;
    double err1, err2, err3;
    double _i, _j;
    double _0;

    acx = (double)(pa[0] - pc[0]);
    bcx = (double)(pb[0] - pc[0]);
    acy = (double)(pa[1] - pc[1]);
    bcy = (double)(pb[1] - pc[1]);

    my_Two_Product(acx, bcy, detleft, detlefttail);
    my_Two_Product(acy, bcx, detright, detrighttail);

    my_Two_Two_Diff(detleft, detlefttail, detright, detrighttail,
        B3, B[2], B[1], B[0]);
    B[3] = B3;

    det = my_estimate(4, B);
    errbound = ccwerrboundB * detsum;
    if ((det >= errbound) || (-det >= errbound)) {
        return det;
    }

    my_Two_Diff_Tail(pa[0], pc[0], acx, acxtail);
    my_Two_Diff_Tail(pb[0], pc[0], bcx, bcxtail);
    my_Two_Diff_Tail(pa[1], pc[1], acy, acytail);
    my_Two_Diff_Tail(pb[1], pc[1], bcy, bcytail);

    if ((acxtail == 0.0) && (acytail == 0.0)
        && (bcxtail == 0.0) && (bcytail == 0.0)) {
        return det;
    }

    errbound = ccwerrboundC * detsum + resulterrbound * my_Absolute(det);
    det += (acx * bcytail + bcy * acxtail)
        - (acy * bcxtail + bcx * acytail);
    if ((det >= errbound) || (-det >= errbound)) {
        return det;
    }

    my_Two_Product(acxtail, bcy, s1, s0);
    my_Two_Product(acytail, bcx, t1, t0);
    my_Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3] = u3;
    C1length = my_fast_expansion_sum_zeroelim(4, B, 4, u, C1);

    my_Two_Product(acx, bcytail, s1, s0);
    my_Two_Product(acy, bcxtail, t1, t0);
    my_Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3] = u3;
    C2length = my_fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2);

    my_Two_Product(acxtail, bcytail, s1, s0);
    my_Two_Product(acytail, bcxtail, t1, t0);
    my_Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3] = u3;
    Dlength = my_fast_expansion_sum_zeroelim(C2length, C2, 4, u, D);

    return(D[Dlength - 1]);
}

double my_orient2d(double* pa, double* pb, double* pc) {
    
    double detleft, detright, det;
    double detsum, errbound;

    detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
    detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
    det = detleft - detright;

    if (detleft > 0.0) {
        if (detright <= 0.0) {
            return det;
        }
        else {
            detsum = detleft + detright;
        }
    }
    else if (detleft < 0.0) {
        if (detright >= 0.0) {
            return det;
        }
        else {
            detsum = -detleft - detright;
        }
    }
    else {
        return det;
    }

    errbound = ccwerrboundA * detsum;
    if ((det >= errbound) || (-det >= errbound)) {
        return det;
    }

    return my_orient2dadapt(pa, pb, pc, detsum);
}
