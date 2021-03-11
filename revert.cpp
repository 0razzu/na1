#include "revert.hpp"
#include <iostream>


int sign(double x) {
    return x > 0? 1 : -1;
}


void swap_str(double a[N][N], unsigned str1, unsigned str2, unsigned from) {
    for (unsigned j = from; j < N; j++)
        std::swap(a[j][str1], a[j][str2]);
}


void up_triangularize(double a[N][N], unsigned col_transp[N], unsigned str_transp[N]) {
    for (unsigned j = 0; j < N - 1; j++) {
        if (j < N - 1) {
            double max_len = 0;
            unsigned max_len_col = j;
            for (unsigned i = 0; i < N; i++)
                max_len += a[j][i] * a[j][i];
            
            for (unsigned k = j + 1; k < N; k++) {
                double len = 0;
                for (unsigned i = 0; i < N; i++)
                    len += a[k][i] * a[k][i];
                
                if (len > max_len) {
                    max_len = len;
                    max_len_col = k;
                }
            }
            
            col_transp[j] = max_len_col;
            std::swap(a[j], a[max_len_col]);
            
            double max_elem = fabs(a[j][j + 1]);
            unsigned max_elem_str = j + 1;
            for (unsigned i = j + 2; i < N; i++)
                if (fabs(a[j][i]) > max_elem) {
                    max_elem = fabs(a[j][i]);
                    max_elem_str = i;
                }
            
            str_transp[j + 1] = max_elem_str;
            swap_str(a, j + 1, max_elem_str, j);
        }
        
        for (unsigned i = j + 1; i < N; i++) {
            double a_i = a[j][i];
            
            if (a_i != 0) {
                double a_j = a[j][j];
                double z = fmax(fabs(a_j), fabs(a_i));
                double a_l = fmin(fabs(a_j), fabs(a_i)) / z;
                double den = sqrt(1 + a_l * a_l);
                double s = -a_i / z / den;
                double c = a_j / z / den;
                
                a[j][j] = c * a_j - s * a_i;
                a[j][i] = c + sign(s);
                
                for (unsigned k = j + 1; k < N; k++) {
                    double kj = a[k][j];
                    
                    a[k][j] = c * kj - s * a[k][i];
                    a[k][i] = s * kj + c * a[k][i];
                }
            }
        }
    }
}


void revert_tr(double a[N][N]) {
    for (unsigned i = N - 1; i >= 0 && i < UINT_MAX; i--) {
        for (unsigned j = N - 1; j >= i && j < UINT_MAX; j--) {
            if (i != j) {
                double ij = 0;
                for (unsigned k = i + 1; k <= j; k++)
                    ij += a[k][i] * a[j][k];
                a[j][i] = -1 / a[i][i] * ij;
            }
            
            else
                a[i][i] = 1 / a[i][i];
        }
    }
}


void de_up_triangularize(double a[N][N], unsigned col_transp[N], unsigned str_transp[N]) {
    for (unsigned j = N - 2; j >= 0 && j < UINT_MAX; j--) {
        double q[N];
        for (unsigned i = j + 1; i < N; i++) {
            q[i] = a[j][i];
            a[j][i] = 0;
        }
        
        for (unsigned i = N - 1; i > j; i--)
            if (q[i] != 0) {
                double s, c;
                
                if (q[i] < 0) {
                    c = q[i] + 1;
                    s = -sqrt(1 - c * c);
                }
                
                else {
                    c = q[i] - 1;
                    s = sqrt(1 - c * c);
                }
                
                for (unsigned k = 0; k < N; k++) {
                    double jk = a[j][k];
                    
                    a[j][k] = c * jk + s * a[i][k];
                    a[i][k] = -s * jk + c * a[i][k];
                }
            }
        
        std::swap(a[j + 1], a[str_transp[j + 1]]);
        swap_str(a, j, col_transp[j], j);
    }
}


void revert(double a[N][N]) {
    unsigned col_transp[N], str_transp[N];
    
    up_triangularize(a, col_transp, str_transp);
    revert_tr(a);
    de_up_triangularize(a, col_transp, str_transp);
}
