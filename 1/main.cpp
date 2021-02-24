#include <iostream>


const unsigned N = 4;
const double EPS = 1E-6;
const unsigned WIDTH = 6;
const unsigned PRECISION = 3;


void print(double a[N][N]) {
    std::cout.precision(PRECISION);
    std::cout.setf(std::ios::fixed);
    
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < N; j++) {
            std::cout.width(WIDTH);
            std::cout << a[j][i] << ' ';
        }
        
        std::cout << std::endl;
    }
    
    std::cout << std::endl;
}


int sign(double x) {
    return x > 0? 1 : -1;
}


void swap_str(double a[N][N], unsigned str1, unsigned str2) {
    for (unsigned j = 0; j < N; j++)
        std::swap(a[j][str1], a[j][str2]);
}


void sort(double a[N][N], unsigned transp[N]) {
    unsigned last_swap = N - 1;
    for (unsigned i = 0; i < N && i != last_swap; i = (i + 1) % N)
        if (i != transp[i]) {
            std::swap(transp[i], transp[transp[i]]);
            swap_str(a, i, transp[i]);
            last_swap = i;
    }
}


void up_triangularize(double a[N][N], unsigned transp[N]) {
    for (unsigned i = 0; i < N; i++)
        transp[i] = i;
    
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
            
            std::swap(transp[j], transp[max_len_col]);
            std::swap(a[j], a[max_len_col]);
        }
        
        for (unsigned i = j + 1; i < N; i++) {
            double a_i = a[j][i];
            
            if (fabs(a_i) > EPS) {
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
    for (int i = N - 1; i >= 0; i--) {
        if (abs(a[i][i]) < EPS)
            throw std::runtime_error("det(A) = 0!");
        
        for (int j = N - 1; j >= i; j--) {
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


void de_up_triangularize(double a[N][N], unsigned transp[N]) {
    for (int j = N - 2; j >= 0; j--) {
        double q[N];
        for (unsigned i = j + 1; i < N; i++) {
            q[i] = a[j][i];
            a[j][i] = 0;
        }
        
        for (unsigned i = N - 1; i > j; i--)
            if (fabs(q[i]) > EPS) {
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
    }
    
    sort(a, transp);
}


int main() {
    double a[N][N] = {
        {2, 2, 1, 8},      // array of columns
        {4, 2, 1, 0},
        {8, 3, 1, 0},
        {1, 2, 1, 8}
    };
    unsigned transp[N];

    print(a);
    up_triangularize(a, transp);
    print(a);
    revert_tr(a);
    print(a);
    de_up_triangularize(a, transp);
    print(a);
    
    return 0;
}
