#include "invert.hpp"
#include <iostream>
#include <cmath>


using namespace std;


void mygen(double a[N][N], double a_inv[N][N], double alpha, double beta, int sign_law, int lambda_law, int type, int schema);
void Q_matrix(double Q[N][N], int schema);
void matr_mul(double a[N][N], double b[N][N], double c[N][N]);
double matr_inf_norm(double a[N][N]);


int main(int argc, const char* argv[]) { // α β a/b ζ_max y/i/j
    if (argc == 6) {
        double alpha = atof(argv[1]), beta = atof(argv[2]), zeta_max = atof(argv[4]);
        char changing_border = argv[3][0], type = argv[5][0];
        
        switch (type) {
            case 'y':
                cout << "Симметричная матрица,,,,,,,,";
                break;
            
            case 'i':
                cout << "Матрица простой структуры,,,,,,,,";
                break;
                
            case 'j':
                cout << "Матрица с жордановой клеткой,,,,,,,,";
                break;
                    
            default:
                return 1;
        }

        cout << "N = " << N << endl;
        cout << "α,β,||A||,||A^-1||,ν(A),||Z||,ζ,||R||,||R_gen||" << endl;

        double zeta = 0;

        while (zeta < zeta_max) {
            cout << alpha << ',';
            cout << beta << ',';

            double a[N][N], a_inv[N][N];
            switch (type) {
                case 'y':
                    mygen(a, a_inv, alpha, beta, 1, 2, 0, 1);
                    break;
                
                case 'i':
                    mygen(a, a_inv, alpha, beta, 1, 2, 1, 1);
                    break;
                
                case 'j':
                    mygen(a, a_inv, alpha, beta, 0, 0, 2, 1);
                    break;
                    
                default:
                    break;
            }

            double a_norm = matr_inf_norm(a);
            cout << a_norm << ',';

            double a_inv_norm = matr_inf_norm(a_inv);
            cout << a_inv_norm << ',';

            // ν(A)
            cout << a_norm * a_inv_norm << ',';

            double a_inv_c[N][N];
            for (int j = 0; j < N; j++)
                for (int i = 0; i < N; i++)
                    a_inv_c[j][i] = a[j][i];
            
            invert(a_inv_c);
            
            // ||Z||
            double z[N][N];
            for (int j = 0; j < N; j++)
                for (int i = 0; i < N; i++)
                    z[j][i] = a_inv_c[j][i] - a_inv[j][i];
            
            double z_norm = matr_inf_norm(z);
            cout << z_norm << ',';

            // ζ
            zeta = z_norm / matr_inf_norm(a_inv);
            cout << zeta << ',';

            // ||R||
            double r[N][N];
            matr_mul(a, a_inv_c, r);
            for (int i = 0; i < N; i++)
                r[i][i] -= 1;

            cout << matr_inf_norm(r) << ',';
            
            // ||R_gen||
            double r_gen[N][N];
            matr_mul(a, a_inv, r_gen);
            for (int i = 0; i < N; i++)
                r_gen[i][i] -= 1;
            
            cout << matr_inf_norm(r_gen) << endl;

            if (changing_border == 'a')
                alpha /= 10;
            
            else
                beta *= 10;
        }
    }

    else {
        cerr << "Формат команды: gen <α> <β> <a/b> <ζ_max> <y/i/j>" << endl;
        cerr << "<a/b> — какую из границ собственнных значений менять (α / β)" << endl;
        cerr << "<ζ_max> — до какого значения относительной ошибки продолжать генерацию" << endl;
        cerr << "<y/i/j> — матрицу какого типа генерировать (симметричную / простой структуры / с жордановой клеткой)" << endl;
    }

    return 0;
}


void mygen(double a[N][N], double a_inv[N][N], double alpha, double beta, int sign_law, int lambda_law, int type, int schema) {
    int i, j, k;
    
    // распределение знаков
    double sign[N];
    for (i = 0; i < N; i++)
        sign[i] = 1;

    switch (sign_law) {
        case -1:
            for (i = 0; i < N; i++)
                sign[i] = -1;
            break;
        
        case 0:
            sign[0] = 1;
            for (i = 1; i < N; i++)
                sign[i] = -sign[i - 1];
            break;
    }

    // распределение собственнных чисел
    double kappa[N];
    for (i = 0; i < N; i++)
        kappa[i] = (double)i / double(N - 1);
    
    switch (lambda_law) {
        case 1: // κ = sqrt()
            for (i = 0; i < N; i++)
                kappa[i] = sqrt(kappa[i]);
            break;
            
        case 2: // κ = sin()
            double pi_half = .5 * M_PI;
            for (i = 0; i < N; i++)
                kappa[i] = sin(pi_half * kappa[i]);
            break;
    }
    
    double J[N];
    for (i = 0; i < N; i++)
        J[i] = sign[i] * ((1 - kappa[i]) * alpha + kappa[i] * beta);

    double J_inv[N];
    for (i = 0; i < N; i++)
        J_inv[i] = 1. / J[i];

    double Q[N][N];
    double aa[3];
    
    switch (type) {
        case 0: // симметричная матрица
            switch (schema) {
                case 1:
                    Q_matrix(Q, schema);
                    
                    a[0][0] = 0;
                    
                    for (k = 0; k < N; k++)
                        a[0][0] += Q[k][0] * J[k] * Q[k][0];
                    
                    for (j = 1; j < N; j++) {
                        a[j][0] = 0;
                        
                        for (k = j - 1; k < N; k++)
                            a[j][0] += Q[k][0] * J[k] * Q[k][j];
                        
                        a[0][j] = a[j][0];
                    }
                    
                    for (i = 1; i < N; i++) {
                        a[i][i] = 0;
                        
                        for (k = i - 1; k < N; k++)
                            a[i][i] += Q[k][i] * J[k] * Q[k][i];
                        
                        for (j = i + 1; j < N; j++) {
                            a[j][i] = 0;
                            
                            for (k = j - 1; k < N; k++)
                                a[j][i] += Q[k][i] * J[k] * Q[k][j];
                            
                            a[i][j] = a[j][i];
                        }
                    }
                    
                    a_inv[0][0] = 0;
                    
                    for (k = 0; k < N; k++)
                        a_inv[0][0] += Q[k][0] * J_inv[k] * Q[k][0];
                    
                    for (j = 1; j < N; j++) {
                        a_inv[j][0] = 0;
                        
                        for (k = j - 1; k < N; k++)
                            a_inv[j][0] += Q[k][0] * J_inv[k] * Q[k][j];
                        
                        a_inv[0][j] = a_inv[j][0];
                    }
                    
                    for (i = 1; i < N; i++) {
                        a_inv[i][i] = 0;
                        
                        for (k = i - 1; k < N; k++)
                            a_inv[i][i] += Q[k][i] * J_inv[k] * Q[k][i];
                        
                        for (j = i + 1; j < N; j++) {
                            a_inv[j][i] = 0;
                            
                            for (k = j - 1; k < N; k++)
                                a_inv[j][i] += Q[k][i] * J_inv[k] * Q[k][j];
                            
                            a_inv[i][j] = a_inv[j][i];
                        }
                    }
                    
                    break;
            }
            
            break;
        
        case 1: // матрица простой структуры
            switch (schema) {
                case 1:
                    // TJ
                    // первая строка
                    a[0][0] = J[0];
                    a[1][0] = -J[1];
                    
                    // до последней
                    for (i = 1; i < N - 1; i++) {
                        a[i - 1][i] = -J[i - 1];
                        a[i][i] = J[i] + J[i];
                        a[i + 1][i] = -J[i + 1];
                    }
                    
                    // последняя (N - 1)
                    a[N - 2][N - 1] = -J[N - 2];
                    a[N - 1][N - 1] = J[N - 1] + J[N - 1];

                    // (TJ)T^{-1}
                    // первая строка
                    aa[1] = a[0][0];
                    aa[2] = a[1][0];
                    a[0][0] = aa[1] * N + aa[2] * (N - 1);
                    
                    double s = aa[1] + aa[2];
                    for (j = 1; j < N; j++)
                        a[j][0] = s * (N - j);
                    
                    // до последней
                    for (i = 1; i < N - 1; i++) {
                        aa[0] = a[i - 1][i];
                        aa[1] = a[i][i];
                        aa[2] = a[i + 1][i];
                        
                        for (j = 0; j < i; j++)
                            a[j][i] = aa[0] * (N - i + 1) + aa[1] * (N - i) + aa[2] * (N - i - 1);
                        
                        s = aa[0] + aa[1];
                        a[i][i] = s * (N - i) + aa[2] * (N - i - 1);
                        s += aa[2];
                        
                        for (j = i + 1; j < N; j++)
                            a[j][i] = s * (N - j);
                    }
                    
                    // последняя (N - 1)
                    aa[0] = a[N - 2][N - 1];
                    aa[1] = a[N - 1][N - 1];
                    s = aa[0] + aa[0] + aa[1];
                    
                    for (j = 0; j < N - 1; j++)
                        a[j][N - 1] = s;
                    
                    a[N - 1][N - 1] = aa[0] + aa[1];

                    // TJ^{-1}
                    // первая строка
                    a_inv[0][0] = J_inv[0];
                    a_inv[1][0] = -J_inv[1];
            
                    // до последней
                    for (i = 1; i < N - 1; i++) {
                        a_inv[i - 1][i] = -J_inv[i - 1];
                        a_inv[i][i] = J_inv[i] + J_inv[i];
                        a_inv[i + 1][i] = -J_inv[i + 1];
                    }
                    
                    // последняя (N - 1)
                    a_inv[N - 2][N - 1] = -J_inv[N - 2];
                    a_inv[N - 1][N - 1] = J_inv[N - 1] + J_inv[N - 1];

                    // (TJ^{-1})T^{-1}
                    // первая строка
                    aa[1] = a_inv[0][0];
                    aa[2] = a_inv[1][0];
                    a_inv[0][0] = aa[1] * N + aa[2] * (N - 1);
                    s = aa[1] + aa[2];
                    
                    for(j = 1; j < N; j++)
                        a_inv[j][0] = s * (N - j);
            
                    // до последней
                    for (i = 1; i < N - 1; i++) {
                        aa[0] = a_inv[i - 1][i];
                        aa[1] = a_inv[i][i];
                        aa[2] = a_inv[i + 1][i];
                        
                        for (j = 0; j < i; j++)
                            a_inv[j][i] = aa[0] * (N - i + 1) + aa[1] * (N - i) + aa[2] * (N - i - 1);
                        
                        s = aa[0] + aa[1];
                        a_inv[i][i] = s * (N - i) + aa[2] * (N - i - 1);
                        s += aa[2];
                        
                        for (j = i + 1; j < N; j++)
                            a_inv[j][i] = s * (N - j);
                    }
                    
                    // последняя (N - 1)
                    aa[0] = a_inv[N - 2][N - 1];
                    aa[1] = a_inv[N - 1][N - 1];
                    s = aa[0] + aa[0] + aa[1];
                    
                    for (j = 0; j < N - 1; j++)
                        a_inv[j][N - 1] = s;
                    
                    a_inv[N - 1][N - 1] = aa[0] + aa[1];
                    
                    break;
            }
            
            break;

        case 2: // одна жорданова клетка 2✕2 при минимальном с. з.
            switch (schema) {
                case 1:
                    // TJ
                    // первая строка
                    a[0][0] = J[0];
                    a[1][0] = 1 - J[0];
            
                    // вторая строка
                    a[0][1] = -J[0];
                    a[1][1] = -1 + J[0] + J[0];
                    a[2][1] = -J[2];
            
                    // третья строка
                    a[1][2] = -J[0];
                    a[2][2] = J[2] + J[2];
                    
                    if (N > 3)
                        a[3][2] = -J[3];
                    
                    // до последней
                    for (i = 3; i < N - 1; i++) {
                        a[i - 1][i] = -J[i - 1];
                        a[i][i] = J[i] + J[i];
                        a[i + 1][i] = -J[i + 1];
                    }
                    
                    // последняя (N - 1)
                    if (N > 3) {
                        a[N - 2][N - 1] = -J[N - 2];
                        a[N - 1][N - 1] = J[N - 1] + J[N - 1];
                    }
                    
                    // (TJ)T^{-1}
                    // первая строка
                    aa[1] = a[0][0];
                    aa[2] = a[1][0];
                    a[0][0] = aa[1] * N + aa[2] * (N - 1);
            
                    double s = aa[1] + aa[2];
                    for (j = 1; j < N; j++)
                        a[j][0] = s * (N - j);
            
                    // до последней
                    for (i = 1; i < N - 1; i++) {
                        aa[0] = a[i - 1][i];
                        aa[1] = a[i][i];
                        aa[2] = a[i + 1][i];
                        
                        for (j = 0; j < i; j++)
                            a[j][i] = aa[0] * (N - i + 1) + aa[1] * (N - i) + aa[2] * (N - i - 1);
                        
                        s = aa[0] + aa[1];
                        a[i][i] = s * (N - i) + aa[2] * (N - i - 1);
                        s += aa[2];
                        
                        for (j = i + 1; j < N; j++)
                            a[j][i] = s * (N - j);
                    }
                    
                    // последняя (N - 1)
                    aa[0] = a[N - 2][N - 1];
                    aa[1] = a[N - 1][N - 1];
                    s = aa[0] + aa[0] + aa[1];
            
                    for (j = 0; j < N - 1; j++)
                        a[j][N - 1] = s;
                    
                    a[N - 1][N - 1] = aa[0] + aa[1];

                    // TJ^{-1}
                    // первая строка
                    a_inv[0][0] = J_inv[0];
                    a_inv[1][0] = -J_inv[0] * J_inv[0] - J_inv[0];
                    
                    // вторая строка
                    a_inv[0][1] = -J_inv[0];
                    a_inv[1][1] = J_inv[0] * J_inv[0] + J_inv[0] + J_inv[0];
                    a_inv[2][1] = -J_inv[2];
                    
                    // третья строка
                    a_inv[1][2] = -J_inv[0];
                    a_inv[2][2] = J_inv[2] + J_inv[2];
                    
                    if (N > 3)
                        a_inv[3][2] = -J_inv[3];
                    
                    // до последней
                    for (i = 3; i < N - 1; i++) {
                        a_inv[i - 1][i] = -J_inv[i - 1];
                        a_inv[i][i] = J_inv[i] + J_inv[i];
                        a_inv[i + 1][i] = -J_inv[i + 1];
                    }
                    
                    // последняя (N - 1)
                    if (N > 3) {
                        a_inv[N - 2][N - 1] = -J_inv[N - 2];
                        a_inv[N - 1][N - 1] = J_inv[N - 1] + J_inv[N - 1];
                    }
                    
                    // (TJ^{-1})T^{-1}
                    // первая строка
                    aa[1] = a_inv[0][0];
                    aa[2] = a_inv[1][0];
                    a_inv[0][0] = aa[1] * N + aa[2] * (N - 1);
                    
                    s = aa[1] + aa[2];
                    for (j = 1; j < N; j++)
                        a_inv[j][0] = s * (N - j);
                    
                    // до последней
                    for (i = 1; i < N - 1; i++) {
                        aa[0] = a_inv[i - 1][i];
                        aa[1] = a_inv[i][i];
                        aa[2] = a_inv[i + 1][i];
                        
                        for (j = 0; j < i; j++)
                            a_inv[j][i] = aa[0] * (N - i + 1) + aa[1] * (N - i) + aa[2] * (N - i - 1);
                        
                        s = aa[0] + aa[1];
                        a_inv[i][i] = s * (N - i) + aa[2] * (N - i - 1);
                        s += aa[2];
                        
                        for (j = i + 1; j < N; j++)
                            a_inv[j][i] = s * (N - j);
                    }
                    
                    // последняя (N - 1)
                    aa[0] = a_inv[N - 2][N - 1];
                    aa[1] = a_inv[N - 1][N - 1];
                    
                    s = aa[0] + aa[0] + aa[1];
                    for (j = 0; j < N - 1; j++)
                        a_inv[j][N - 1] = s;
                    
                    a_inv[N - 1][N - 1] = aa[0] + aa[1];
                    
                    break;
            }
            
            break;
    }
}


void Q_matrix(double Q[N][N], int schema) {
    int i, j;
    double q, curr, next = 1;
    
    for (j = 0; j < N - 1; j++) {
        curr = next;
        next += 1;
        
        q = 1. / sqrt(curr * next);
        
        for (i = 0; i <= j; i++)
            Q[j][i] = q;
        
        Q[j][j + 1] = -sqrt(curr / next);
        
        for (i = j + 2; i < N; i++)
            Q[j][i] = 0;
    }
    
    q = 1. / sqrt(N);
    
    for (i = 0; i < N; i++)
        Q[N - 1][i] = q;
}


void matr_mul(double a[N][N], double b[N][N], double c[N][N]) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            c[j][i] = 0;
            
            for (int k = 0; k < N; k++)
                c[j][i] += a[k][i] * b[j][k];
    }
}


double matr_inf_norm(double a[N][N]) {
    double s, norm = 0;
    
    for (int i = 0; i < N; i++) {
        s = 0;
        
        for (int j = 0; j < N; j++)
            s += fabs(a[j][i]);
        
        if (s > norm)
            norm = s;
    }
    
    return norm;
}
