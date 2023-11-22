#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <windows.h>


using namespace std;

const string path = "d: / port / OUT / ";
const double f = 0.0001;
const double rb = 0.0001;
const double dt = 0.01;
const int nx = 385;
const int ny = 420;
const double g = 9.8;
const double dx = 2;
const int tn = 100000;
const int kx = 1000;
const int ro = 1000;

typedef double tt1[nx][ny];
int j0_0; //j0 существует в c++
double ww, hs, k0, rs, mt, bb, c, l, k, T0, w00, w0, h01, hh, w1, r;
int nx1, i, j, t;
string ts;
int Xg[ny];
tt1 Lm0, Lm, Ks, Cs, h, h0, u, v, fx, fy; //сет
double gr1[nx][101], gr2[nx][101], gr3[nx][101];

//гиперболический тангенс
double th(double a1) {
	return (exp(a1) - exp(-a1)) / (exp(a1) + exp(-a1));
}
//распространением волн в закрытой акватории
double KW(double TT, double H0) {
    double KK1, a, b;
    double KK = 1;

    do {
        b = KK * H0;
        KK1 = 4 * M_PI * M_PI / (TT * TT * g * th(b));
        a = abs((KK - KK1) / KK);
        KK = KK1;
    } while (a > 0.01);

    return KK;
}

double plan(tt1 zz) {
    int i, j;
    double d, hmax, hmin;
    do {
        hmax = 0;
        hmin = 0;
        for (i = 2; i <= (nx - 1); i++) {
            for (j = 2; j <= (ny - 1); j++) {
                if (zz[i][j] > hmax) {
                    hmax = zz[i][j];
                }
                if (zz[i][j] < hmin) {
                    hmin = zz[i][j];
                }
            }
        }
        d = round(10 * (hmax - hmin));
    } while (d < 0.0000001);
    return d;
}

void MAX(int t, tt1 Lm0, tt1 Lm, tt1 h0){
    int i,j;
    double tt;

    tt = t % 10000;
    if ((abs(tt) < 0.001) || (t == 0)) {
        for (i = 2; i <= nx - 1; i++) {
            for (j = 2; j <= ny - 1; j++) {
                Lm0[i][j] = Lm[i][j];
                Lm[i][j] = 0;
            }
        }
    }
    for (i = 2; i <= nx - 1; i++) {
        for (j = 2; j <= ny - 1; j++) {
            if (Lm[i][j] < abs(h0[i][j])) {
                Lm[i][j] = abs(h0[i][j]);
            }
        }
    }
    cout << Lm[0][0] << endl;
}

int main() {
    double a1 = 1.2;
    double th_res = th(a1);
    cout << th_res << endl;
    double TT = 10.0;  // Замените значение на необходимое
    double H0 = 1.0;   // Замените значение на необходимое
    double KW_res = KW(TT, H0);
    cout << KW_res << endl;
    // генерация случайных значений h0(Позже из БД)
    for (int i = 0; i <= nx; i++)
    {
        for (int j = 0; j <= ny; j++)
        {
            //sw
            h0[i][j] = rand();
        }
    }
    double plan_res = plan(h0);
    //cout << plan_res << endl;
    // генерация случайных значений Lm(Позже из БД)
    for (int i = 0; i <= nx; i++)
    {
        for (int j = 0; j <= ny; j++)
        {
            //SW
            Lm[i][j] = rand();
        }
    }
    int t1 = 5;
    MAX(t1, Lm0, Lm, h0);
    //период
    cin >> T0;
    //амплитуда колебания
    cin >> w00;
    //глубина
    cin >> hh;
    //волновое число
    cin >> k;
    //фазовая скорость
    c = 2 * M_PI / T0 / k;
    cout << c;
    //отражение правой границы
    cin >> rs;
    //нижняя граница мола
    cin >> j0_0;
    //cout << hh << endl;
    r = rb / hh;
    for (int j = 1; j <= ny; j++) {
        if (j < 280) {
            Xg[j] = 105 + (j - 1);
        }
        else {
            Xg[j] = 385;
        }
    }
    //cout << Xg[3] <<endl;
    for (int i = 1; i <= nx; i++) {
        for (int j = 1; j <= ny; j++) {
            h0[i][j] = 0;
            Lm[i][j] = 0;
            //cout << h[i][j] << endl;
            //cout << Xg[j] << endl;
            //cout << i << endl;
            if (i <= Xg[j]) {
                //cout << "HH" << hh << endl;
                h[i][j] = hh * (1 - pow((j - 210), 2) / (pow(210, 2)) * (1 - pow(i, 2) / pow(Xg[j], 2))) + 1;
                //cout << h[i][j];
            }
            else {
                h[i][j] = 0;
            }
            if ((j == 1) || (j == ny) || (i == nx)) {
                h[i][j] = 0;
            }
            if ((nx - i) + j < 245) {
                h[i][j] = 0;
            }
            //cout << h[i][j] << endl;
            u[i][j] = 0;
            v[i][j] = 0;
            hs = h[i][j];
            if (hs > 0) {
                Ks[i][j] = KW(T0,hs);
                Cs[i][j] = 2 * M_PI / (T0 * Ks[i][j]);
            }
            else {
                Cs[i][j] = 0;
            }
            //cout << Ks[i][j] << endl;
            //cout << Cs[i][j] << endl;
        }
    }
    for (int t = 0; t < tn; t++) {
        mt = t % 100;
        if (mt < 0.01) {
            MAX(t1, Lm0, Lm, h0);
        }
        for (int i = 2; i <= nx - 1; i++) {
            for (int j = 2; j <= ny - 1; j++) {
                //cout << h0[i][j];
                //ЗАМЕНИТЬ h0 на h
                fx[i][j] = -g * (h[i][j] - h[i - 1][j]) / dx;
                fy[i][j] = -g * (h[i][j] - h[i][j - 1]) / dx;
                u[i][j] = u[i][j] * exp(-r * dt) + fx[i][j] / r * (1 - exp(-r * dt));
                v[i][j] = v[i][j] * exp(-r * dt) + fy[i][j] / r * (1 - exp(-r * dt));
                //cout << fx[i][j];
                //cout << fy[i][j];
                //cout << u[i][j];
                //cout << v[i][j];
            }
        }
        for (int j = 1; j <= ny; j++) {
            //cout << -g/ c;
            u[1][j] = -g / c * h[1][j];
            //cout << u[1][j];
            u[nx][j] = 0;
            //cout << u[1][j];
        }
        for (int i = 1; i <= nx; i++) {
            v[i][1] = 0;
            v[i][ny] = 0;
            //cout << v[i][1];
        }
        for (int j = j0_0; j <= ny - 1; j++) {
            //cout << u[20][j];
            u[20][j] = 0;
            u[21][j] = 0;
            //cout << u[20][j];
        }
        for (int i = 1; i <= nx - 1; i++) {
            //cout << '1';
            for (int j = 1; j <= ny - 1; j++) {
                //cout << '2';
                if (h[i][j] > 0) {
                    //cout << '3';
                    w0 = 0;
                    //cout << bb;
                    bb = Ks[i][j] * h[i][j];
                    //cout << bb;
                    if (i == 3) {
                        //отклонение колеблющейся величины в текущий момент времени t от среднего за период значения
                        w0 = w00 * sin(2 * M_PI / T0 * (t - 1) * dt);
                    }
                    k0 = Ks[i][j];
                    ww = w0 - (u[i + 1][j] - u[i][j]) / (pow(dx * th(bb), 2) * (k0 - (v[i][j + 1] - v[i][j])) * k0);
                    h[i][j] = h[i][j] + ww * dt;
                    //cout << h[i][j];
                }
            }
        }
        cout << h[5][5];
        //Градиент - от 0 до 99, х - 100, у - 101
        for (int i = 1; i <= nx; i++) {
            //cout << i << endl;
            gr1[i][100] = 20 + 2 * (i - 1);
            //cout << gr1[i][100] << endl;
            gr2[i][100] = 20 + 2 * (i - 1);
            gr3[i][100] = 20 + 2 * (i - 1);
            gr1[i][101] = 200 - round(30 * h[i][25]);
            gr2[i][101] = 200 + 1 - round(30 * h[i][50]);
            gr3[i][101] = 200 + 2 - round(30 * h[i][75]);
            //cout << gr3[i][101];
        }
        //cout << "end_1" << endl;
    }
    cout << "end" << endl;

}