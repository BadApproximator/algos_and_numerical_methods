#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#define eps 0.00001

using namespace std;

ofstream out("out.txt");

double f(double x)
{
    return exp(2*x*x);
    //return sin(x);
}

vector<double> shooting(double c, double N)
{
    vector<double> y(N + 1, 0);
    double h = 1. / N;
    y[N] = 1;
    y[N - 1] = c; //shooting method
    //cout << endl << "c = " << c << endl;
    for (int i = N - 1; i >= 1; --i) {
        y[i - 1] = -y[i + 1] + 2 * y[i] - (f(i * h) - i * h * pow(y[i], 3)) * h * h;
        cout << "i = " << i - 1 << " " << y[i - 1] << endl;
    }
    return y;
}

double calculate_integral(double c, double N)
{
    double res = 0;
    vector<double> y(N + 1, 0);
    y = shooting(c, N);
    double h = 1. / N;
    for (int i = 0; i <= N - 1; ++i) {
        res += sin(i * h) * y[i] + sin(i * h + h) * y[i + 1];
    }
    res = res * h / 2;
    return res;
}

double root(double a, double b, double N)
{
    double t, f1, f2, x;
    do {
        f1 = calculate_integral(a, N) - 1;
        t = (a + b) / 2.0;
        f2 = calculate_integral(t, N) - 1;
        if (f1 * f2 <= 0) b = t;
        else a = t;
    } while (fabs(b - a) > eps);
    x = (a + b) / 2.0;
    f1 = calculate_integral(x, N) - 1;
    if (fabs(f1) <= eps) {
        cout << "\nRoot with error ";
        cout << fixed << eps;
        cout << ", X = ";
        cout << x;
        cout << "\nFunction value F(X) = " << f1 << endl;
    } else cout << "No roots!" << endl;

    return x;
}

int main()
{
    double N = 20;
    double c; // for shooting method
    double h;
    //cin >> N;
    h = 1. / N;
    vector<double> x;
    vector<double> y(N + 1);
    for (int i = 0; i < N + 1; ++i) {
        x.push_back(i * h);
    }
    y[N] = 1;/*
    y = shooting(1./2, N);
    for (int j = 0; j <= N; ++j) {
        cout << y[j] << " ";
    }
    cout << endl;*/
    //cout << "integral = " << calculate_integral(, N);
    c = root(0.9, 1.4, N);
    /*cout << "const in u[1-ih] is " << c <<endl;
    y = shooting(c, N);
    for (int i = 0; i < N+1; ++i) {
        cout << y[i] << endl;
    }*/
    y = shooting(c, N);
    cout << "All good, result of shooting = " << c;
    for (int i = 0; i < N + 1; ++i) {
        out << x[i] << " " << y[i] << endl;
    }
    return 0;
}