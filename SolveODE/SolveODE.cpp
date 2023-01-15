#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

#define eps 1e-9
#define PI 3.14159265358979323846

using namespace std;

vector<double> operator*(const vector<double>& v, double alpha)
{
    vector<double> result;
    for (int i = 0; i < v.size(); ++i) {
        result.push_back(alpha * v[i]);
    }
    return result;
}

vector<double> operator+(const vector<double>& v1, const vector<double>& v2)
{
    vector<double> sum;
    for (int i = 0; i < v1.size(); ++i) {
        sum.push_back(v1[i] + v2[i]);
    }
    return sum;
}

vector<double> operator-(const vector<double>& v1, const vector<double>& v2)
{
    vector<double> sum;
    for (int i = 0; i < v1.size(); ++i) {
        sum.push_back(v1[i] - v2[i]);
    }
    return sum;
}

double operator*(const vector<double>& v1, const vector<double>& v2)
{
    double sum = 0;
    if (v1.size() != v2.size()) cout << "Error! Vectors have different dimentions!" << endl;
    for (int i = 0; i < v1.size(); ++i) {
        sum += v1[i] * v2[i];
    }
    return sum;
}

vector <double> multiple_mat_vec(const vector < vector <double>>& a, vector < double > x) {
    int N = a[0].size();
    vector <double> y(N, 0);
    // i - number of coord, j - number of column
    for (int i = 0; i < N; ++i) {
        double curr_coord = 0;
        for (int j = 0; j < N; j++) {
            curr_coord += (a[i][j] * x[j]);
        }
        y[i] = curr_coord;
    }
    return y;
}

void output_solution(const vector <double>& x, double L) {
    string filename;
    filename = "solution.out";
    ofstream out(filename);
    int N = x.size();
    double step = L / N;
    int k = 0;
    for (double i = 0, pos = 0; i < N; ++i, pos+=step) {
        out << pos << " " << x[i] << endl;
    }
}

vector <double> f(const vector < vector <double>>& M, const vector <double> y) {
    return multiple_mat_vec(M, y);
}

vector <double> find_y_star(const vector < vector <double>>& M, const vector <double> y, const double h) {
    vector <double> y_star;
    y_star = y + f(M, y) * h;
    return y_star;
}

vector <double> find_next_y(const vector < vector <double>>& M, const vector <double> y, double h) {
    vector <double> y_next;
    vector <double> y_star = find_y_star(M, y, h);
    y_next = y + (f(M, y + y_star)) * (h / 2); // trapezoid method O(h^2)
    //y_next = y + f(M, y) * h; // Euler's method O(h)
    return y_next;
}

vector <double> find_solution(const vector < vector <double>>& M, const vector <double> y0, int N, double L) {
    const double h = L / N;
    vector <double> solution;
    vector <double> y = find_next_y(M, y0, h);
    solution.push_back(y0[0]);
    solution.push_back(y[0]);
    for (int i = 1; i < N; ++i) {
        y = find_next_y(M, y, h);
        solution.push_back(y[0]);
    }
    return solution;
}

int main() {
    vector < vector <double>> M = {
        {0,1,0,0,0},
        {0,0,1,0,0},
        {0,0,0,1,0},
        {0,0,0,0,1},
        {-243,-405,-270,-90,-15} };
    vector < double > y0 = {0,3,-9,-8,0};
    unsigned const int N = 1000;
    const double L = 5.;
    double h = L / N;
    cout << "Grid N, where N = " << N << endl;
    cout << "Step h = " << h << endl;

    vector <double> y = find_solution(M, y0, N, L);

    //print_vector(y);
    output_solution(y, L);
    return 0;
}