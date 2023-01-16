#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <time.h>
#include <string>

#define eps 1e-6

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
    if(v1.size() != v2.size()) cout << "Error! Vectors have different dimentions!" << endl;
    for (int i = 0; i < v1.size(); ++i) {
        sum += v1[i] * v2[i];
    }
    return sum;
}

double squarenorm(const vector<double>& v)
{
    double sum = 0;
    for (int i = 0; i < v.size(); ++i) {
        sum += v[i]*v[i];
    }
    return sum;
}

double func(double x, double y){
    return 1./(x+1) * y;
    //return x*x*y*25;
    //return sin(pow(x,3)*y)*exp(10*y*y);
    //return pow(sin((1./5-x)*y),3)*10;
    //return x*y*1000;
}

void print_matrix(const vector < vector <double>> &A){
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].size(); ++j) {
            printf("%3.0lf ", A[i][j]);
        }
        cout << endl;
    }
}

void print_vector(const vector <double> &v){
    int M = v.size();
    int N = sqrt(M) + 1;
    /*cout<< "(";
    for (int i = 0; i < v.size(); ++i) {
        printf("%8.2lf", v[i]);
    }
    cout<< ")" << endl;*/
    int l = 0;
    cout << endl;
    for (int i = 0; i < N - 1; ++i) {
        for (int j = 0; j < N - 1; ++j) {
            printf("%8.2lf", v[l++]);
        }
        cout << endl;
    }
}

void create_matrix(vector < vector <double>> &A, const int &N){
    int M = (N-1)*(N-1);
    vector< vector<double> > a(M, vector<double>(M, 0)); // два размера
    int d, dup, ddown; //d - diagonal
    for (int i = 0; i < N - 1; ++i) { //i - number of block
        for (int j = 0; j < N - 1; ++j) { //j - number of row of the block
            d = (N-1)*i+j;
            dup = d + N-1;
            ddown = d - (N-1);
            //first row in the block
            if(j == 0){
                a[d][d] = 4;
                a[d][d + 1] = -1;
            }
                //last row in the block
            else if(j == N-2){
                a[d][d] = 4;
                a[d][d - 1] = -1;
            }
            else{
                a[d][d] = 4;
                a[d][d + 1] = -1;
                a[d][d - 1] = -1;
            }
            /* filling pseudodiagonals by -1 */
            if(i == 0){
                a[d][dup] = -1;
            }
            else if(i == N-2){
                a[d][ddown] = -1;
            }
            else{
                a[d][dup] = -1;
                a[d][ddown] = -1;
            }
        }
    }
    A = a;
    a.clear();
}

void create_rightpart(vector <double> &b, const int &N, double alpha){
    int M = (N-1)*(N-1);
    vector <double> v (M, 0);
    for (int i = 0; i < N-1; ++i) {
        for (int j = 0; j < N - 1; ++j) {
            v[i*(N-1) + j] = func((double)i/N, (double)j/N);
        }

    }
    double h = 1./N;
    v = v*h*h;
    for (int i = 0; i < N - 1; ++i) {
        for (int j = 0; j < N - 1; ++j) {
            if (i==0) v[j] += alpha;
            if (j==0) v[i*(N-1) + j] += alpha;
            if (i==N-2) v[i*(N-1) + j] += alpha;
            if (j==N-2) v[i*(N-1) + j] += alpha;
        }
    }
    b = v;
    v.clear();
}

vector <double> mutiple_mat_vec(const vector < vector <double>> &a, vector < double > x){
    int N = sqrt(a[0].size()) + 1;
    //cout << "N = "<< N << endl;
    int M = (N-1)*(N-1);
    double h = 1./N;
    int d = 0, dup, ddown; //d - diagonal
    vector <double> y (M, 0);
    // i - number of block, j - number of row of the block
    for (int i = 0; i < N - 1; ++i) {
        for (int j = 0; j < N - 1; j++, d++) {
            // d - current number of row of all matrix a
            dup = d + N-1;
            ddown = d - (N-1);

            //first row in the block
            if(j == 0){
                y[d] += a[d][d] * x[d];
                y[d] += a[d][d + 1] * x[d+1];
            }
                //last row in the block
            else if(j == N-2){
                y[d] += a[d][d] * x[d];
                y[d] += a[d][d - 1] * x[d-1];
            }
            else{
                y[d] += a[d][d] * x[d];
                y[d] += a[d][d+1] * x[d+1];
                y[d] += a[d][d-1] * x[d-1];
            }
            /* filling pseudodiagonals by -1 */
            if(i == 0){
                y[d] += a[d][dup] * x[dup];
            }
            else if(i == N-2){
                y[d] += a[d][ddown] * x[ddown];
            }
            else{
                y[d] += a[d][dup] * x[dup];
                y[d] += a[d][ddown] * x[ddown];
            }
            //cout << y[d] << endl;
        }
    }
    return y;
}

vector <double> gradient_descent(const vector <vector <double>> &A, vector <double> x, const vector <double> &f){
    double value_of_step; // норма разности текущего и предыдущего приближения
    int k = 0; // номер итерации
    vector <double> xi; //вектор кси = вектор невязки
    double alpha; // коэффициент в методе наискорейшего градиентного спуска
    vector <double> A_xi, x_prev; // вспомогательные векторы
    xi = mutiple_mat_vec(A, x) - f;
    //print_vector(xi);
    do{
        x_prev = x;
        A_xi = mutiple_mat_vec(A, xi); // умножение матрицы A на вектор кси
        alpha = squarenorm(xi)/(A_xi*xi); // просто формула из лекций
        x = x - (xi*alpha); // тоже просто формула, задающая итерационный процесс
        xi = xi - (A_xi*alpha); // нужно вычислять, чтобы было только одно умножение матрицы на вектор

        if (squarenorm(xi) < eps){ // проверяем, что вектор невязки ненулевой
            cout << "vector xi is null" << endl;
            break;
        }
        value_of_step = squarenorm(x - x_prev); // вычисляем насколько изменилось приближение
        //cout << "squarenorm(x - x_prev) = " << value_of_step << endl;
        k++; // увеличиваем номер итерации
    } while (value_of_step > eps);
    cout << "Number of interations in [Gradient Descent Method] = " << k << endl;
    return x;
}

vector <double> seidel_method(const vector <vector <double>> &A, vector <double> x, const vector <double> &f){
    double value_of_step; // норма разности текущего и предыдущего приближения
    int k = 0; // номер итерации
    int M = x.size();
    int N = sqrt(M) + 1;
    vector <double> x_next (M, 0);
    do {
        //cout << "Iteration [" << k << "]:" << endl;
        for (int j = 0; j < M; ++j) {
            /* Это медленное вычисление
             * for (int i = 0; i < j; ++i) {
                x_next[j] -= A[j][i]*x_next[i];
            }
            for (int i = j+1; i < M; ++i) {
                x_next[j] -= A[j][i]*x[i];
            }*/
            if (j <= N-2){
                if (j == 0) x_next[j] = -A[j][j+1]*x[j+1] - A[j][j+N-1] * x[j+N-1];
                else x_next[j] = -A[j][j-1]*x_next[j-1] - A[j][j+1]*x[j+1] - A[j][j+N-1] * x[j+N-1];
            }
            else if (j >= N-1 && j < M-(N-1)){
                x_next[j] = -A[j][j-(N-1)]*x_next[j-(N-1)]-A[j][j-1]*x_next[j-1] - A[j][j+1]*x[j+1] - A[j][j+N-1] * x[j+N-1];
            }
            else if (j < M){
                if (j == M-1) x_next[j] = -A[j][j-(N-1)]*x_next[j-(N-1)]-A[j][j-1]*x_next[j-1];
                else x_next[j] = -A[j][j-(N-1)]*x_next[j-(N-1)]-A[j][j-1]*x_next[j-1] - A[j][j+1]*x[j+1];
            }
            x_next[j] = (x_next[j]+f[j])/A[j][j];
        }
        value_of_step = squarenorm(x_next - x);
        //cout << "squarenorm(x_next - x) = " << value_of_step << endl;
        //print_vector(x_next);
        x = x_next;
        k++;
    }
    while (value_of_step > eps);
    cout << "Number of interations in [Seidel Method] = " << k << endl;
    return x;
}

void output_solution(const vector <double> &x, double alpha, int flag){
    string filename;
    if(flag == 0) filename = "gradient.out";
    if(flag == 1) filename = "seidel.out";
    ofstream out(filename);
    int M = x.size();
    double N = sqrt(M) + 1;
    int k = 0;
    for (int i = 0; i < N-1; ++i) {
        for (int j = 0; j < N - 1; ++j) {
            //if(i == 0 || j == 0 || i == N) out << i/N << " " << j/N << " " << alpha << endl;
            out << (i+1)/N << " " << (j+1)/N << " " << x[k] << endl;
            k++;
        }
    }
    for (int i = 0; i <= N; i++){
        out << '0' << " " << i/N << " " << alpha << endl;
    }
    for (int i = 1; i <= N; i++){
        out << i/N << " " << '0' << " " << alpha << endl;
    }
    for (int i = 1; i <= N; i++){
        out << '1' << " " << i/N << " " << alpha << endl;
    }
    for (int i = 1; i <= N-1; i++){
        out << i/N << " " << '1' << " " << alpha << endl;
    }

}

vector <double> residual (const vector < vector <double>> &a, const vector < double > &x, const vector < double > &b){
    int N = x.size();
    vector <double > res;
    res = mutiple_mat_vec(a, x);
    res = res - b;
    return res;
}

int main() {
    vector < vector <double>> A;
    vector < double > b;
    double alpha = 1;
    unsigned int N = 20;
    //cin >> N;
    int M = (N-1)*(N-1);
    vector < double > x (M, 1/3);
    vector <double> grad (M), seidel(M);
    double resgrad (M), resseidel(M); //residuals
    double h;
    h = 1./N;
    cout << "Grid N x N, where N = " << N << endl;
    cout << "Step h = "<< h << endl;
    create_matrix(A, N);
    create_rightpart(b, N, alpha);
    int start;
    int end;
    int t;
    cout << "--------------" << endl << "GRADIENT DESCENT METHOD..." << endl;
    start = clock();
    grad = gradient_descent(A, x, b);
    end = clock(); // засекаем время окончания
    t = (end - start)/CLOCKS_PER_SEC;// команда CLOCKS_PER_SEC нужна для перевода результата функции clock в секунды
    cout<<"Time: "<< t << " sec" << endl;
    resgrad = sqrt(squarenorm(residual(A, grad, b)));
    cout << "Residual = " << resgrad << endl;
    cout << endl << "--------------" << endl << "SEIDEL METHOD..." << endl;
    start = clock();
    seidel = seidel_method(A, x, b);
    resseidel = sqrt(squarenorm(residual(A, seidel, b)));
    cout << "Residual = " << resseidel << endl;
    end = clock();
    t = (end - start)/CLOCKS_PER_SEC;// команда CLOCKS_PER_SEC нужна для перевода результата функции clock в секунды
    cout<<"Time: "<<t << " sec" << endl;
    cout << "--------------" << endl;
    cout << "Norm(seidel - GDM) = " << squarenorm(seidel - grad) << endl;
    output_solution(grad, alpha, 0);
    output_solution(seidel, alpha, 1);
    return 0;
}