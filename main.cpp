#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

std::vector<int> make_vector(int, int, int);
std::vector<double> make_vector_f(std::vector<double> &);
std::vector<double> make_vector_f_tild(std::vector<double> &, double);
std::vector<double> make_vector_f_result(std::vector<double> &, std::vector<double> &, int, int);
std::vector<double> make_vector_x(double, double, std::vector<int> &, int);
std::vector<double> make_vector_alpha(int, int);
std::vector<double> make_vector_lyambda(std::vector<int> &, int);

void explore(double, double, double, double, double);
double func_expr(double);
double harmonic_mean(std::vector<double> &, std::vector<double> &, int, int);
double my_random(double, double);
void print(std::ofstream&, const std::vector<double>&, const std::vector<double>&);

double e = 0.01;
double x_min = 0;
double x_max = 3.1415926535897932384626433832795;
double q = e / (x_max - x_min);
double P = 0.95;
double a = 0.25; // амплитуда 2a = 0.5, a = +- 0.25

int K = 100;
int L = 10;
int r = 5; // 3;
int M = (r - 1) / 2;
int min_lyambda;

std::vector<int> vector_k;
std::vector<double> vector_lyambda;

std::vector<double> vector_x;
std::vector<double> vector_f;
std::vector<double> vector_f_tild;

std::vector<double> vector_f_result; // заполняются после explore
std::vector<double> vector_alpha; // заполняются после explore

int main()
{
    std::srand(time(0));
    std::ofstream outstream("output_last.txt");
    // 1) получение x
    vector_k = make_vector(0, K, 1);
    vector_x = make_vector_x(x_min, x_max, vector_k, K);
    // 2) получение f
    vector_f = make_vector_f(vector_x);
    // 3) зашумлённой f
    vector_f_tild = make_vector_f_tild(vector_x, a);
    // 4) проведение испытаний и получение новой f
    explore(P, e, x_min, x_max, L);
    // 5) вывод
    print(outstream, vector_x, vector_f);
    print(outstream, vector_x, vector_f_tild);
    print(outstream, vector_x, vector_f_result);
    print(outstream, vector_alpha, vector_alpha);
    outstream.close();
    std::cin.get();
    return 0;
}

double func_expr(double x) {
    return std::sin(x) + 0.5;
}

double my_random(double a, double b) {
    return (rand() % (int((b - a) * 10000) + 1) + a * 10000) / 10000.f;
}

double harmonic_mean(std::vector<double> & f_tild, std::vector<double> & alpha, int k, int M) {
    double result = 0;
    for (int j = k - M; j <= k + M; ++j) {
        if (j - 1 < 0 || j - 1 > 100)
            continue;
        result += alpha.at(j + M + 1 - k - 1) / f_tild.at(j - 1);
    }
    return std::pow(result, -1);
}

std::vector<double> make_vector_f_result(std::vector<double> & vector_f, std::vector<double> & vector_alpha, int K, int M) {
    std::vector<double> result;
    for (int k = 0; k <= K; ++k) {
        result.push_back(harmonic_mean(vector_f, vector_alpha, k, M));
    }
    return result;
}

std::vector<double> make_vector_f_tild(std::vector<double> & vector_x, double amplitude) {
    std::vector<double> result;
    for (auto i : vector_x) {
        auto r = my_random(-amplitude, amplitude);
        result.push_back(func_expr(i) + r);
    }
    return result;
}

std::vector<double> make_vector_f(std::vector<double> & vector_x) {
    std::vector<double> result;
    for (auto i : vector_x) {
        result.push_back(func_expr(i));
    }
    return result;
}

std::vector<double> make_vector_x(double x_min, double x_max, std::vector<int> & vector_k, int K) {
    std::vector<double> result;
    for (auto i : vector_k) {
        result.push_back(x_min + i * (x_max - x_min) / K);
    }
    return result;
}

std::vector<double> make_vector_alpha(int r, int M) {
    M -= 1;
    r -= 1;
    std::vector<double> alpha(r + 1);
    alpha.at(M + 1) = my_random(0, 1);
    alpha.at(M) = alpha.at(M + 2) = 0.5 * my_random(0, 1 - alpha.at(M + 1));

    if (r < 4)
        return alpha;

    for (int m = 2; m < M + 1; ++m) {
        double sum = 0;
        for (int s = m + 1; s < r - m + 1; ++s) {
            sum += alpha.at(s);
        }
        alpha.at(m) = alpha.at(r - m + 1) = 0.5 * my_random(0, 1 - sum);
    }

    double sum = 0;
    for (int s = 2; s < r; ++s) {
        sum += alpha.at(s);
    }
    alpha.at(1) = alpha.at(r) = 0.5 * (1 - sum);
    return alpha;
}

std::vector<int> make_vector(int a, int b, int k) {
    std::vector<int> result;
    for (int i = a; i <= b; i += k) {
        result.push_back(i);
    }
    return result;
}

////
double dist(double w, double d) {
    return abs(w) + abs(d);
}

double get_w(std::vector<double> vector_f, int K) {
    double result = 0;
    for (int k = 1; k < K + 1; ++k) {
        result += abs(vector_f.at(k) - vector_f.at(k - 1));
    }
    return result;
}

double get_d(std::vector<double> vector_f, std::vector<double> vector_f_tild, int K) {
    double result = 0;
    for (int k = 0; k < K + 1; ++k) {
        result += abs(vector_f.at(k) - vector_f_tild.at(k));
    }
    return result / K;
}

int get_N(double P, double e, double x_min, double x_max) {
    return std::round(std::log(1 - P) / std::log(1 - (e / (x_max - x_min))));
}

double get_J(double & lyambda, double w, double d) {
    return (lyambda * w + (1 - lyambda) * d * 10000.f) / 10000.f; //round()
}

void explore(double P, double e, double a, double b, double L) {
    // определение количества испытаний
    auto N = get_N(P, e, a, b);
    // создание вектора лямбда
    vector_lyambda = make_vector_lyambda(make_vector(0, L, 1), L);
    // создание вектора альфа
    std::vector<std::vector<double>> alphas;
    int min_elem = 0; // в векторе лямбда
    int min_alphas = 0; // среди альфа

    for (int i = 0; i < N; ++i) {
        //random vector_alpha
        alphas.push_back(make_vector_alpha(r, M));
        // min
        auto vector_f_min = make_vector_f_result(vector_f_tild, alphas.at(min_alphas), K, M);
        double w_min = get_w(vector_f_min, K);
        double d_min = get_d(vector_f_min, vector_f_tild, K);
        double min_result = get_J(vector_lyambda.at(min_elem), w_min, d_min);
        //
        auto vector_s_result = make_vector_f_result(vector_f_tild, alphas.at(alphas.size() - 1), K, M);
        double w = get_w(vector_s_result, K);
        double d = get_d(vector_s_result, vector_f_tild, K);

        for (auto j = 0; j < vector_lyambda.size(); ++j) {
            if (dist(get_J(vector_lyambda.at(j), w, d), 0) < dist(min_result, 0)) {
                min_elem = j;
                min_alphas = alphas.size() - 1;
            }
        }
    }
    vector_alpha = alphas.at(min_alphas);
    vector_f_result = make_vector_f_result(vector_f_tild, alphas.at(min_alphas), K, M);
    min_lyambda = min_elem;
}

std::vector<double> make_vector_lyambda(std::vector<int> & l, int L) {
    std::vector<double> result;
    for (auto i : l) {
        result.push_back(double(i) / double(L));
    }
    return result;
}

void print(std::ofstream& outstream, const std::vector<double>& vector_x, const std::vector<double>& vector_y) {
    outstream << "\nprint...\n";
    for (auto i = 0; i < vector_x.size(); ++i) {
        outstream << '(' << vector_x.at(i) << ';' << vector_y.at(i) << ')';
    }
    outstream << std::endl;
}
