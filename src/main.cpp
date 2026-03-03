#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

const double rho = 2.0;

double eps_r(double x){
    if (x>= 0 && x <=1) return 1.0;
    else if (x > 1 && x <= 2) return 0.001;
    return 5.0;
}

// Funkcja bazowa e_i(x)
double e(int i, double x, double h) {
    if (x < (i - 1) * h || x > (i + 1) * h) return 0.0;
    if (x < i * h) return (x - (i - 1) * h) / h;
    return ((i + 1) * h - x) / h;
}

// Pochodna funkcji bazowej e'_i(x)
double e_prim(int i, double x, double h) {
    if (x < (i - 1) * h || x > (i + 1) * h) return 0.0;
    if (x < i * h) return 1.0 / h;
    return -1.0 / h;
}

// Całkowanie numeryczne (Kwadratura Gaussa-Legendre'a 2-punktowa)
double integrate_gauss(auto func, double a, double b) {
    if (a >= b) return 0.0;
    static const double nodes[2] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
    static const double weights[2] = {1.0, 1.0};
    
    double mid = (a + b) / 2.0;
    double half_len = (b - a) / 2.0;
    double area = 0.0;
    for(int i = 0; i < 2; ++i) {
        area += weights[i] * func(mid + half_len * nodes[i]);
    }
    return area * half_len;
}

double b_func(int i, int j, double h) {
    if (abs(i - j) > 1) return 0.0;
    
    double left = max({0.0, (i - 1) * h, (j - 1) * h});
    double right = min({(i + 1) * h, (j + 1) * h, 4.0});
    double mid = (i*h + j*h) / 2.0;
    
    auto integrand = [&](double x) { return e_prim(i, x, h) * e_prim(j, x, h); };
    double integral = 0.0;
    if (left < mid) {
        integral += integrate_gauss(integrand, left, mid);
    }
    if (mid < right) {
        integral += integrate_gauss(integrand, mid, right);
    }

    return e(i, 0, h) * e(j, 0, h) + integral;
}

double l_func(int j, double h) {
    double left = max(0.0, (j - 1) * h);
    double right = min((j + 1) * h, 4.0);
    
    auto integrand = [&](double x) { return (rho / eps_r(x)) * e(j, x, h); };
    double integral = integrate_gauss(integrand, left, j*h) + integrate_gauss(integrand, j*h, right);
    return -3.0 * e(j, 0, h) + integral;
}

// Algorytm Thomasa
vector<double> solveLinearSystem(vector<vector<double>>& B, vector<double>& L) {
    int n = L.size();
    vector<double> c_prime(n, 0.0);
    vector<double> d_prime(n, 0.0);
    vector<double> x(n, 0.0);

    c_prime[0] = B[0][1] / B[0][0];
    d_prime[0] = L[0] / B[0][0];

    for (int i = 1; i < n; i++) {
        double a = B[i][i-1];
        double b = B[i][i];
        double c = (i < n - 1) ? B[i][i+1] : 0.0;

        double m = b - a * c_prime[i-1];
        c_prime[i] = c / m;
        d_prime[i] = (L[i] - a * d_prime[i-1]) / m;
    }

    x[n-1] = d_prime[n-1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }

    return x;
}

string format(double x){
    string s = to_string(x);
    replace(s.begin(), s.end(), '.', ',');
    return s;
}

void saveResults(vector<double>& c, double h, int n) {
    ofstream file("wyniki.csv");
    file << "x;Phi\n";

    for (int i = 0; i < n; ++i) {
        double x_i = i * h; 
        double y_i = c[i] + 1; 
        file << format(x_i) << ";" << format(y_i) << "\n";
    }
    
    file << format(4.0) << ";" << format(1.0) << "\n";
    file.close();
}

int main() {
    int n;
    cout << "Podaj n: ";
    cin >> n;

    double h = 4.0 / n;
    vector<vector<double>> B(n, vector<double>(n));
    vector<double> L(n);

    // Budowanie macierzy 
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B[i][j] = b_func(j, i, h);
        }
        L[i] = l_func(i, h);
    }

    // Rozwiązanie
    vector<double> alphas = solveLinearSystem(B, L);

    // Wypisanie wyników 
    cout << fixed << setprecision(6);
    cout << "x\t\tphi(x)" << endl;
    cout << "------------------------" << endl;
    
    for (int i = 0; i < n; i++) {
        double x_val = i * h;
        double y_val = alphas[i] + 1.0; 
        cout << x_val << "\t" << y_val << endl;
    }
    // Ostatni punkt z warunku Dirichleta
    cout << 4.0 << "\t" << 1.0 << endl;

    saveResults(alphas, h, n);

    return 0;
}
