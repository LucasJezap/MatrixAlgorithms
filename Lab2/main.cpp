#include <iostream>
#include <random>
#include <iomanip>
#include <chrono>

#define TIMING

#ifdef TIMING
#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now();
#define START_TIMER  start = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name)  std::cout << "RUNTIME of " << name << ": " << \
    std::chrono::duration_cast<std::chrono::milliseconds>( \
            std::chrono::high_resolution_clock::now()-start \
    ).count() << " ms " << std::endl;
#else
#define INIT_TIMER
#define START_TIMER
#define STOP_TIMER(name)
#endif

bool show_matrixes;

void show_matrix(double **M, int m, int n) {
    if (!show_matrixes)
        return;
    std::cout << std::fixed;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::setprecision(3) << std::setw(10) << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

double **create_matrix(int m, int n, const std::string &gen) {
    auto **matrix = new double *[m];
    std::default_random_engine  generator{static_cast<std::bernoulli_distribution::result_type>(std::chrono::system_clock::now().time_since_epoch().count())};
    std::uniform_int_distribution<int> distribution(1,100);
    for (int i = 0; i < m; i++)
        matrix[i] = new double[n];

    if (gen != "y") {
        std::cout << "Please enter consecutive matrix rows "
                     "(the size of the matrix is " << m << "x" << n << ")" << std::endl;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                std::cin >> matrix[i][j];
    } else {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                matrix[i][j] = distribution(generator) - distribution(generator) % 100;
    }
    std::cout << "Generated matrix:" << std::endl;
    show_matrix(matrix, m, n);
    return matrix;
}

double **gauss_elimination(double **A, int m, int n) {
    INIT_TIMER
    START_TIMER
    std::cout << "Step " << 0 << ":" << std::endl;
    show_matrix(A, m, n);

    for (int k = 0; k < m - 1; k++) {
        int ind = k;
        for (int j = k + 1; j < m; j++) {
            if (std::abs(A[j][k]) > std::abs(A[ind][k]))
                ind = j;
        }
        for (int p = k; p < n; p++)
            std::swap(A[k][p], A[ind][p]);

        std::cout << "Step " << k + 1 << ":" << std::endl;
        std::cout << "After pivoting:" << std::endl;
        show_matrix(A, m, n);

        double a_kk = A[k][k];
        for (int j = k + 1; j < m; j++) {
            double ratio = A[j][k] / a_kk;
            for (int p = k; p < n; p++) {
                A[j][p] -= 1.0 * ratio * A[k][p];
            }
        }
        std::cout << "After elimination" << std::endl;
        show_matrix(A, m, n);
    }
    STOP_TIMER("the Gauss elimination (with pivoting)")

    std::cout << "Matrix A after the Gauss elimination (with pivoting):" << std::endl;
    show_matrix(A, m, n);

    return A;
}

double **raw_gauss_elimination(double **A, int m, int n) {
    INIT_TIMER
    START_TIMER

    for (int k = 0; k < m - 1; k++) {
        int ind = k;
        for (int j = k + 1; j < m; j++) {
            if (std::abs(A[j][k]) > std::abs(A[ind][k]))
                ind = j;
        }
        for (int p = k; p < n; p++)
            std::swap(A[k][p], A[ind][p]);

        double a_kk = A[k][k];
        for (int j = k + 1; j < m; j++) {
            double ratio = A[j][k] / a_kk;
            for (int p = k; p < n; p++) {
                A[j][p] -= 1.0 * ratio * A[k][p];
            }
        }
    }
    STOP_TIMER("the Gauss elimination (with pivoting)")

    std::cout << "Matrix A after the Gauss elimination (with pivoting):" << std::endl;
    show_matrix(A, m, n);

    return A;
}

int main() {
    std::cout << "Do you want to show matrixes (problematic with big ones)? (y/n)" << std::endl;
    std::string gen;
    std::cin >> gen;
    if (gen == "y")
        show_matrixes = true;
    else
        show_matrixes = false;
    std::cout << "Do you want to generate the matrixes? (y/n)" << std::endl;
    std::cin >> gen;
    if (gen != "y" && gen != "n") {
        std::cout << "Wrong answer";
        return -1;
    }
    int m, n;
    std::cout << "Matrix A:" << std::endl;
    std::cout << "Please enter number of rows: " << std::endl;
    std::cin >> m;
    std::cout << "Please enter number of columns: " << std::endl;
    std::cin >> n;

    double **A = create_matrix(m, n, gen);

    std::cout << "Do you want to see steps? (y/n)" << std::endl;
    std::cin >> gen;
    if (gen == "y")
        A = gauss_elimination(A, m, n);
    else
        A = raw_gauss_elimination(A, m, n);

    for (int i = 0; i < m; i++)
        delete A[i];
    delete A;

    return 0;
}
