#include <iostream>
#include <random>
#include <iomanip>
#include <chrono>

#define TIMING

#ifdef TIMING
#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now();
#define START_TIMER  start = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name)  std::cout << "RUNTIME of " << (name) << ": " << \
    std::chrono::duration_cast<std::chrono::milliseconds>( \
            std::chrono::high_resolution_clock::now()-start \
    ).count() << " ms " << std::endl;
#else
#define INIT_TIMER
#define START_TIMER
#define STOP_TIMER(name)
#endif

struct Coordinate_format {
    int nnz;
    int n;
    std::vector<int> irn;
    std::vector<int> jcn;
    std::vector<double> val;
};

struct Csr_format {
    int nnz;
    int n;
    std::vector<int> icl;
    std::vector<int> rowptr;
    std::vector<double> val;
};

bool show_matrixes;

Coordinate_format create_new_coordinate_format(int nnz, int n) {
    Coordinate_format format{};
    format.nnz = nnz;
    format.n = n;

    return format;
}

Csr_format create_new_csr_format(int nnz, int n) {
    Csr_format format{};
    format.nnz = nnz;
    format.n = n;
    for (int i = 0; i < nnz; i++) {
        format.icl.push_back(0);
        format.val.push_back(0);
    }
    for (int i = 0; i < n + 1; i++) {
        format.rowptr.push_back(0);
    }

    return format;
}

void print_coordinate_format(Coordinate_format format) {
    if (!show_matrixes)
        return;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "COORDINATE FORMAT" << std::endl;
    std::cout << "NNZ = " << format.nnz << std::endl;
    std::cout << "N = " << format.n << std::endl;
    std::cout << "VAL = (";
    for (int i = 0; i < format.nnz; i++) {
        std::cout << format.val[i];
        if (i != format.nnz - 1) std::cout << ",";
    }
    std::cout << ")\nIRN = (";
    for (int i = 0; i < format.nnz; i++) {
        std::cout << format.irn[i];
        if (i != format.nnz - 1) std::cout << ",";
    }
    std::cout << ")\nJCN = (";
    for (int i = 0; i < format.nnz; i++) {
        std::cout << format.jcn[i];
        if (i != format.nnz - 1) std::cout << ",";
    }
    std::cout << ")\n---------------------------------------------------";
}

void print_csr_format(Csr_format format) {
    if (!show_matrixes)
        return;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "CSR FORMAT" << std::endl;
    std::cout << "NNZ = " << format.nnz << std::endl;
    std::cout << "N = " << format.n << std::endl;
    std::cout << "VAL = (";
    for (int i = 0; i < format.nnz; i++) {
        std::cout << format.val[i];
        if (i != format.nnz - 1) std::cout << ",";
    }
    std::cout << ")\nICL = (";
    for (int i = 0; i < format.nnz; i++) {
        std::cout << format.icl[i];
        if (i != format.nnz - 1) std::cout << ",";
    }
    std::cout << ")\nROWPTR = (";
    for (int i = 0; i < format.n + 1; i++) {
        std::cout << format.rowptr[i];
        if (i != format.n) std::cout << ",";
    }
    std::cout << ")\n---------------------------------------------------";
}

Csr_format create_csr_format_from_matrix(double **A, int a, int b) {
    int nnz = 0;
    int n = b;
    for (int i = 0; i < a; i++)
        for (int j = 0; j < b; j++)
            if (A[i][j] != 0)
                nnz++;
    Csr_format format = create_new_csr_format(nnz, n);

    int idx1 = 0, idx2 = 0;
    format.rowptr[idx2++] = 1;
    for (int i = 0; i < a; i++) {
        for (int j = 0; j < b; j++) {
            if (A[i][j] != 0) {
                format.val[idx1] = A[i][j];
                format.icl[idx1++] = j + 1;
            }
        }
        format.rowptr[idx2++] = idx1 + 1;
    }

    return format;
}

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

double **create_matrix_from_csr_format(Csr_format format) {
    auto **A = new double *[format.n];
    for (int i = 0; i < format.n; i++) {
        A[i] = new double[format.n];
        for (int j = 0; j < format.n; j++)
            A[i][j] = 0;
    }

    int idx = 0;
    for (int row = 0; row < format.n; row++) {
        int columns = format.rowptr[row + 1] - format.rowptr[row];
        for (int i = 0; i < columns; i++) {
            A[row][format.icl[idx] - 1] = format.val[idx];
            idx++;
        }
    }

    return A;
}

double **create_matrix(int m, int n, const std::string &gen, int percentage) {
    auto **matrix = new double *[m];
    std::default_random_engine generator{
            static_cast<std::bernoulli_distribution::result_type>(std::chrono::system_clock::now().time_since_epoch().count())};
    std::uniform_int_distribution<int> distribution(1, 100);
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
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    while ((matrix[i][j] = distribution(generator) - distribution(generator)) == 0);
                } else if (distribution(generator) > percentage)
                    matrix[i][j] = 0;
                else
                    matrix[i][j] = distribution(generator) - distribution(generator) % 100;
            }
    }
    std::cout << "Generated matrix:" << std::endl;
    show_matrix(matrix, m, n);
    return matrix;
}


void hello() {
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "|       CSR LU DECOMPOSITION        |" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
}

void lu_decomposition(double **M, int a, int b) {
    int n = a;
    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++) {
            if (M[k][k] != 0 && M[i][k] != 0)
                M[i][k] /= M[k][k];
        }
        for (int j = k + 1; j < n; j++) {
            for (int i = k + 1; i < n; i++) {
                if (M[k][j] != 0 && M[i][k] != 0)
                    M[i][j] -= M[i][k] * M[k][j];
            }
        }
    }
    std::cout << "\nLU Decomposition (in matrix A): " << std::endl;
    show_matrix(M, a, b);
}

Coordinate_format create_coordinate_from_csr(Csr_format csr_format) {   // złożoność O(NNZ^2)
    Coordinate_format format = create_new_coordinate_format(csr_format.nnz, csr_format.n);

    for (int column = 0; column < csr_format.n; column++) {
        for (int i = 0; i < csr_format.n; i++) {
            for (int j = csr_format.rowptr[i]; j < csr_format.rowptr[i + 1]; j++) {
                if (csr_format.icl[j - 1] == column + 1) {
                    format.irn.push_back(i + 1);
                    format.jcn.push_back(column + 1);
                    format.val.push_back(csr_format.val[j - 1]);
                }
            }
        }
    }

    return format;
}

Csr_format create_csr_from_coordinate(Coordinate_format coordinate_format) {   // złożoność O(NNZ^2)
    Csr_format format = create_new_csr_format(coordinate_format.nnz, coordinate_format.n);

    int idx = 0;
    format.rowptr[0] = 1;
    for (int row = 0; row < coordinate_format.n; row++) {
        for (int i = 0; i < coordinate_format.nnz; i++) {
            if (coordinate_format.irn[i] == row + 1) {
                format.icl[idx] = coordinate_format.jcn[i];
                format.val[idx++] = coordinate_format.val[i];
            }
        }
        format.rowptr[row + 1] = idx + 1;
    }

    return format;
}

void lu_decomposition_csr(Csr_format csr_format) {     // O(N * NNZ^2)
    INIT_TIMER
    START_TIMER
    Coordinate_format coordinate_format = create_coordinate_from_csr(std::move(csr_format));

    for (int k = 0; k < coordinate_format.n; k++) {       // O(N)
        double a_kk = 0.0;
        for (int i = 0; i < coordinate_format.nnz and coordinate_format.jcn[i] <= k + 1; i++) {
            if (coordinate_format.irn[i] == k + 1 && coordinate_format.jcn[i] == k + 1) {
                a_kk = coordinate_format.val[i];
                continue;
            }
            if (a_kk != 0.0) {
                coordinate_format.val[i] /= (1.0 * a_kk);
            }                                                   // O(NNZ)
        }
        for (int i = 0; i < coordinate_format.nnz and coordinate_format.jcn[i] <= k + 1; i++) {  // O(NNZ)
            // k-ta kolumna
            if (coordinate_format.jcn[i] == k + 1 && coordinate_format.irn[i] > k + 1) {
                int j = i;
                int col = coordinate_format.jcn[i];
                while (j < coordinate_format.nnz) {       // O(NNZ)
                    while (coordinate_format.jcn[j] <= col) { j++; }  // pierwszy element w nowej kolumnie
                    col = coordinate_format.jcn[j];
                    while (coordinate_format.jcn[j] == col &&   
                           coordinate_format.irn[j] < k + 1) { j++; } // element w k-tym wierszu
                    // sprawdzam czy to odpowiedni element - wiersz k
                    if (coordinate_format.irn[j] == k + 1) {
                        int tmp = j;
                        // szukam czy istnieje niezerowy element na przecięciu k x k
                        // jeśli nie - trzeba go stworzyć
                        while (coordinate_format.jcn[j] == col &&
                               coordinate_format.irn[j] < coordinate_format.irn[i]) { j++; }
                        // jeśli to ten element
                        if (coordinate_format.jcn[j] == col &&
                            coordinate_format.irn[j] == coordinate_format.irn[i]) {
                            coordinate_format.val[j] -= 1.0 * coordinate_format.val[i] * coordinate_format.val[tmp];
                        } // w przeciwnym razie
                        else {
                            coordinate_format.nnz++;
                            coordinate_format.irn.insert(coordinate_format.irn.begin() + j,
                                                         coordinate_format.irn[i]);
                            coordinate_format.jcn.insert(coordinate_format.jcn.begin() + j,
                                                         coordinate_format.jcn[tmp]);
                            coordinate_format.val.insert(coordinate_format.val.begin() + j,
                                                         -1.0 * coordinate_format.val[i] * coordinate_format.val[tmp]);
                        }
                    }
                }

            }
        }
    }
    STOP_TIMER("LU FACTORIZATION (FROM CSR)")
    Csr_format format = create_csr_from_coordinate(coordinate_format);
    std::cout << "\nLU Decomposition using CSR Format (in matrix A): " << std::endl;
    show_matrix(create_matrix_from_csr_format(format), format.n, format.n);
    print_csr_format(format);
}

int main() {
    hello();
    std::cout << "Do you want to enter CSR matrix format or generate it from given matrix? (1/2)" << std::endl;
    std::string type;
    bool from_matrix;
    std::cin >> type;
    if (type == "2")
        from_matrix = true;
    else
        from_matrix = false;
    std::cout << "Do you want to show matrixes (problematic with big ones)? (y/n)" << std::endl;
    std::string gen;
    std::cin >> gen;
    if (gen == "y")
        show_matrixes = true;
    else
        show_matrixes = false;

    if (from_matrix) {
        std::cout << "Do you want to generate the matrixes? (y/n)" << std::endl;
        std::cin >> gen;
        if (gen != "y" && gen != "n") {
            std::cout << "Wrong answer";
            return -1;
        }
        int percentage = 100;
        if (gen == "y") {
            std::cout << "Please enter the percentage of non-zero elements in matrix (0-100)" << std::endl;
            std::cin >> percentage;
        }
        int m, n;
        std::cout << "Matrix A:" << std::endl;
        std::cout << "Please enter number of rows: " << std::endl;
        std::cin >> m;
        std::cout << "Please enter number of columns: " << std::endl;
        std::cin >> n;

        double **A = create_matrix(m, n, gen, percentage);

        Csr_format csr_format = create_csr_format_from_matrix(A, m, n);
        print_csr_format(csr_format);

        lu_decomposition(A, m, n);
        lu_decomposition_csr(csr_format);

        for (int i = 0; i < m; i++)
            delete A[i];
        delete A;
    } else {
        int nnz, n;
        std::cout << "Please enter NNZ (number of non-zeros): " << std::endl;
        std::cin >> nnz;
        std::cout << "Please enter n (size of matrix): " << std::endl;
        std::cin >> n;
        Csr_format csr_format = create_new_csr_format(nnz, n);
        std::cout << "Please enter " << nnz << " ICL values (indices of columns)" << std::endl;
        for (int i = 0; i < csr_format.nnz; i++) {
            std::cin >> csr_format.icl[i];
        }
        std::cout << "Please enter " << n + 1 << " ROWPTR values (extents of rows)" << std::endl;
        for (int i = 0; i < csr_format.n + 1; i++) {
            std::cin >> csr_format.rowptr[i];
        }
        std::cout << "Please enter " << nnz << " VAL values (values of consecutive elements)" << std::endl;
        for (int i = 0; i < csr_format.nnz; i++) {
            std::cin >> csr_format.val[i];
        }

        double **A = create_matrix_from_csr_format(csr_format);
        std::cout << "Matrix A:" << std::endl;
        show_matrix(A, n, n);

        print_csr_format(csr_format);
        Coordinate_format coordinate_format = create_coordinate_from_csr(csr_format);
        print_coordinate_format(coordinate_format);

        lu_decomposition(A, n, n);
        lu_decomposition_csr(csr_format);

        for (int i = 0; i < n; i++)
            delete A[i];
        delete A;
    }

    return 0;
}
