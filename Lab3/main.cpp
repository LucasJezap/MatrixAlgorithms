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

struct Csc_format {
    int nnz;
    int n;
    int *irn;
    int *colptr;
    double *val;
};

struct Csr_format {
    int nnz;
    int n;
    int *icl;
    int *rowptr;
    double *val;
};

bool show_matrixes;

Csr_format create_new_csr_format(int nnz, int n) {
    Csr_format format{};
    format.nnz = nnz;
    format.n = n;
    format.icl = new int[nnz];
    format.val = new double[nnz];
    format.rowptr = new int[n + 1];

    return format;
}

void print_csr_format(Csr_format format) {
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

Csr_format create_csr_from_csc_format(Csc_format csc_format) {
    INIT_TIMER
    START_TIMER
    Csr_format csr_format = create_new_csr_format(csc_format.nnz, csc_format.n);

    // first i get rowptr array -> it's fairly easy
    for (int i = 0; i < csc_format.n + 1; i++)
        csr_format.rowptr[i] = 0;
    for (int i = 0; i < csc_format.n; i++) {
        csr_format.rowptr[csc_format.irn[i]]++;
    }
    csr_format.rowptr[0] = 1;
    for (int i = 1; i < csc_format.n + 1; i++) {
        csr_format.rowptr[i] += csr_format.rowptr[i - 1];
    }

    // now i can read values by each column (CSC) and put them in appropriate place
    int row_indexes[csc_format.n];
    row_indexes[0] = 0;
    for (int i = 1; i < csc_format.n; i++) {
        row_indexes[i] = row_indexes[i - 1] + (csr_format.rowptr[i] - csr_format.rowptr[i - 1]);
    }

    int column = 1;
    int count = 1;
    for (int i = 0; i < csc_format.nnz; i++) {
        int row = csc_format.irn[i];
        csr_format.val[row_indexes[row - 1]] = csc_format.val[i];
        csr_format.icl[row_indexes[row - 1]] = column;
        if (++count == csc_format.colptr[column])
            column++;
        row_indexes[row - 1]++;
    }

    STOP_TIMER("CSC to CSR conversion")
    return csr_format;
}

Csc_format create_new_csc_format(int nnz, int n) {
    Csc_format format{};
    format.nnz = nnz;
    format.n = n;
    format.irn = new int[nnz];
    format.val = new double[nnz];
    format.colptr = new int[n + 1];

    return format;
}

Csc_format create_csc_format_from_matrix(double **A, int a, int b) {
    int nnz = 0;
    int n = b;
    for (int i = 0; i < a; i++)
        for (int j = 0; j < b; j++)
            if (A[i][j] != 0)
                nnz++;
    Csc_format format = create_new_csc_format(nnz, n);

    int idx1 = 0, idx2 = 0;
    format.colptr[idx2++] = 1;
    for (int j = 0; j < b; j++) {
        for (int i = 0; i < a; i++) {
            if (A[i][j] != 0) {
                format.val[idx1] = A[i][j];
                format.irn[idx1++] = i + 1;
            }
        }
        format.colptr[idx2++] = idx1 + 1;
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

double **create_matrix_from_csc_format(Csc_format format) {
    auto **A = new double *[format.n];
    for (int i = 0; i < format.n; i++) {
        A[i] = new double[format.n];
        for (int j = 0; j < format.n; j++)
            A[i][j] = 0;
    }

    int idx = 0;
    for (int column = 0; column < format.n; column++) {
        int rows = format.colptr[column + 1] - format.colptr[column];
        for (int i = 0; i < rows; i++) {
            A[format.irn[idx] - 1][column] = format.val[idx];
            idx++;
        }
    }

    return A;
}

void print_csc_format(Csc_format format) {
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "CSC FORMAT" << std::endl;
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
    std::cout << ")\nCOLPTR = (";
    for (int i = 0; i < format.n + 1; i++) {
        std::cout << format.colptr[i];
        if (i != format.n) std::cout << ",";
    }
    std::cout << ")\n---------------------------------------------------\n";
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
                if (distribution(generator) > percentage)
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
    std::cout << "|    CSC to CSR format calculator   |" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
}

int main() {
    hello();
    std::cout << "Do you want to enter CSC matrix format or generate it from given matrix? (1/2)" << std::endl;
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

        Csc_format csc_format = create_csc_format_from_matrix(A, m, n);
        print_csc_format(csc_format);

        Csr_format csr_format = create_csr_from_csc_format(csc_format);
        print_csr_format(csr_format);

        for (int i = 0; i < m; i++)
            delete A[i];
        delete A;
    } else {
        int nnz, n;
        std::cout << "Please enter NNZ (number of non-zeros): " << std::endl;
        std::cin >> nnz;
        std::cout << "Please enter n (size of matrix): " << std::endl;
        std::cin >> n;
        Csc_format csc_format = create_new_csc_format(nnz, n);
        std::cout << "Please enter " << nnz << " IRN values (indices of rows)" << std::endl;
        for (int i = 0; i < csc_format.nnz; i++) {
            std::cin >> csc_format.irn[i];
        }
        std::cout << "Please enter " << n + 1 << " COLPTR values (extents of columns)" << std::endl;
        for (int i = 0; i < csc_format.n + 1; i++) {
            std::cin >> csc_format.colptr[i];
        }
        std::cout << "Please enter " << nnz << " VAL values (values of consecutive elements)" << std::endl;
        for (int i = 0; i < csc_format.nnz; i++) {
            std::cin >> csc_format.val[i];
        }
        double **A = create_matrix_from_csc_format(csc_format);
        show_matrix(A, n, n);

        print_csc_format(csc_format);
        Csr_format csr_format = create_csr_from_csc_format(csc_format);
        print_csr_format(csr_format);

        for (int i = 0; i < n; i++)
            delete A[i];
        delete A;
    }

    return 0;
}
