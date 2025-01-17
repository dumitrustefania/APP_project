#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <chrono>
#include <filesystem>
#include <mpi.h>
#include <omp.h>

using namespace std;

void print_matrix(float **matrix, int size, ofstream &output_file);
void LU_decomposition(float **a, float **l, float **u, int size); // Only threads
void forward_substitution(float **L, float *b, float *y, int size);
void backward_substitution(float **U, float *y, float *x, int size);
void compute_inverse(float **a, float **l, float **u, float **a_inv, int size, int rank, int num_procs); // MPI + threads
string get_output_filename(const string &input_filename, const string &matrix_type);

// Forward substitution (sequential for each vector)
void forward_substitution(float **L, float *b, float *y, int size)
{
    for (int i = 0; i < size; ++i)
    {
        y[i] = b[i];
        for (int j = 0; j < i; ++j)
        {
            y[i] -= L[i][j] * y[j];
        }
        y[i] /= L[i][i];
    }
}

// Backward substitution (sequential for each vector)
void backward_substitution(float **U, float *y, float *x, int size)
{
    for (int i = size - 1; i >= 0; --i)
    {
        x[i] = y[i];
        for (int j = i + 1; j < size; ++j)
        {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}

// Compute inverse using MPI + threads
// Distribute columns among ranks (MPI), each rank uses OpenMP to process its subset of columns in parallel
void compute_inverse(float **a, float **l, float **u, float **a_inv, int size, int rank, int num_procs)
{
    int cols_per_proc = size / num_procs;
    int remainder = size % num_procs;

    int start_col = rank * cols_per_proc + min(rank, remainder);
    int num_cols = cols_per_proc + (rank < remainder ? 1 : 0);

// Parallelize the loop over columns this rank is responsible for
#pragma omp parallel
    {
        float *y = new float[size];
        float *x = new float[size];

#pragma omp for schedule(static)
        for (int i = start_col; i < start_col + num_cols; ++i)
        {
            float *e = new float[size]();
            e[i] = 1.0f;

            // Each column solve is sequential but handled by different threads in parallel
            forward_substitution(l, e, y, size);
            backward_substitution(u, y, x, size);

            for (int j = 0; j < size; ++j)
            {
                a_inv[j][i] = x[j];
            }
            delete[] e;
        }

        delete[] y;
        delete[] x;
    }

    // Send columns to root if rank != 0
    if (rank != 0)
    {
        for (int i = 0; i < size; ++i)
        {
            MPI_Send(&a_inv[i][start_col], num_cols, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        // Root receives data from other ranks
        for (int p = 1; p < num_procs; ++p)
        {
            int p_start_col = p * cols_per_proc + min(p, remainder);
            int p_num_cols = cols_per_proc + (p < remainder ? 1 : 0);
            for (int i = 0; i < size; ++i)
            {
                MPI_Recv(&a_inv[i][p_start_col], p_num_cols, MPI_FLOAT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
}

// LU decomposition using only threads (on rank 0)
void LU_decomposition(float **a, float **l, float **u, int size)
{
    for (int i = 0; i < size; i++)
    {
// Compute L column i
#pragma omp parallel for schedule(static)
        for (int j = 0; j < size; j++)
        {
            if (j < i)
            {
                l[j][i] = 0.0f;
            }
            else
            {
                float sum = a[j][i];
                for (int k = 0; k < i; k++)
                {
                    sum -= l[j][k] * u[k][i];
                }
                l[j][i] = sum;
            }
        }

// Compute U row i
#pragma omp parallel for schedule(static)
        for (int j = 0; j < size; j++)
        {
            if (j < i)
            {
                u[i][j] = 0.0f;
            }
            else if (j == i)
            {
                u[i][j] = 1.0f;
            }
            else
            {
                float sum = a[i][j];
                for (int k = 0; k < i; k++)
                {
                    sum -= (l[i][k] * u[k][j]);
                }
                u[i][j] = sum / l[i][i];
            }
        }
    }
}

// Printing and filename generation
void print_matrix(float **matrix, int size, ofstream &output_file)
{
    output_file << size << "\n";
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            output_file << setw(10) << matrix[i][j] << " ";
        }
        output_file << "\n";
    }
}

string get_output_filename(const string &input_filename, const string &matrix_type)
{
    string dir = "outputs/hybrid/" + input_filename.substr(7, input_filename.find(".") - 7);
    filesystem::create_directories(dir);
    return dir + "/" + matrix_type + ".out";
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Set number of OpenMP threads
    omp_set_num_threads(stoi(argv[2]));

    if (argc < 2)
    {
        if (rank == 0)
            cerr << "Usage: " << argv[0] << " <input_file>" << '\n';
        MPI_Finalize();
        return 1;
    }

    double start = MPI_Wtime();

    string input_filename = argv[1];

    int size;
    float **a;
    if (rank == 0)
    {
        ifstream input_file(input_filename);
        if (!input_file.is_open())
        {
            cerr << "Error opening input file." << '\n';
            MPI_Finalize();
            return 1;
        }

        input_file >> size;
        a = new float *[size];
        for (int i = 0; i < size; i++)
        {
            a[i] = new float[size];
            for (int j = 0; j < size; j++)
            {
                input_file >> a[i][j];
            }
        }
        input_file.close();
    }

    // Broadcast matrix size to all ranks
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // For non-root ranks, allocate a(matrix a)
    if (rank != 0)
    {
        a = new float *[size];
        for (int i = 0; i < size; i++)
        {
            a[i] = new float[size];
        }
    }

    // Broadcast A
    for (int i = 0; i < size; i++)
    {
        MPI_Bcast(a[i], size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    float **l = new float *[size];
    float **u = new float *[size];
    float **a_inv = new float *[size];

#pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        l[i] = new float[size];
        u[i] = new float[size];
        a_inv[i] = new float[size];
    }

    // Only rank 0 performs LU decomposition using threads only
    double t1, t2;
    if (rank == 0)
    {
        t1 = MPI_Wtime();
        LU_decomposition(a, l, u, size); // Threads only, no MPI inside
        t2 = MPI_Wtime();
        cout << "LU time: " << t2 - t1 << '\n';
    }

    // Broadcast L and U to all other ranks after LU decomposition is done
    // No need for barrier since broadcast is blocking
    for (int i = 0; i < size; i++)
    {
        MPI_Bcast(l[i], size, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(u[i], size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    // Compute inverse using MPI + threads
    t1 = MPI_Wtime();
    compute_inverse(a, l, u, a_inv, size, rank, num_procs);
    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();

    if (rank == 0)
    {
cout << "Compute inverse time: " << t2 - t1 << '\n';

    // Write matrices to output files
#pragma omp parallel sections
        {
#pragma omp section
            {
                ofstream l_out(get_output_filename(input_filename, "L"));
                print_matrix(l, size, l_out);
            }
#pragma omp section
            {
                ofstream u_out(get_output_filename(input_filename, "U"));
                print_matrix(u, size, u_out);
            }
#pragma omp section
            {
                ofstream a_inv_out(get_output_filename(input_filename, "A_inv"));
                print_matrix(a_inv, size, a_inv_out);
            }
        }

        // Cleanup
#pragma omp parallel for schedule(static)
        for (int i = 0; i < size; i++)
        {
            delete[] a[i];
            delete[] l[i];
            delete[] u[i];
            delete[] a_inv[i];
        }
        delete[] a;
        delete[] l;
        delete[] u;
        delete[] a_inv;
    }

    double end = MPI_Wtime();
    if (rank == 0)
        cout << "Total time: " << end - start << '\n';

    MPI_Finalize();
    return 0;
}
