#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>
#include <ctime>
#include <mpi.h>
#include <filesystem>

using namespace std;

// Define the parallelization threshold
const int PARALLEL_THRESHOLD = 500; 

void print_matrix(float **matrix, int size, ofstream &output_file);
void LU_decomposition(float **a, float **l, float **u, int size, int rank, int num_procs);
void forward_substitution(float **L, float *b, float *y, int size);
void backward_substitution(float **U, float *y, float *x, int size);
void compute_inverse(float **a, float **l, float **u, float **a_inv, int size, int rank, int num_procs);
string get_output_filename(const string &input_filename, const string &matrix_type);

// Forward substitution function (L * y = b)
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

// Backward substitution function (U * x = y)
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

// Function to compute the inverse of matrix 'a' using LU decomposition
void compute_inverse(float **a, float **l, float **u, float **a_inv, int size, int rank, int num_procs)
{
    float *y = new float[size]; // Temporary array for forward substitution
    float *x = new float[size]; // Temporary array for backward substitution

    // Compute the number of columns per process
    int cols_per_proc = size / num_procs;
    int remainder = size % num_procs;

    int start_col = rank * cols_per_proc + min(rank, remainder);
    int num_cols = cols_per_proc + (rank < remainder ? 1 : 0);
    int end_col = start_col + num_cols - 1;

    // Each process computes its assigned columns
    for (int i = start_col; i <= end_col; ++i)
    {
        float *e = new float[size](); // Unit vector with 1 at position i
        e[i] = 1;
        forward_substitution(l, e, y, size);  // Solve L * y = e
        backward_substitution(u, y, x, size); // Solve U * x = y
        for (int j = 0; j < size; ++j)
        {
            a_inv[j][i] = x[j]; // Store the solution x
        }
        delete[] e;
    }

    // Send the computed columns to the root process
    if (rank != 0)
    {
        for (int i = 0; i < size; ++i)
        {
            MPI_Send(&a_inv[i][start_col], num_cols, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        // Receive data from other processes
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

    delete[] y;
    delete[] x;
}


// LU decomposition function (L * U = A)
void LU_decomposition(float **a, float **l, float **u, int size, int rank, int num_procs)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            // Compute the lower triangular matrix L
            if (j < i)
            {
                if (rank == 0)
                    l[j][i] = 0; // Set the lower triangular matrix L
            }
            else
            {
                if (rank == 0)
                    l[j][i] = a[j][i]; // Copy elements from A to L

                // Synchronize before starting the loop
                MPI_Bcast(&l[j][i], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

                int count = i;

                if (count > PARALLEL_THRESHOLD)
                {
                    // Parallelize the inner loop over k
                    float sum = 0.0;
                    int chunk_size = count / num_procs;
                    int remainder = count % num_procs;
                    int start_k = rank * chunk_size + min(rank, remainder);
                    int end_k = start_k + chunk_size - 1;
                    if (rank < remainder)
                        end_k += 1;

                    for (int k = start_k; k <= end_k && k < i; k++)
                    {
                        sum += l[j][k] * u[k][i]; // Compute partial sum
                    }

                    // Reduce the partial sums to calculate the total sum
                    float total_sum = 0.0;
                    MPI_Allreduce(&sum, &total_sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

                    if (rank == 0)
                        l[j][i] -= total_sum; // Subtract the total sum from l[j][i]
                }
                else
                {
                    // Execute sequentially on root
                    if (rank == 0)
                    {
                        for (int k = 0; k < i; k++)
                        {
                            l[j][i] -= l[j][k] * u[k][i];
                        }
                    }
                }

                // Broadcast the updated l[j][i] to all processes
                MPI_Bcast(&l[j][i], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
            }

            // Compute the upper triangular matrix U
            if (j < i)
            {
                if (rank == 0)
                    u[i][j] = 0; // Set the upper triangular matrix U
            }
            else if (j == i)
            {
                if (rank == 0)
                    u[i][j] = 1; // Diagonal of U is set to 1 (since L is lower triangular)
            }
            else
            {
                if (rank == 0)
                    u[i][j] = a[i][j] / l[i][i];

                // Synchronize before starting the loop
                MPI_Bcast(&u[i][j], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

                int count = i;

                if (count > PARALLEL_THRESHOLD)
                {
                    // Parallelize the inner loop over k
                    float sum = 0.0;
                    int chunk_size = count / num_procs;
                    int remainder = count % num_procs;
                    int start_k = rank * chunk_size + min(rank, remainder);
                    int end_k = start_k + chunk_size - 1;
                    if (rank < remainder)
                        end_k += 1;

                    for (int k = start_k; k <= end_k && k < i; k++)
                    {
                        sum += (l[i][k] * u[k][j]) / l[i][i]; // Compute partial sum
                    }

                    // Reduce the partial sums to calculate the total sum
                    float total_sum = 0.0;
                    MPI_Allreduce(&sum, &total_sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

                    if (rank == 0)
                        u[i][j] -= total_sum; // Subtract the total sum from u[i][j]
                }
                else
                {
                    // Execute sequentially on root
                    if (rank == 0)
                    {
                        for (int k = 0; k < i; k++)
                        {
                            u[i][j] -= (l[i][k] * u[k][j]) / l[i][i];
                        }
                    }
                }

                // Broadcast the updated u[i][j] to all processes
                MPI_Bcast(&u[i][j], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
            }
        }
    }
}

// Function to print the matrix in a readable format to the output file
void print_matrix(float **matrix, int size, ofstream &output_file)
{
    output_file << size << '\n';
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            output_file << setw(10) << matrix[i][j] << " ";
        }
        output_file << '\n';
    }
}

// Function to generate the output filename for each matrix (L, U, A_inv)
string get_output_filename(const string &input_filename, const string &matrix_type)
{
    string dir = "outputs/mpi/" + input_filename.substr(7, input_filename.find(".") - 7);
    filesystem::create_directories(dir); // Create the directory if it doesn't exist
    return dir + "/" + matrix_type + ".out";
}

int main(int argc, char *argv[])
{
    double start, end;

    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    start = MPI_Wtime();
    
    if (argc < 2)
    {
        if (rank == 0)
            cerr << "Usage: " << argv[0] << " <input_size>" << '\n';
        MPI_Finalize();
        return 1;
    }

    string input_size = argv[1];

    string input_filename = argv[1];
    ifstream input_file(input_filename);

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

        input_file >> size; // Read the size of the matrix

        // Allocate and read matrix a
        a = new float *[size];
        for (int i = 0; i < size; ++i)
        {
            a[i] = new float[size];
            for (int j = 0; j < size; ++j)
            {
                input_file >> a[i][j];
            }
        }
        input_file.close();
    }

    // Broadcast the size to all processes
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Now allocate a on other processes
    if (rank != 0)
    {
        a = new float *[size];
        for (int i = 0; i < size; ++i)
        {
            a[i] = new float[size];
        }
    }

    // Now broadcast the matrix a
    for (int i = 0; i < size; ++i)
    {
        MPI_Bcast(a[i], size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    // Allocate l, u, a_inv
    float **l = new float *[size];
    float **u = new float *[size];
    float **a_inv = new float *[size];
    for (int i = 0; i < size; ++i)
    {
        l[i] = new float[size];
        u[i] = new float[size];
        a_inv[i] = new float[size];
    }

    // Measure the runtimes
    double t1, t2;
    
    t1 = MPI_Wtime();

    // Perform LU decomposition
    LU_decomposition(a, l, u, size, rank, num_procs);
    MPI_Barrier(MPI_COMM_WORLD);

    t2 = MPI_Wtime(); // Calculate elapsed time

    if (rank == 0)
        cout << "LU time: " << t2 - t1 << '\n'; 


    // bcast l and u to all processes
    for (int i = 0; i < size; ++i)
    {
        MPI_Bcast(l[i], size, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(u[i], size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    t1 = MPI_Wtime();

    // Compute inverse
    compute_inverse(a, l, u, a_inv, size, rank, num_procs);
    MPI_Barrier(MPI_COMM_WORLD);
 
    t2 = MPI_Wtime(); // Calculate elapsed time

    if (rank == 0)
        cout << "Compute inverse time: " << t2 - t1 << '\n'; 

    // Root process writes the output
    if (rank == 0)
    {
        ofstream l_out(get_output_filename(input_filename, "L"));
        ofstream u_out(get_output_filename(input_filename, "U"));
        ofstream a_inv_out(get_output_filename(input_filename, "A_inv"));

        print_matrix(l, size, l_out);
        print_matrix(u, size, u_out);
        print_matrix(a_inv, size, a_inv_out);
    }

    // Clean up
    for (int i = 0; i < size; ++i)
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

    end =  MPI_Wtime();
    if (rank == 0) {
        cout << "Total time: " << end - start << '\n';
    }

    MPI_Finalize();
    return 0;
}
