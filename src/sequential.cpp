#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <chrono>

using namespace std;

void print_matrix(float **matrix, int size, ofstream &output_file);
void LU_decomposition(float **a, float **l, float **u, int size);
void initialize_matrices(float **a, float **l, float **u, int size);
void forward_substitution(float **L, float *b, float *y, int size);
void backward_substitution(float **U, float *y, float *x, int size);
void compute_inverse(float **a, float **l, float **u, float **a_inv, int size);
string get_output_filename(const string &input_filename, const string &matrix_type);

// Forward substitution function (L * y = b)
void forward_substitution(float **L, float *b, float *y, int size)
{
    for (int i = 0; i < size; ++i)
    {
        y[i] = b[i]; // Initialize y[i] with b[i]
        for (int j = 0; j < i; ++j)
        {
            y[i] -= L[i][j] * y[j]; // Subtract the previous values in the row
        }
        y[i] /= L[i][i]; // Divide by L[i][i] to get the solution for y[i]
    }
}

// Backward substitution function (U * x = y)
void backward_substitution(float **U, float *y, float *x, int size)
{
    for (int i = size - 1; i >= 0; --i)
    {
        x[i] = y[i]; // Initialize x[i] with y[i]
        for (int j = i + 1; j < size; ++j)
        {
            x[i] -= U[i][j] * x[j]; // Subtract the previous values in the row
        }
        x[i] /= U[i][i]; // Divide by U[i][i] to get the solution for x[i]
    }
}

// Function to compute the inverse of matrix 'a' using LU decomposition
void compute_inverse(float **a, float **l, float **u, float **a_inv, int size)
{
    float *y = new float[size]; // Temporary array for forward substitution
    float *x = new float[size]; // Temporary array for backward substitution

    for (int i = 0; i < size; ++i)
    {
        float *e = new float[size](); // Create the unit vector (e) with a 1 at position i
        e[i] = 1;
        forward_substitution(l, e, y, size);  // Solve L * y = e
        backward_substitution(u, y, x, size); // Solve U * x = y
        for (int j = 0; j < size; ++j)
        {
            a_inv[j][i] = x[j]; // Store the solution x (the i-th column of the inverse)
        }
        delete[] e;
    }
    delete[] y;
    delete[] x;
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

// LU decomposition function (L * U = A)
void LU_decomposition(float **a, float **l, float **u, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            // Compute the lower triangular matrix L
            if (j < i)
            {
                l[j][i] = 0; // Set the lower triangular matrix L
            }
            else
            {
                l[j][i] = a[j][i]; // Copy elements from A to L
                for (int k = 0; k < i; k++)
                {
                    l[j][i] -= l[j][k] * u[k][i]; // Subtract the previously computed elements in the row
                }
            }

            // Compute the upper triangular matrix U
            if (j < i)
            {
                u[i][j] = 0; // Set the upper triangular matrix U
            }
            else if (j == i)
            {
                u[i][j] = 1; // Diagonal of U is set to 1 (since L is lower triangular)
            }
            else
            {
                u[i][j] = a[i][j] / l[i][i]; // Compute the upper triangular matrix U
                for (int k = 0; k < i; k++)
                {
                    u[i][j] -= (l[i][k] * u[k][j]) / l[i][i]; // Subtract the previously computed elements
                }
            }
        }
    }
}

// Function to generate the output filename for each matrix (L, U, A_inv)
string get_output_filename(const string &input_filename, const string &matrix_type)
{
    string dir = "outputs/sequential/" + input_filename.substr(7, input_filename.find(".") - 7);
    filesystem::create_directories(dir); // Create the directory if it doesn't exist
    return dir + "/" + matrix_type + ".out";
}

int main(int argc, char *argv[])
{
    auto start = chrono::high_resolution_clock::now();

    if (argc < 2)
    {
        cerr << "Usage: " << argv[0] << " <input_file>" << '\n';
        return 1;
    }

    string input_filename = argv[1];
    ifstream input_file(input_filename);

    if (!input_file.is_open())
    {
        cerr << "Error opening input file." << '\n';
        return 1;
    }

    int size;
    input_file >> size; // Read the size of the matrix

    float **a = new float *[size];     // Matrix A
    float **l = new float *[size];     // Lower triangular matrix L
    float **u = new float *[size];     // Upper triangular matrix U
    float **a_inv = new float *[size]; // Inverse matrix of A
    for (int i = 0; i < size; ++i)
    {
        a[i] = new float[size];
        l[i] = new float[size];
        u[i] = new float[size];
        a_inv[i] = new float[size];
    }

    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            input_file >> a[i][j]; // Read the elements of matrix A from the input file
        }
    }
    input_file.close();

    // Measure the runtimes
    double runtime;
    runtime = clock() / (double)CLOCKS_PER_SEC;

    LU_decomposition(a, l, u, size); // Perform LU decomposition

    runtime = (clock() / (double)CLOCKS_PER_SEC) - runtime; 
    cout << "LU time:" << runtime << '\n';

    runtime = clock() / (double)CLOCKS_PER_SEC;

    compute_inverse(a, l, u, a_inv, size); // Compute the inverse of matrix A

     runtime = (clock() / (double)CLOCKS_PER_SEC) - runtime; 

         cout << "Compute inverse time: " << runtime << '\n';


    ofstream l_out(get_output_filename(input_filename, "L"));
    ofstream u_out(get_output_filename(input_filename, "U"));
    ofstream a_inv_out(get_output_filename(input_filename, "A_inv"));
    print_matrix(l, size, l_out);
    print_matrix(u, size, u_out);
    print_matrix(a_inv, size, a_inv_out);

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

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    cout << "Total time: " << elapsed.count() << " seconds" << '\n';
    return 0;
}
