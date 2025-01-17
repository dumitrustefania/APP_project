import numpy as np

def load_matrix(filename):
    """
    Load a matrix from a file. Assumes the first line contains the size.
    """
    matrix = []
    with open(filename, 'r') as f:
        size = int(f.readline())
        for line in f:
            if line.strip():  # Skip empty lines
                row = [float(num) for num in line.split()]
                matrix.append(row)
    # Convert to a NumPy array
    try:
        return np.array(matrix)
    except ValueError:
        print(f"Error: Inconsistent matrix dimensions in file {filename}.")
        return None

def verify_inverse(a, a_inv):
    """
    Verify if the product of a matrix and its inverse equals the identity matrix.
    """
    if a is None or a_inv is None:
        return False
    identity = np.eye(a.shape[0])
    product = np.dot(a, a_inv)
    return np.allclose(product, identity, atol=1e-3)

def verify_lu_decomposition(a, l, u):
    """
    Verify if the product of L and U equals the original matrix A.
    """
    if a is None or l is None or u is None:
        return False
    reconstructed_a = np.dot(l, u)
    return np.allclose(reconstructed_a, a, atol=1e-2)

if __name__ == "__main__":
    # Prompt the user for the implementation type (e.g. sequential, openmp, mpi, hybrid)
    alg_type = input("Enter the implementation type (sequential, openmp, mpi, hybrid): ")

    # Prompt the user for the matrix size (e.g., small, medium, large)
    matrix_size = input("Enter the matrix size (small, medium, large): ")

    # Construct the file paths
    a_filename = f"../inputs/{matrix_size}.in"
    a_inv_filename = f"../outputs/{alg_type}/{matrix_size}/A_inv.out"
    l_filename = f"../outputs/{alg_type}/{matrix_size}/L.out"
    u_filename = f"../outputs/{alg_type}/{matrix_size}/U.out"

    try:
        # Load matrices
        a = load_matrix(a_filename)
        a_inv = load_matrix(a_inv_filename)
        l = load_matrix(l_filename)
        u = load_matrix(u_filename)

        # Verify the inverse
        if a is None or a_inv is None:
            print("Matrix loading failed due to input file issues.")
        elif verify_inverse(a, a_inv):
            print("Verification successful: A * A^-1 = I")
        else:
            print("Verification failed: A * A^-1 != I")

        # Verify the LU decomposition
        if a is None or l is None or u is None:
            print("Matrix loading failed due to input file issues.")
        elif verify_lu_decomposition(a, l, u):
            print("Verification successful: L * U = A")
        else:
            print("Verification failed: L * U != A")
    except FileNotFoundError as e:
        print(f"Error: {e}")
