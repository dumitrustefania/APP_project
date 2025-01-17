import random

# Fill the matrix with random values
def random_fill(size):
    # Fill the matrix with random values
    matrix = [[0 for _ in range(size)] for _ in range(size)]
    for i in range(size):
        for j in range(size):
            matrix[i][j] = random.randint(-5000, 5000)  # Generates numbers from -5000 to 5000

    # Ensure the matrix is diagonal dominant to guarantee invertible-ness
    # diagCount will help keep track of which column the diagonal is in
    diag_count = 0
    sum_values = 0
    for i in range(size):
        for j in range(size):
            # Sum all column values
            sum_values += abs(matrix[i][j])
        # Remove the diagonal value from the sum
        sum_values -= abs(matrix[i][diag_count])
        # Add a random value to the sum and place in diagonal position
        matrix[i][diag_count] = sum_values + (random.randint(1, 5))
        diag_count += 1
        sum_values = 0

    return matrix

if __name__ == "__main__":
    # Random seed
    random.seed(1)

    # Size of matrix
    size = int(input())
    print(size)

    matrix = random_fill(size)

    for i in range(size):
        for j in range(size):
            print(f"{matrix[i][j]} ", end="")
        print()
