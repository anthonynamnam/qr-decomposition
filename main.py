from Matrix import Matrix, QR_GS, QR_MGS, QR_Householder, QR_GivensRotations

# ----- Example of QR Decomposition -----
# Create the matrix
A = Matrix(
    n_row=3,
    n_col=3,
    two_d_array=[[2, -1, -2], [-4, 6, 3], [-4, -2, 8]],
)

# main.py
# ----- Example of Gram Schmdit Process -----
print("----- Gram-Schmidt Process -----")
A.print_values(matrix_name="A")
Q, R = QR_GS(A)

Q.print_values(matrix_name="Q", num_digits=4)

R.print_values(matrix_name="R", num_digits=4)

print("----- Validate Result -----")
new_A = Q.multiply(R)
new_A.print_values(matrix_name="new_A", num_digits=4)


# print("===============================================================")

# main.py
# ----- Example of Modified Gram Schmdit Process -----
print("----- Modified Gram-Schmidt Process -----")
A.print_values(matrix_name="A")
Q, R = QR_MGS(A)

Q.print_values(matrix_name="Q", num_digits=4)

R.print_values(matrix_name="R", num_digits=4)

print("----- Validate Result -----")
new_A = Q.multiply(R)
new_A.print_values(matrix_name="new_A", num_digits=4)

# print("===============================================================")

# main.py
# ----- Example of Householder Transformation -----
print("----- Householder Transformation -----")
A.print_values(matrix_name="A")
Q, R = QR_Householder(A)

Q.print_values(matrix_name="Q", num_digits=4)

R.print_values(matrix_name="R", num_digits=4)

print("----- Validate Result -----")
new_A = Q.multiply(R)
new_A.print_values(matrix_name="new_A", num_digits=4)

# print("===============================================================")

# main.py
# ----- Example of Givens Rotations -----
print("----- Givens Rotations -----")

A.print_values(matrix_name="A")
Q, R = QR_GivensRotations(A)

Q.print_values(matrix_name="Q", num_digits=4)

R.print_values(matrix_name="R", num_digits=4)

print("----- Validate Result -----")
new_A = Q.multiply(R)
new_A.print_values(matrix_name="new_A", num_digits=4)
