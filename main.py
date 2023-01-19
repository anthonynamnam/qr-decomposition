from Matrix import Matrix

# ----- Example of QR Decomposition -----
# Create the matrix
A = Matrix(
    n_row=3,
    n_col=3,
    two_d_array=[[2, -1, -2], [-4, 6, 3], [-4, -2, 8]],
)

# ----- Example of Modified Gram Schmdit -----
# Matrix.QR_MGS(A, method="MGS", print_step=False)
