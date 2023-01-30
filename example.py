from Matrix import (
    Matrix,
    QR_CGS,
    QR_MGS,
    QR_Householder,
    QR_Givens_Rotation,
    backward_substitution,
)


# ----- Example of Gram Schmdit Process -----
def cgs():

    A = Matrix(
        n_row=3,
        n_col=3,
        two_d_array=[[2, -1, -2], [-4, 6, 3], [-4, -2, 8]],
    )

    print("----- Example of Gram-Schmidt Process -----")
    A.print_values(matrix_name="A")
    Q, R = QR_CGS(A)

    Q.print_values(matrix_name="Q", num_digits=4)

    R.print_values(matrix_name="R", num_digits=4)

    print("----- Validate Result -----")
    new_A = Q.multiply(R)
    new_A.print_values(matrix_name="new_A", num_digits=4)


# ----- Example of Modified Gram Schmdit Process -----
def mgs():

    A = Matrix(
        n_row=3,
        n_col=3,
        two_d_array=[[2, -1, -2], [-4, 6, 3], [-4, -2, 8]],
    )

    print("----- Example of Modified Gram-Schmidt Process -----")
    A.print_values(matrix_name="A")
    Q, R = QR_MGS(A)

    Q.print_values(matrix_name="Q", num_digits=4)

    R.print_values(matrix_name="R", num_digits=4)

    print("----- Validate Result -----")
    new_A = Q.multiply(R)
    new_A.print_values(matrix_name="new_A", num_digits=4)


def householder():

    A = Matrix(
        n_row=3,
        n_col=3,
        two_d_array=[[2, -1, -2], [-4, 6, 3], [-4, -2, 8]],
    )

    print("----- Example of Householder Transformation -----")
    A.print_values(matrix_name="A")
    Q, R = QR_Householder(A)

    Q.print_values(matrix_name="Q", num_digits=4)

    R.print_values(matrix_name="R", num_digits=4)

    print("----- Validate Result -----")
    new_A = Q.multiply(R)
    new_A.print_values(matrix_name="new_A", num_digits=4)


def givens_rotation():

    A = Matrix(
        n_row=3,
        n_col=3,
        two_d_array=[[2, -1, -2], [-4, 6, 3], [-4, -2, 8]],
    )

    print("----- Givens Rotations -----")
    A.print_values(matrix_name="A")
    Q, R = QR_Givens_Rotation(A)

    Q.print_values(matrix_name="Q", num_digits=4)

    R.print_values(matrix_name="R", num_digits=4)

    print("----- Validate Result -----")
    new_A = Q.multiply(R)
    new_A.print_values(matrix_name="new_A", num_digits=4)


def least_squares_problem():

    print(f"\n ----- Example of Solving Least Square Problem by QR Decomposition -----")
    A = Matrix(
        n_row=3,
        n_col=3,
        two_d_array=[[2, -1, -2], [-4, 6, 3], [-4, -2, 8]],
    )
    A.print_values(matrix_name="A")
    b = Matrix(
        n_row=3,
        n_col=1,
        two_d_array=[[9], [-2], [5]],
    )
    b.print_values(matrix_name="b")

    print(f"\n----- QR Decompostition Result -----")
    Q, R = QR_Givens_Rotation(A)
    Q.print_values(matrix_name="Q")
    R.print_values(matrix_name="R")

    print(f"\n----- Compute y by multiply transpose of Q to b -----")
    Q.transpose()
    y = Q.multiply(b)
    y.print_values(matrix_name="y")

    print(f"\n----- Compute x by backward substitution -----")
    R.print_values(matrix_name="R")
    x = backward_substitution(R, y)
    x.print_values(matrix_name="x")

    print(f"\n----- Compute Ax and validate if it is equal to b -----")
    A.print_values(matrix_name="A")
    x.print_values(matrix_name="x")
    validate = A.multiply(x)
    validate.print_values(matrix_name="validate")
