import example

# example.cgs()

# example.mgs()

# example.householder()

# example.givens_rotation()

# example.least_squares_problem()


from Matrix import Matrix
import math

A = Matrix(two_d_array=[[2, 3, 6], [-1, -3, 4], [9, -3, 5]])
A.print_values(matrix_name="A")
G = Matrix(n_col=3, n_row=3)
G.set_diagonal()
i_col = 1
j_col = 2
rotate_theta = 90  # in degree, anti-clockwise

# Swap it, make sure i is smaller than j
if i_col > j_col:
    i_col, j_col = j_col, i_col

# Construct Rotation matrix
if i_col != j_col:
    min_dim = min(G.n_col, G.n_row)
    if i_col <= min_dim and i_col <= min_dim:
        theta = math.radians(rotate_theta)  # in radians
        G.values[i_col - 1][i_col - 1] = math.cos(theta)
        G.values[j_col - 1][j_col - 1] = math.cos(theta)

        G.values[i_col - 1][j_col - 1] = -1 * math.sin(theta)
        G.values[j_col - 1][i_col - 1] = 1 * math.sin(theta)

G.print_values(matrix_name="G")
new_A = G.multiply(A)
new_A.print_values(matrix_name="new_A")
