import math
import numpy as np


# Set up of Matrix Class
class Matrix:
    values = []
    determinant = None

    def __init__(self, n_row=1, n_col=1, all_zero=True, two_d_array=None):
        self.values = []
        self.n_row = n_row
        self.n_col = n_col
        if two_d_array is not None:
            assert type(two_d_array) == list and type(two_d_array[0]) == list
            n_row = max(len(two_d_array), n_row)
            n_col = max(max([len(row) for row in two_d_array]), n_col)
            self.n_row = n_row
            self.n_col = n_col
            for i in range(n_row):
                row_val = []
                for j in range(n_col):

                    try:
                        row_val.append(float(two_d_array[i][j]))
                    except:
                        row_val.append(0)
                self.values.append(row_val)
        elif all_zero:
            for i in range(n_row):
                row_val = []
                for j in range(n_col):
                    row_val.append(0.0)
                self.values.append(row_val)

        self.calculate_determinant(self.values)

    def deep_copy(self):
        new_M = Matrix(self.n_row, self.n_col, two_d_array=self.values)
        return new_M

    def set_diagonal(self, values=1):
        assert self.n_col == self.n_row, "set_diagonal only available to square matrix"
        for i in range(self.n_col):
            self.values[i][i] = values

    @staticmethod
    # Return the cofactor for calculating the determinant
    def sub_matrix_value(two_d_array, remove_row, remove_col):
        new_matrix = []
        if type(remove_row) == str:
            remove_row = [int(remove_row)]
        elif type(remove_row) == int:
            remove_row = [remove_row]

        if type(remove_col) == str:
            remove_col = [int(remove_col)]
        elif type(remove_col) == int:
            remove_col = [remove_col]

        if type(remove_row) == list:
            for i in range(len(two_d_array)):
                if i in remove_row:
                    continue
                new_row = []
                for j in range(len(two_d_array[i])):
                    if j in remove_col:
                        continue
                    new_row.append(two_d_array[i][j])
                new_matrix.append(new_row)

        return new_matrix

    def sub_matrix(self, remove_row, remove_col):
        sub_matrix_value = self.sub_matrix_value(self.values, remove_row, remove_col)
        return Matrix(
            self.n_row - len(remove_row),
            self.n_col - len(remove_col),
            two_d_array=sub_matrix_value,
        )

    def calculate_determinant(self, two_d_array):
        n_row = len(two_d_array)
        n_col = len(two_d_array[0])
        if n_row != n_col:
            # print(f">> calculate_determinant only available to square matrix")
            self.determinant = None
            return

        if n_row == 1:
            return np.sum(two_d_array)

        sum = 0
        for i in range(n_row):
            if i % 2 == 0:
                factor = 1
            else:
                factor = -1
            first_num_in_row = two_d_array[i][0]
            sub_matrix_value = self.sub_matrix_value(two_d_array, i, 0)
            sum += (
                factor * first_num_in_row * self.calculate_determinant(sub_matrix_value)
            )

        self.determinant = sum
        return sum

    def transpose(self):
        new_values = []
        for i in range(self.n_col):
            new_row = []
            for j in range(self.n_row):
                new_row.append(self.values[j][i])
            new_values.append(new_row)
        self.values = new_values

    def print_shape(self):
        print(f"Matrix Shape: ({self.n_row},{self.n_col})")

    def print_determinant(self):
        print(f"Determinant: {self.determinant}")

    def print_values(self, matrix_name="", num_digits=4):
        if len(matrix_name) > 0:
            print(f"{matrix_name}:")
        for row in self.values:
            for row_val in row:
                print(f"{round(row_val,num_digits)}\t", end="")
            print("")
        print("")

    # Get column vectors
    def get_column_vectors(self, column_index):
        assert column_index < self.n_col
        return [row[column_index] for row in self.values]

    # Get the projection of vector on another vector
    def proj_on(self, vector, on_vector):
        assert len(vector) == len(on_vector), f"Vectors must be same length"
        dot_prod = self.dot_product(vector, on_vector)
        on_vector_length = self.length(on_vector) ** 2
        # print(f"Dot: {dot_prod}")
        # print(f"on_vector length: {on_vector_length}")
        res = (dot_prod / on_vector_length) * np.array(on_vector)
        # print(f"Project of {vector} on {on_vector} RES: {res}")
        return res

    @staticmethod
    # Get the length of a vector
    def length(vector):
        return np.sqrt(np.sum([i**2 for i in vector]))

    @staticmethod
    # Dot Product function
    def dot_product(vector_1, vector_2):
        return np.sum([i * j for i, j in zip(vector_1, vector_2)])

    def multiply(self, new_M):
        assert type(new_M) == Matrix
        assert self.n_col == new_M.n_row
        res_M = Matrix(self.n_row, new_M.n_col)
        for i in range(self.n_row):
            for j in range(new_M.n_col):
                left_vector = self.values[i]
                right_vector = [row[j] for row in new_M.values]
                res_M.values[i][j] = self.dot_product(left_vector, right_vector)
        return res_M

    def isUpperTriangular(self):
        for i in range(self.n_row):
            for j in range(self.n_col):
                if i > j:
                    if abs(self.values[i][j] - 0) >= 1e-12:
                        return False
                else:
                    continue
        return True

    def isLowerTriangular(self):
        for i in range(self.n_row):
            for j in range(self.n_col):
                if i < j:
                    if abs(self.values[i][j] - 0) >= 1e-12:
                        return False
                else:
                    continue
        return True


def backward_substitution(A: Matrix, b: Matrix):
    assert b.n_col == 1
    assert A.n_col == b.n_row
    assert A.isUpperTriangular()
    result = []
    for i in range(A.n_row):
        result.append(0)
    for j in range(A.n_col):
        result[2 - j] = (
            b.values[2 - j][0] - A.dot_product(A.values[2 - j], result)
        ) / (A.values[2 - j][2 - j])

    result = [[k] for k in result]
    x = Matrix(n_row=A.n_row, n_col=b.n_col, two_d_array=result)
    return x


# QR Decomposition by Gram-Schmidt Process
def QR_CGS(M: Matrix) -> tuple():
    Q_vectors = []
    for i in range(M.n_col):
        # Get the column vector
        curr_column = [row[i] for row in M.values]

        # Removing Projection of current vector on calculated orthogonal vector
        q = curr_column
        if i > 0:
            for k in range(i):
                q -= M.proj_on(curr_column, Q_vectors[k])

        # Normalise the column vector
        q = q / M.length(q)

        # Store the vector in the row, need to do transpose when completed all iteration
        Q_vectors.append(list(q))

    # Create a new matrix of the Q_vectors
    Q = Matrix(n_col=M.n_col, n_row=M.n_row, two_d_array=Q_vectors)

    # As we perform the iteration column-by-column, we need to do transpose to get the correct Q matrix
    Q.transpose()

    # --- Here we have the correct Q Matrix ---

    # Calculate R by tranpose of Q times A
    Q.transpose()
    R = Q.multiply(M)

    # Convert back to the correct Q
    Q.transpose()

    return (Q, R)


# QR Decomposition by Modified Gram-Schmidt Process
def QR_MGS(M: Matrix) -> tuple():

    # Get all column vectors in M
    M_column_vector = []
    for i in range(M.n_col):
        M_column_vector.append(M.get_column_vectors(i))

    # Start iteration
    Q_vectors = []
    for i in range(M.n_col):
        # Get the current column vector
        q = M_column_vector[i]

        # Normalise the column vector
        q = q / M.length(q)

        # Store the vector in the row, need to do transpose when completed all iteration
        Q_vectors.append(list(q))

        # Removing Projection of current vector on calculated orthogonal vector
        for j in range(i, M.n_col):
            M_column_vector[j] -= M.proj_on(M_column_vector[j], q)

    # Create a new matrix of the Q_vectors
    Q = Matrix(n_col=M.n_col, n_row=M.n_row, two_d_array=Q_vectors)

    # As we perform the iteration column-by-column, we need to do transpose to get the correct Q matrix
    Q.transpose()

    # --- Here we have the correct Q Matrix ---

    # Calculate R by tranpose of Q times A
    Q.transpose()
    R = Q.multiply(M)

    # Convert back to the correct Q
    Q.transpose()

    return (Q, R)


# QR Decomposition by Householder Transformation
def QR_Householder(M: Matrix) -> tuple():

    Qs = []

    curr_M = M
    # ===========================================================
    for k in range(0, M.n_row):
        # Get the sub-matrix, nothing happened then k = 0
        minor_matrix_values = []
        for i in range(M.n_row):
            row = []
            for j in range(M.n_col):
                if i < k or j < k:
                    continue
                row.append(curr_M.values[i][j])
            if len(row) > 0:
                minor_matrix_values.append(row)

        # Get Sub-Matrix for next iteration
        minor_matrix = Matrix(
            n_row=M.n_row - k, n_col=M.n_col - k, two_d_array=minor_matrix_values
        )

        # ==================== Start Calculation ============================
        sign_a11 = 1 if minor_matrix.values[0][0] >= 0 else -1

        # Get the current column vector
        q = minor_matrix.get_column_vectors(0)
        e = [1 if j == 0 else 0 for j in range(minor_matrix.n_col)]

        # Normalise the column vector
        b = q - (-sign_a11) * minor_matrix.length(q) * np.array(e)
        q = b / minor_matrix.length(b)

        # Calculate matrix Q
        q_values = []
        for i in range(minor_matrix.n_row):
            row = []
            for j in range(minor_matrix.n_col):
                if i == j:
                    row.append(1 - 2 * q[i] * q[j])
                else:
                    row.append(-2 * q[i] * q[j])
            q_values.append(row)
        Q = Matrix(
            n_row=minor_matrix.n_row, n_col=minor_matrix.n_col, two_d_array=q_values
        )

        # ===========================================================
        # Check if need extend, nothing to do if k = 0
        if Q.n_row < M.n_row and Q.n_col < M.n_col:
            # Extend the Matrix to original size
            dim_this = Q.n_row
            dim_first = M.n_row
            dim_diff = dim_first - dim_this
            if dim_diff > 0:
                extended_Q_values = []
                for i in range(dim_first):
                    row = []
                    for j in range(dim_first):
                        # Add the extended part
                        if i < dim_diff or j < dim_diff:
                            if i == j:
                                row.append(1)
                            else:
                                row.append(0)
                        # Add back the result of the Q
                        else:
                            row.append(Q.values[i - dim_diff][j - dim_diff])
                    extended_Q_values.append(row)
                Q = Matrix(
                    n_row=dim_first, n_col=dim_first, two_d_array=extended_Q_values
                )

        # ===========================================================
        # Save the Q matrix
        Qs.append(Q)
        curr_M = Q.multiply(M)

    # Calculate Q and R with all the {Q_1 ... Q_k}
    final_R = M
    final_Q = Matrix(n_row=M.n_row, n_col=M.n_col)
    final_Q.set_diagonal()

    for matrix in Qs:
        final_R = matrix.multiply(final_R)
        matrix.transpose()
        final_Q = final_Q.multiply(matrix)

    return final_Q, final_R


def QR_Givens_Rotation(M: Matrix) -> tuple():
    G_s = []
    R = M.deep_copy()

    for j in range(R.n_col):
        for i in range(R.n_row - 1, j, -1):
            if i > j:

                # Construct Givens Rotation Matrix
                col_to_remove = list(
                    filter(lambda s: s != i and s != j, [i for i in range(R.n_row)])
                )
                sub = R.sub_matrix(col_to_remove, col_to_remove)
                col_vector = sub.get_column_vectors(0)

                # Calculate the angle to rotate
                theta = math.atan2(
                    -col_vector[1], col_vector[0]
                )  # y-coordinate first, then x-coordinate

                # Construct Givens Rotation Matrix
                G_mat = Matrix(n_row=R.n_row, n_col=R.n_col)
                G_mat.set_diagonal()
                for row in [i, j]:
                    for col in [i, j]:
                        if row == col:
                            G_mat.values[row][col] = math.cos(theta)
                        else:
                            G_mat.values[row][col] = (
                                1 if row > col else -1
                            ) * math.sin(theta)

                # Save the Givens Rotation Matrix
                G_s.append(G_mat)

                # Calculate the R Matrix for next iteration
                R = G_mat.multiply(R)

    # Calculate the Q matrix by multiplying all the transpose of Givens Matrix
    Q = Matrix(M.n_row, M.n_col)
    Q.set_diagonal(1)
    for G in G_s:
        G.transpose()
        Q = Q.multiply(G)
        G.transpose()

    return (Q, R)
