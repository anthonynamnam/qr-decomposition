import numpy as np
import copy

# Set up of Matrix Class
class Matrix:
    values = []
    determinant = None

    def __init__(self, n_row, n_col, all_zero=True, two_d_array=None):
        self.values = []
        self.n_row = n_row
        self.n_col = n_col
        if two_d_array is not None:
            assert type(two_d_array) == list and type(two_d_array[0]) == list
            n_row = max(len(two_d_array), n_row)
            n_col = max(max([len(row) for row in two_d_array]), n_col)
            for i in range(n_row):
                row_val = []
                for j in range(n_col):

                    try:
                        row_val.append(float(two_d_array[i][j]))
                    except:
                        row_val.append(0)
                self.values.append(row_val)
        elif all_zero:
            for i in range(n_col):
                row_val = []
                for j in range(n_row):
                    row_val.append(0.0)
                self.values.append(row_val)

        self.calculate_determinant(self.values)

    def set_diagonal(self, values=1):
        assert self.n_col == self.n_row, "set_diagonal only available to square matrix"
        for i in range(self.n_col):
            self.values[i][i] = values

    @staticmethod
    # Return the cofactor for calculating the determinant
    def sub_matrix(two_d_array, remove_row, remove_col):
        new_matrix = []
        for i in range(len(two_d_array)):
            if remove_row == i:
                continue
            new_row = []
            for j in range(len(two_d_array[i])):
                if remove_col == j:
                    continue
                new_row.append(two_d_array[i][j])
            new_matrix.append(new_row)
        return new_matrix

    def calculate_determinant(self, two_d_array):
        n_row = len(two_d_array)
        n_col = len(two_d_array[0])
        assert n_row == n_col, "calculate_determinant only available to square matrix"

        if n_row == 1:
            return np.sum(two_d_array)

        sum = 0
        for i in range(n_row):
            if i % 2 == 0:
                factor = 1
            else:
                factor = -1
            first_num_in_row = two_d_array[i][0]
            sub_matrix = self.sub_matrix(two_d_array, i, 0)
            sum += factor * first_num_in_row * self.calculate_determinant(sub_matrix)

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
        print(f"Matrix Shape: ({self.n_col},{self.n_row})")

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


# QR Decomposition by Gram-Schmidt Process
def QR_GS(M: Matrix) -> tuple():
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
        M_column_vector.append([row[i] for row in M.values])

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
def QR_householder(M: Matrix) -> tuple():

    Qs = []

    curr_M = M
    # ===========================================================
    for k in range(0, M.n_row):
        # Get the sub-matrix, nothing happened then k = 0
        sub_matrix_values = []
        for i in range(M.n_row):
            row = []
            for j in range(M.n_col):
                if i < k or j < k:
                    continue
                row.append(curr_M.values[i][j])
            if len(row) > 0:
                sub_matrix_values.append(row)

        # Get Sub-Matrix for next iteration
        sub_matrix = Matrix(
            n_row=M.n_row - k, n_col=M.n_col - k, two_d_array=sub_matrix_values
        )

        # ==================== Start Calculation ============================
        sign_a11 = 1 if sub_matrix.values[0][0] >= 0 else -1

        # Get the current column vector
        q = [row[0] for row in sub_matrix.values]
        e = [1 if j == 0 else 0 for j in range(sub_matrix.n_col)]

        # Normalise the column vector
        b = q - (-sign_a11) * sub_matrix.length(q) * np.array(e)
        q = b / sub_matrix.length(b)

        # Calculate matrix Q
        q_values = []
        for i in range(sub_matrix.n_row):
            row = []
            for j in range(sub_matrix.n_col):
                if i == j:
                    row.append(1 - 2 * q[i] * q[j])
                else:
                    row.append(-2 * q[i] * q[j])
            q_values.append(row)
        Q = Matrix(n_row=sub_matrix.n_row, n_col=sub_matrix.n_col, two_d_array=q_values)

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
                Q = Matrix(n_row=dim_first, n_col=dim_first, two_d_array=extended_Q_values)
        
        # ===========================================================
        # Save the Q matrix
        Qs.append(Q)
        curr_M = Q.multiply(M)
    
    # Calculate Q and R with all the {Q_1 ... Q_k}
    final_R = M
    final_Q = Matrix(n_row=M.n_row,n_col=M.n_col)
    final_Q.set_diagonal()
    
    for matrix in Qs:
        final_R = matrix.multiply(final_R)
        matrix.transpose()
        final_Q = final_Q.multiply(matrix)
    
    return final_Q, final_R