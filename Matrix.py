import numpy as np
import copy


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

    # def get_numpy_determinant(self):
    #     print(f"numpy determinant:{np.linalg.det(np.array(self.values))}")

    def transpose(self):
        new_values = []
        for i in range(self.n_col):
            new_row = []
            for j in range(self.n_row):
                new_row.append(self.values[j][i])
            new_values.append(new_row)
        self.values = new_values

    def inverse(self):
        pass

    def print_shape(self):
        print(f"Matrix Shape: ({self.n_col},{self.n_row})")

    def print_determinant(self):
        print(f"Determinant: {self.determinant}")

    def print_values(self, matrix_name=""):
        if len(matrix_name) > 0:
            print(f"{matrix_name}:")
        for row in self.values:
            for row_val in row:
                print(f"{row_val}\t", end="")
            print("")
        print("")


def QR_MGS(M: Matrix, print_step=False) -> tuple():

    M.print_values()
    return
