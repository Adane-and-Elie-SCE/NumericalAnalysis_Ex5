#Adane Adgo 315721969
#Elie Bracha 204795900

# github: https://github.com/Adane-and-Elie-SCE/NumericalAnalysis_Ex5.git


def determinant_calc(mat):
    if len(mat) == 2:
        ans = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]
        return ans

    minor = [[0 for i in range(len(mat) - 1)] for j in range(len(mat) - 1)]
    determinant = 0

    for k in range(len(mat)):
        i, j = 0, 0
        while i < len(mat):
            if i != k:
                minor[j] = mat[i][1:]
                j += 1
            i += 1
        determinant += ((-1) ** k) * mat[k][0] * determinant_calc(minor)
    return determinant


def elementary_delete(mat, tag):
    i, j = tag
    ans = [[int(i == j) for j in range(len(mat))] for i in range(len(mat))]
    ans[i][j] = -1 * (mat[i][j] / mat[j][j])
    return ans


def elementary_switch_rows(mat, tag):
    ans = [[int(i == j) for j in range(len(mat))] for i in range(len(mat))]
    ans[tag[0]], ans[tag[1]] = ans[tag[1]], ans[tag[0]]
    return ans


def elementary_mul_row(mat, index, scalar):
    e = [[int(i == j) for j in range(len(mat))] for i in range(len(mat))]
    mat[index][index] = scalar
    return mat


def matrix_pivoting(mat, index):
    size = len(mat)
    max_value = mat[index][index]
    max_index = index

    for i in range(index + 1, size):
        if mat[i][index] > max_value:
            max_value = mat[i][index]
            max_index = i

    return matrix_mul(elementary_switch_rows(mat, [index, max_index]), mat)


def fix_approximation(mat):
    epsilon = 2 ** (-26)
    for i in range(len(mat)):
        for j in range(len(mat)):
            if abs(mat[i][j]) < epsilon:
                mat[i][j] = 0
    return mat


def matrix_mul(a, b):
    size = len(a)
    ans = [[0 for i in range(size)] for j in range(size)]
    for i in range(size):
        for j in range(size):
            for k in range(size):
                ans[i][j] += a[i][k] * b[k][j]
    return fix_approximation(ans)


def matrix_vector_mul(mat, b):
    size = len(mat)
    ans = [0 for i in range(size)]

    for i in range(size):
        for j in range(size):
            ans[i] += mat[i][j] * b[j]
    return ans


def print_matrix(mat):
    for row in mat:
        print(row)
    print()


def print_vector(v):
    size = len(v)
    for i in range(size):
        print('[' + str(v[i]) + ']')


def inverse_by_gauss(mat):
    if not determinant_calc(mat):
        print('no inverse')

    size = len(mat)
    temp = mat
    ans = [[int(i == j) for j in range(size)] for i in range(size)]

    # below diagonal
    for j in range(0, size - 1):

        for i in range(j + 1, size):
            e = elementary_delete(temp, [i, j])
            temp = matrix_mul(e, temp)
            ans = matrix_mul(e, ans)

    # above diagonal
    for j in range(size - 1, -1, -1):

        for i in range(j - 1, -1, -1):
            e = elementary_delete(temp, [i, j])
            temp = matrix_mul(e, temp)
            ans = matrix_mul(e, ans)

    # final step
    for i in range(size):
        for j in range(size):
            ans[i][j] /= temp[i][i]

    return ans


def lu_decomposition(mat):
    size = len(mat)
    u = mat
    l = [[int(i == j) for j in range(size)] for i in range(size)]

    # below diagonal
    for j in range(0, size):
        for i in range(j + 1, size):
            e = elementary_delete(u, [i, j])
            u = matrix_mul(e, u)
            e[i][j] = -1 * e[i][j]
            l = matrix_mul(l, e)

    print("Triangle Matrix L:")
    print_matrix(l)
    print("Triangle Matrix U:")
    print_matrix(u)


# ------------------------------------------------------------------


list1 = [[0, 0], [1, 0.8415], [2, 0.9093], [3, 0.1411], [4, -0.7568], [5, -0.9589], [6, -0.2794]]

list2 = [[1, 0.8415], [2, 0.9093], [3, 0.1411]]

list3 = [[1, 1], [2, 0], [4, 1.5]]

#list4 = [[1, 0.7651], [1.3, 0.62], [1.6, 0.4554], [1.9, 0.2818], [2.2, 0.1103]]

x = 2.5


def linear(list, x):
    for i in range(len(list)):
        if list[i][0] < x < list[i + 1][0]:
            result = (list[i][1] - list[i + 1][1]) / (list[i][0] - list[i + 1][0]) * x + (
                    (list[i + 1][1] * list[i][0]) - (list[i][1] * list[i + 1][0])) / (list[i][0] - list[i + 1][0])
            return result


def polynomial(list, x):
    sum = 0
    for i in range(len(list)):
        sum += list[i] * pow(x, i)
    return sum


def Lagrange(list, x):
    L = []
    fact = 1
    for i in range(len(list)):
        for j in range(len(list)):
            if i != j:
                fact *= (x - list[j][0]) / (list[i][0] - list[j][0])
        L.append(fact)
        fact = 1
    sum = 0
    for i in range(len(list)):
        sum += L[i] * list[i][1]

    return sum


def Neville(list, x):
    n = len(list)
    p = n * [0]
    for k in range(n):
        for i in range(n - k):
            if k == 0:
                p[i] = list[i][1]
            else:
                p[i] = ((x - list[i + k][0]) * p[i] + (list[i][0] - x) * p[i + 1]) / (list[i][0] - list[i + k][0])
    return p[0]


def Amatrix(list):
    matrix = []
    row = []
    for i in range(len(list)):
        for j in range(len(list)):
            row.append(pow(list[i][0], j))
        matrix.append(row)
        row = []
    return matrix


def bvector(list):
    vector = []
    for i in range(len(list)):
        vector.append(list[i][1])
    return vector

def Driver(list, x):

    a = input("X = 2,5 \nWhat Method do you want ? \n1- Linear \n2- Polynomial \n3- Lagrange \nElse- Neville\n")


    if a == "1":
        print("Linear: y = " + str(linear(list, x)))
    elif a == "2":
        print("Polynomial: y = " + str(polynomial(matrix_vector_mul(inverse_by_gauss(Amatrix(list)), bvector(list)), x)))
    elif a == "3":
        print("Lagrange: y = " + str(Lagrange(list, x)))
    else:
        print("Neville: y = " + str(Neville(list, x)))


Driver(list1, x)
