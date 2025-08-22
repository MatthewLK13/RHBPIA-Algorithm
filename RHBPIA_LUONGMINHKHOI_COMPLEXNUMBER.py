from sympy import symbols, diff, factorial, simplify, expand, sympify, lambdify, I
import matplotlib.pyplot as plt
import numpy as np

# Khai báo biến x là biến phức
x = symbols('x', complex=True)

def single_diff(expr, order, at_point):
    """
    Tính đạo hàm bậc `order` của biểu thức tại điểm `at_point`
    """
    return diff(expr, x, order).evalf(subs={x: at_point})

def rhbpia_simplified(Zn, alpha, R):
    """
    Zn     : Danh sách toàn bộ z_i (tọa độ các điểm)
    alpha  : Danh sách các danh sách bậc đạo hàm tương ứng tại mỗi z_i
    R      : Danh sách các danh sách giá trị tương ứng với từng alpha[i][k]
    """
    pi = sympify(0)
    N = []
    w = sympify(1)

    for i, zi in enumerate(Zn):
        print(f'voi i = {i} va zi = {zi}')
        alphas = alpha[i]
        values = R[i]
        d = w.evalf(subs={x: zi})
        print(f'd = {d} \n')
        for k in range(len(alphas)):
            print(f' voi k = {k}')
            a_k = alphas[k]
            print(' alphaik = ', a_k)
            if i == 0 and k == 0:
                i_k = -1 + 1 + k
            else:
                i_k = N[i - 1] + 1 + k
            print(' ik = ', i_k)
            if k == 0:
                g_k = simplify(w * (x - zi) ** a_k)
            else:
                g_k = simplify(g_k * (x - zi) ** (a_k - alphas[k - 1]))
            print(' gik = ', g_k)
            pi_val = single_diff(pi, a_k, zi)
            print(' daoham = ', pi_val)
            coeff = simplify((values[k] - pi_val) / (factorial(a_k) * d))
            print(' Aik = ', coeff)
            pi = simplify(pi + coeff * g_k)
            print(' piik = ', pi)
            if k == len(alphas) - 1:
                w = simplify(g_k * (x - zi))
        if i == 0:
            N.append(-1 + 1 + (len(alphas) - 1))
        else:
            N.append(N[i - 1] + 1 + (len(alphas) - 1))
        print('w = ', w)
        print('Ni = ', N[i], '\n')
    if i == len(Zn) - 1:
        print("\nInterpolated Polynomial π(x):", expand(pi))
    return expand(pi)

def plot_polynomial(poly_expr, xmin=-1, xmax=3):
    f = lambdify(x, poly_expr, modules='numpy')
    xs = np.linspace(xmin, xmax, 400)
    ys = f(xs)
    plt.figure(figsize=(8, 4))
    plt.plot(xs, ys.real, label="Re(π(x))", color="blue")
    plt.plot(xs, ys.imag, label="Im(π(x))", color="orange", linestyle='dashed')
    plt.axhline(0, color='black', lw=0.5)
    plt.axvline(0, color='black', lw=0.5)
    plt.grid(True)
    plt.xlabel('x')
    plt.ylabel('π(x)')
    plt.title('Graph of Interpolated Polynomial (Real and Imag)')
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    j = I
    # Dữ liệu ví dụ: nội suy phức
    # Zn = [1 + j, 2 - j]
    # alpha = [
    #     [0],
    #     [0, 1]
    # ]
    # R = [
    #     [1 + j],
    #     [1, j]
    # ]

    Zn = [1 + j, 2 - j, 3 + 2j]
    alpha = [
        [2],
        [0, 2],
        [1]
    ]
    R = [
        [2+j],
        [-1+2j,5j],
        [1+3j]
    ]
    
    interpolated_poly = rhbpia_simplified(Zn, alpha, R)
    plot_polynomial(interpolated_poly)
