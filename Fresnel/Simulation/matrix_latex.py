import sympy as sp

def multiply_and_export_latex():
    # 1. Define the algebraic symbols you plan to use
    a, b, c, d, x, y = sp.symbols('a b c d x y')

    # 2. Define your first matrix (Matrix A)
    # This example is a 2x2 matrix, but you can add more rows/columns
    matrix_A = sp.Matrix([
        [a, b],
        [c, d]
    ])

    # 3. Define your second matrix (Matrix B)
    # Ensure inner dimensions match for multiplication (e.g., 2x2 * 2x2)
    matrix_B = sp.Matrix([
        [x, y],
        [x**2, y**2]  # You can include powers, fractions, etc.
    ])

    # 4. Perform the matrix multiplication
    result_matrix = matrix_A * matrix_B

    # Optional: Simplify the algebraic terms in the result
    simplified_result = sp.simplify(result_matrix)

    # 5. Convert the result to a LaTeX formatted string
    latex_string = sp.latex(simplified_result)

    # Print the result so you can copy it
    print("--- Copy the text below into your LaTeX editor ---")
    print(latex_string)

if __name__ == "__main__":
    multiply_and_export_latex()