import numpy as np

def read_rvecs_from_chk(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    rvecs_floats = []
    found_rvecs = False
    for i, line in enumerate(lines):
        if line.strip().startswith('rvecs'):
            found_rvecs = True
            # Start reading lines after 'rvecs' to get 9 float values
            for next_line in lines[i+1:]:
                try:
                    numbers = list(map(float, next_line.strip().split()))
                    rvecs_floats.extend(numbers)
                    if len(rvecs_floats) >= 9:
                        break
                except ValueError:
                    continue
            break

    if not found_rvecs:
        raise ValueError("rvecs not found in the file.")
    if len(rvecs_floats) < 9:
        raise ValueError(f"Could not read 9 float values for rvecs. Got {len(rvecs_floats)}.")

    rvecs = np.array(rvecs_floats[:9]).reshape((3, 3))
    return rvecs

def compute_cell_parameters(rvecs_bohr, nr_ucells):
    bohr_to_angstrom = 0.529177210903
    rvecs = rvecs_bohr * bohr_to_angstrom

    a_vec, b_vec, c_vec = rvecs
    a = np.linalg.norm(a_vec)
    b = np.linalg.norm(b_vec)
    c = np.linalg.norm(c_vec)

    alpha = np.degrees(np.arccos(np.dot(b_vec, c_vec) / (b * c)))
    beta  = np.degrees(np.arccos(np.dot(a_vec, c_vec) / (a * c)))
    gamma = np.degrees(np.arccos(np.dot(a_vec, b_vec) / (a * b)))

    volume = np.abs(np.dot(a_vec, np.cross(b_vec, c_vec)))

    return a/nr_ucells, b/nr_ucells, c/nr_ucells, alpha, beta, gamma, volume/nr_ucells

if __name__ == "__main__":
    import sys
    filename = sys.argv[1] if len(sys.argv) > 1 else 'init.chk'
    nr_ucells = 2
    rvecs_bohr = read_rvecs_from_chk(filename)
    a, b, c, alpha, beta, gamma, volume = compute_cell_parameters(rvecs_bohr, nr_ucells)

    print(f"Cell lengths (Å): a = {a:.4f}, b = {b:.4f}, c = {c:.4f}")
    print(f"Cell angles (°): α = {alpha:.2f}, β = {beta:.2f}, γ = {gamma:.2f}")
    print(f"Cell volume: {volume:.3f} Å³")
