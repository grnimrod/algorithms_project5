import numpy as np

"""
We need to calculate the Q matrix at every iteration of the NJ algorithm
We need to have Ri and Rj calculated in order to obtain the Q matrix
Ri and Rj: the sum of all the distances of the remaining "labels"/"sequences" leading into Ri/Rj
    Qij = (n - 2) * Dij - Ri - Rj
The non-identical pair with the lowest net divergence are joined to form a new node u
Then we need to calculate the distances from u to the rest of the nodes

Steps:
1. Having D, calculate Q matrix
2. Pick the pair with lowest score
3. Calculate distances from original nodes to newly created node
    d(f, u) = 1/2 (Dfg) + 1/(2 * (n - 2)) * [Rf - Rg]
    d(g, u) = d(f, g) - d(f, u)
4. Calculate new distance matrix using
    Duk = 1/2 [Dfk + Dgk - Dfg]
5. Repeat steps 1-4 until we have three remaining nodes in D
"""

def initiate_dist_matrix(file):
    raw = open(file, "r").read()
    lines_sep = raw.split('\n')

    dim = int(lines_sep[0])
    dist_matrix = np.zeros((dim, dim))

    lines_sep = lines_sep[1:]
    
    labels = []

    for index, line in enumerate(lines_sep):
        parsed_line = line.split()
        label = parsed_line[0]
        labels.append(label)
        dist_matrix[index, :] = parsed_line[1:]

    return dist_matrix, labels


def calculate_Q(dist_matrix):
    dim = dist_matrix.shape[0]
    Q = np.zeros((dim, dim))

    for i in range(dim):
        for j in range(dim):
            if i == j:
                Q[i][j] = 0
            else:
                r_i, r_j = 0, 0
                for k in range(dim):
                    r_i += dist_matrix[i][k]
                    r_j += dist_matrix[j][k]
                Q[i][j] = (dim - 2) * dist_matrix[i][j] - r_i - r_j
    
    return Q


def find_lowest_pair(Q_matrix):
    dim = Q_matrix.shape[0]
    min_value = float('inf')

    for i in range(dim):
        for j in range(dim):
            if Q_matrix[i][j] < min_value:
                min_value = Q_matrix[i][j]
                min_index = (i, j)
    
    return min_index


def calc_dist_from_orig_to_joined_node(index, dist_matrix):
    # We use these distances as branch length for the tree that we're building

    dim = dist_matrix.shape[0]
    i, j = index
    r_i, r_j = 0, 0
    for k in range(dim):
        r_i += dist_matrix[i][k]
        r_j += dist_matrix[j][k]
    
    d_ki = (1/2) * (dist_matrix[i][j] + r_i - r_j)
    d_kj = (1/2) * (dist_matrix[i][j] - r_i + r_j)

    return d_ki, d_kj


def update_dist_matrix(index, dist_matrix):
    dim = dist_matrix.shape[0]
    new_dist_matrix = np.zeros((dim - 1, dim - 1)) # -1 because we are deleting 2 rows and columns and adding 1

    # Deleting the rows and columns corresponding to the nodes that we are joining together
    f, g = index
    to_delete = [f, g]

    mask = np.ones(dist_matrix.shape[0], dtype = bool)
    mask[to_delete] = False
    temp_matrix = dist_matrix[mask, :][:, mask]

    """
    Copying remaining values from the original distance matrix to the updated one,
    Calculating distances from new node to rest of the nodes
    """
    dim_temp = temp_matrix.shape[0]

    for n in range(dim_temp):
        for m in range(dim_temp):
            new_dist_matrix[n][m] = temp_matrix[n][m]
    
    index_in_new_dist_matrix = 0

    for k in range(dim):
        if k == f or k == g:
            continue
        else:
            new_dist_matrix[dim - 2][index_in_new_dist_matrix] = (1/2) * (dist_matrix[f][k] + dist_matrix[k][g] - dist_matrix[f][g])
            new_dist_matrix[index_in_new_dist_matrix][dim - 2] = (1/2) * (dist_matrix[f][k] + dist_matrix[k][g] - dist_matrix[f][g])
            index_in_new_dist_matrix += 1

    return new_dist_matrix


D, labels = initiate_dist_matrix("./data/example_slide4.phy")
Q = calculate_Q(D)
min_index = find_lowest_pair(Q)
dist_to_joined = calc_dist_from_orig_to_joined_node(min_index, D)
nd = update_dist_matrix(min_index, D)
print(nd)