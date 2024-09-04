import random
from itertools import combinations, product

# Parametri
n = 15
n_minus_k = 9
w_r = 5
w_c = 3
seed = 78  # Indeks

def generate_random_permutation(seed, size):
    random.seed(seed)
    permutation = list(range(size))
    random.shuffle(permutation)
    return permutation

def create_matrix_H(n, n_minus_k, w_r, w_c, seed):
    H = [[0] * n for _ in range(n_minus_k)]
    
    # Postavljanje početnih vrednosti u matrici H
    for i in range(w_c):
        for j in range(w_r):
            H[i][(i * w_r + j) % n] = 1

    permutations = [generate_random_permutation(seed + i, n) for i in range(1, w_c)]
    
    for i in range(w_c, n_minus_k):
        perm = permutations[(i - w_c) % (w_c - 1)]
        for j in range(n):
            if H[i % w_c][j] == 1:
                H[i][j] = 0
                H[i][perm[j]] = 1

    return H

def compute_syndrome(H, error_vector):
    syndrome = [0] * len(H)
    for i, row in enumerate(H):
        syndrome[i] = sum(row[j] * error_vector[j] for j in range(len(error_vector))) % 2
    return tuple(syndrome)

def list_all_error_patterns(max_weight, n):
    patterns = []
    for weight in range(1, max_weight + 1):
        for positions in combinations(range(n), weight):
            pattern = [0] * n
            for pos in positions:
                pattern[pos] = 1
            patterns.append(pattern)
    return patterns

def create_syndrome_and_correction_table(H):
    n = len(H[0])
    syndrome_table = {}
    max_weight = 2 * (n - len(H)) + 1
    error_patterns = list_all_error_patterns(max_weight, n)
    
    for error_pattern in error_patterns:
        syndrome = compute_syndrome(H, error_pattern)
        weight = sum(error_pattern)
        
        if syndrome not in syndrome_table or weight < syndrome_table[syndrome][0]:
            syndrome_table[syndrome] = (weight, error_pattern)

    return syndrome_table

def display_syndrome_table(syndrome_table):
    print("Tabela sindroma i korektora:")
    for syndrome, (weight, error) in sorted(syndrome_table.items()):
        syndrome_str = ''.join(map(str, syndrome))
        error_str = ''.join(map(str, error))
        print(f'Sindrom: {syndrome_str} -> Korektor: {error_str} (težina: {weight})')

def hamming_distance(x, y):
    return sum(e1 != e2 for e1, e2 in zip(x, y))

def parity_check(H, x):
    return [(sum(H[i][j] * x[j] for j in range(len(x))) % 2) for i in range(len(H))]

def update_check_to_variable_messages(H, messages, variable_values):
    m = len(H)
    n = len(H[0])
    new_messages = [row[:] for row in messages]
    
    for i in range(m):
        for j in range(n):
            if H[i][j] == 1:
                neighbors = [k for k in range(n) if k != j and H[i][k] == 1]
                product = 1
                for neighbor in neighbors:
                    product *= (2 * messages[i][neighbor] - 1) / 2
                new_messages[i][j] = (1 if product > 0 else -1) if abs(product) < 1 else 1
    
    return new_messages

def update_variable_nodes(H, messages, variable_values, threshold_0=0.5, threshold_1=0.5):
    m = len(H)
    n = len(H[0])
    new_variable_values = [0] * n
    
    for j in range(n):
        num_positive = sum(1 for i in range(m) if messages[i][j] > 0)
        num_negative = m - num_positive
        num_connected = sum(1 for i in range(m) if H[i][j] == 1)
        
        if num_negative >= threshold_0 * num_connected:
            new_variable_values[j] = 0
        elif num_positive >= threshold_1 * num_connected:
            new_variable_values[j] = 1
        else:
            new_variable_values[j] = variable_values[j]
    
    return new_variable_values

def gallager_b_decoder(H, received_vector, threshold_0=0.5, threshold_1=0.5, max_iterations=100):
    m = len(H)
    n = len(H[0])
    
    messages = [[0] * n for _ in range(m)]
    current_estimate = received_vector[:]
    
    for _ in range(max_iterations):
        messages = update_check_to_variable_messages(H, messages, current_estimate)
        new_estimate = update_variable_nodes(H, messages, current_estimate, threshold_0, threshold_1)
        
        if sum(parity_check(H, new_estimate)) == 0:
            return new_estimate
        
        current_estimate = new_estimate
    
    raise ValueError("Presegnut je maksimalan broj iteracija.")

def display_decoding_results(received_vector, decoded_vector):
    print("Primljeni vektor =", received_vector)
    print("Dekodirani vektor =", decoded_vector)
    print()

def generate_code_words(H):
    from itertools import product
    n = len(H[0])
    code_words = []
    for vector in product([0, 1], repeat=n):
        if sum(parity_check(H, vector)) == 0:
            code_words.append(list(vector))
    return code_words

def calculate_code_distance(code_words):
    min_distance = float('inf')
    num_words = len(code_words)
    
    for i in range(num_words):
        for j in range(i + 1, num_words):
            dist = hamming_distance(code_words[i], code_words[j])
            if dist < min_distance:
                min_distance = dist
    
    return min_distance

# Kreirajte matricu H
H = create_matrix_H(n, n_minus_k, w_r, w_c, seed)
print("Matrica H:")
for row in H:
    print(row)

# Kreirajte tabelu sindroma
syndrome_table = create_syndrome_and_correction_table(H)
display_syndrome_table(syndrome_table)

# Generišite sve kodne reči
code_words = generate_code_words(H)
code_distance = calculate_code_distance(code_words)
print(f"Minimalno kodno rastojanje: {code_distance}")

received_vector = [random.randint(0, 1) for _ in range(n)]

# Dekodiranje koristeći Gallager B
decoded_vector = gallager_b_decoder(H, received_vector)
display_decoding_results(received_vector, decoded_vector)

# Pronađi n-torku greške koja ne može biti ispravljena
def find_minimum_error_pattern(H, threshold_0=0.5, threshold_1=0.5):

    min_errors = float('inf')
    best_pattern = None
    
    # Tražimo n-torku greške sa najmanje jedinica
    for weight in range(1, len(H[0]) + 1):
        for pattern in combinations(range(len(H[0])), weight):
            error_vector = [0] * len(H[0])
            for pos in pattern:
                error_vector[pos] = 1
            
            # Dekodiraj grešku
            decoded = gallager_b_decoder(H, error_vector, threshold_0, threshold_1)
            if sum(decoded) != len(H[0]): 
                if weight < min_errors:
                    min_errors = weight
                    best_pattern = pattern
    
    return min_errors, best_pattern

min_errors, best_pattern = find_minimum_error_pattern(H)
print(f"N-torka greške sa najmanje jedinica koja ne može biti ispravljena: {min_errors}")
print(f"Greška obrazac: {best_pattern}")
