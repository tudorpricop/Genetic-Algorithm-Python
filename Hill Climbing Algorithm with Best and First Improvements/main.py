import random
import numpy as np
import time
import math

execution_times = []
solutions_values = []
avg_vec_4100 = []


# Constants for the lower and upper bounds of the function

FUNCTION_LOWER = 0
FUNCTION_UPPER = 31

veciniVizitati = 0
vecVizitati4100 = 0
iteratii = 0

multimap = {}


def function(x):
    return x ** 3 - 60 * x ** 2 + 900 * x + 100


def evaluate(candidate_coordinates, l, function, lower, upper):
    num_elements = len(candidate_coordinates) // l
    evaluated_values = []

    for i in range(num_elements):
        start_idx = i * l
        end_idx = (i + 1) * l
        binary_substring = candidate_coordinates[start_idx:end_idx]
        candidate_as_int = int(''.join(map(str, binary_substring)), 2)
        value = lower + (candidate_as_int / ((2 ** l) - 1)) * (upper - lower)
        # evaluated_values.append(value)

    # print(f"representation in number {evaluated_values} .")
    results = function(value)

    return results


def first_improvement(actual_no, neighbors, no_of_bits, l, function, lower, upper):
    num_neighbors = len(neighbors) // no_of_bits  # Calculate the number of neighbors
    neighbor_indices = list(range(num_neighbors))  # Create a list of neighbor indices

    random.shuffle(neighbor_indices)  # Shuffle the order in which neighbors are evaluated

    for neighbor_index in neighbor_indices:

        global veciniVizitati  # Specificam ca vrem sa accesam variabila globala
        veciniVizitati += 1

        neighbor_start = neighbor_index * no_of_bits
        neighbor_end = (neighbor_index + 1) * no_of_bits

        current_neighbor = neighbors[neighbor_start:neighbor_end]
        current_value = evaluate(current_neighbor, l, function, lower, upper)

        if current_value > actual_no:  # Check if the current value is better (modificat)
            return current_value, current_neighbor

    return current_value, current_neighbor


def best_improvement(actual_no, neighbors, no_of_bits, l, function, lower, upper):
    best_neighbor = neighbors[:no_of_bits]
    best_value = evaluate(best_neighbor, l, function, lower, upper)

    for i in range(1, len(neighbors) // no_of_bits):

        global veciniVizitati  # Specificam ca vrem sa accesam variabila globala
        veciniVizitati += 1

        neighbor_start = i * no_of_bits
        neighbor_end = (i + 1) * no_of_bits
        current_neighbor = neighbors[neighbor_start:neighbor_end]
        current_value = evaluate(current_neighbor, l, function, lower, upper)

        if current_value > best_value:  # (modificat)
            best_value = current_value
            best_neighbor = current_neighbor

    return best_value, best_neighbor


def generate_neighbors(candidate_coordinates, l):
    neighbors = []

    for i in range(len(candidate_coordinates) // l):
        for j in range(l):
            neighbor = candidate_coordinates.copy()
            neighbor[i * l + j] = 1 - neighbor[i * l + j]
            neighbors.extend(neighbor)

    return neighbors


def num_to_bin_list(number, length=5):
    binary_str = bin(number)[2:].zfill(length)
    return [int(bit) for bit in binary_str]


def hill_climbing(i, improve, function, lower, upper):
    start = time.time()
    best = -np.inf  # modif

    l = math.ceil(math.log2((upper - lower)))

    no_of_bits = l

    local = False

    candidate_coordinates = num_to_bin_list(i)

    # print (i, function(i))
    # print(candidate_coordinates)

    actual_no = function(i)

    while not local:

        neighbors = generate_neighbors(candidate_coordinates, l)

        chosen_neighbor, byte_repr_of_chosen_neighbor = improve(actual_no, neighbors, no_of_bits, l, function, lower,
                                                                upper)

        if chosen_neighbor > actual_no:  # (modificat)
            global vecVizitati4100
            vecVizitati4100 += 1
            actual_no = chosen_neighbor
            candidate_coordinates = byte_repr_of_chosen_neighbor
            # print(actual_no, candidate_coordinates)
        else:
            local = True

    if actual_no > best:  # (modificat)
        best = actual_no

    if best in multimap:
        multimap[best].append(i)
    else:
        multimap[best] = [i]

    print(actual_no)

    # https://www.geeksforgeeks.org/how-to-check-the-execution-time-of-python-script/

    return actual_no


def print_time_statistics():
    avg_time = sum(execution_times) / len(execution_times)
    print("Average execution time:", avg_time)


def print_function_values_statistics():
    avg_solution = sum(solutions_values) / len(solutions_values)
    best_solution = min(solutions_values)
    # print("Average solution:", avg_solution)
    print("Best solution:", best_solution)


if __name__ == "__main__":

    improvement_methods = []

    RUNS = 1

    for j in range(RUNS):
        veciniVizitati = 0
        multimap = {}
        # functia data
        for i in range(32):
            value = hill_climbing(i, best_improvement, function, FUNCTION_LOWER, FUNCTION_UPPER)
            if value == 3988:
                iteratii += 1
                avg_vec_4100.append(vecVizitati4100)
            vecVizitati4100 = 0
            # print(f"Functia de maxim pornim din punctul {i}: with {improvement.__name__} ")
            # print_function_values_statistics();
            # print()
        solutions_values.append(veciniVizitati)

    avg_solution = sum(solutions_values) / len(solutions_values)

    print(f"In medie, in {RUNS} RUNS, am scanat {avg_solution} vecini cu hill_climbing first improvement")
    print(
        f"In medie trecem prin {sum(avg_vec_4100)/ len(avg_vec_4100)} vecini pana sa ajungem cu hill_climbing la punctul de maxim local (4100)")
    print ("iteratii", iteratii)
    # print(f"Am vizitat {veciniVizitati} vecini")

    for key, values in multimap.items():
        values_as_str = ", ".join(str(value) for value in values)
        print(
            f"Bazinul de atractie este format din punctele {values_as_str} pentru care Hill Climbing a gasit punctul de maxim local f(x) = {key}")
