import random
import numpy as np
import time
import math

execution_times = []
solutions_values = []



# Constants for the lower and upper bounds of the functions
RASTRIGIN_LOWER = -5.12
RASTRIGIN_UPPER = 5.12

MICHALEWICZ_LOWER = 0
MICHALEWICZ_UPPER = np.pi

SCHWEFEL_LOWER = -500
SCHWEFEL_UPPER = 500

DE_JONG_1_LOWER = -5.12
DE_JONG_1_UPPER = 5.12


def rastrigin(x):
    # Rastrigin function
    sum_val = sum([(xi ** 2 - 10 * np.cos(2 * np.pi * xi)) for xi in x])
    return 10 * len(x) + sum_val


def michalewicz(x):
    d = len(x)
    _sum = 0

    for ii in range(1, d + 1):
        xi = x[ii - 1]
        new = np.sin(xi) * (np.sin(ii * xi ** 2 / np.pi)) ** (20)
        _sum += new

    y = -_sum
    return y


def de_jong_1(x):
    return np.sum(np.square(x))


def schwefel(x):
    return -np.sum(x * np.sin(np.sqrt(np.abs(x))))


def evaluate(candidate_coordinates, l, function, lower, upper):
    num_elements = len(candidate_coordinates) // l
    evaluated_values = []

    for i in range(num_elements):
        start_idx = i * l
        end_idx = (i + 1) * l
        binary_substring = candidate_coordinates[start_idx:end_idx]
        candidate_as_int = int(''.join(map(str, binary_substring)), 2)
        value = round(lower + (candidate_as_int / ((2 ** l) - 1)) * (upper - lower), 5)
        evaluated_values.append(value)

    # print(f"representation in number {evaluated_values} .")
    results = round(function(evaluated_values), 5)

    return results


# Worst Improvement
def worst_improvement(actual_no, neighbors, no_of_bits, l, function, lower, upper):
    best_value = - np.inf

    for i in range(0, len(neighbors) // no_of_bits):
        neighbor_start = i * no_of_bits
        neighbor_end = (i + 1) * no_of_bits
        current_neighbor = neighbors[neighbor_start:neighbor_end]
        current_value = evaluate(current_neighbor, l, function, lower, upper)

        if best_value < current_value < actual_no:
            best_value = current_value
            best_neighbor = current_neighbor

    if best_value == -np.inf:
        best_neighbor = neighbors[:no_of_bits]
        best_value = evaluate(best_neighbor, l, function, lower, upper)

    return best_value, best_neighbor


# First Improvement
# def first_improvement(actual_no, neighbors, no_of_bits, l, function, lower, upper):
#
#     random_index = random.randint(0, len(neighbors) - 1)
#     print (random_index)
#
#     neighbor_index = random_index // no_of_bits
#
#     neighbor_start = neighbor_index * no_of_bits
#     neighbor_end = (neighbor_index + 1) * no_of_bits
#
#     current_neighbor = neighbors[neighbor_start:neighbor_end]
#     current_value = evaluate(current_neighbor, l, function, lower, upper)
#
#     return current_value, current_neighbor


def first_improvement(actual_no, neighbors, no_of_bits, l, function, lower, upper):
    num_neighbors = len(neighbors) // no_of_bits  # Calculate the number of neighbors
    neighbor_indices = list(range(num_neighbors))  # Create a list of neighbor indices

    random.shuffle(neighbor_indices)  # Shuffle the order in which neighbors are evaluated

    for neighbor_index in neighbor_indices:
        neighbor_start = neighbor_index * no_of_bits
        neighbor_end = (neighbor_index + 1) * no_of_bits

        current_neighbor = neighbors[neighbor_start:neighbor_end]
        current_value = evaluate(current_neighbor, l, function, lower, upper)

        if current_value < actual_no:  # Check if the current value is better
            return current_value, current_neighbor

    return current_value, current_neighbor


def best_improvement(actual_no, neighbors, no_of_bits, l, function, lower, upper):
    best_neighbor = neighbors[:no_of_bits]
    best_value = evaluate(best_neighbor, l, function, lower, upper)

    for i in range(1, len(neighbors) // no_of_bits):
        neighbor_start = i * no_of_bits
        neighbor_end = (i + 1) * no_of_bits
        current_neighbor = neighbors[neighbor_start:neighbor_end]
        current_value = evaluate(current_neighbor, l, function, lower, upper)

        if current_value < best_value:
            best_value = current_value
            best_neighbor = current_neighbor

    if  best_value > actual_no:

        num_neighbors = len(neighbors) // no_of_bits  # Calculate the number of neighbors

        neighbor_indices = random.randint(0, num_neighbors - 1)

        neighbor_start = neighbor_indices * no_of_bits
        neighbor_end = (neighbor_indices + 1) * no_of_bits

        current_neighbor = neighbors[neighbor_start:neighbor_end]
        current_value = evaluate(current_neighbor, l, function, lower, upper)

        return current_value, current_neighbor

    return best_value, best_neighbor


def generate_neighbors(candidate_coordinates, l):
    neighbors = []

    for i in range(len(candidate_coordinates) // l):
        for j in range(l):
            neighbor = candidate_coordinates.copy()
            neighbor[i * l + j] = 1 - neighbor[i * l + j]
            neighbors.extend(neighbor)

    return neighbors


def hill_climbing(improve, dimensions, function, lower, upper):
    start = time.time()

    t = 0
    best = np.inf

    l = math.ceil(math.log2((upper - lower) * epsilon))

    # print(f"has exactly {l} bits.")

    no_of_bits = dimensions * l

    # print(f"has in total {no_of_bits} bits.")

    while t != MAX:
        local = False

        # Generate a random binary array of length no_of_bits
        candidate_coordinates = [random.choice([0, 1]) for _ in range(no_of_bits)]
        # print(candidate_coordinates)

        actual_no = evaluate(candidate_coordinates, l, function, lower, upper)

        while not local:

            neighbors = generate_neighbors(candidate_coordinates, l)

            chosen_neighbor, byte_repr_of_chosen_neighbor = improve(actual_no, neighbors, no_of_bits, l, function, lower, upper);

            if chosen_neighbor < actual_no:
                actual_no = chosen_neighbor
                candidate_coordinates = byte_repr_of_chosen_neighbor
            else:
                local = True

        t = t + 1

        if actual_no < best:
            best = actual_no

        # print(actual_no)

    # https://www.geeksforgeeks.org/how-to-check-the-execution-time-of-python-script/
    end = time.time()
    execution_time = (end - start)
    execution_times.append(round(execution_time, 3))
    solutions_values.append(round(best, 5))

    return best


def simulated_annealing_iterative(improve, dimensions, function, lower, upper):
    start = time.time()

    t = 0
    best = np.inf
    l = math.ceil(math.log2((upper - lower) * epsilon))
    no_of_bits = dimensions * l

    while t != MAX:

        T = 100
        j = 0
        # Generate a random binary array of length no_of_bits
        candidate_coordinates = [random.choice([0, 1]) for _ in range(no_of_bits)]

        actual_no = evaluate(candidate_coordinates, l, function, lower, upper)

        while T > 0.00000001 and j < 40:
            local = False
            i=0
            while not local:
                neighbors = generate_neighbors(candidate_coordinates, l)
                chosen_neighbor, byte_repr_of_chosen_neighbor = improve(actual_no, neighbors, no_of_bits, l, function, lower, upper);
                # print("chosen : ", chosen_neighbor, " si actual : ", actual_no)
                if chosen_neighbor < actual_no:
                    actual_no = chosen_neighbor
                    candidate_coordinates = byte_repr_of_chosen_neighbor

                elif random.choice([0, 1]) < math.exp(-abs(chosen_neighbor - actual_no) / T):
                    j=0
                    actual_no = chosen_neighbor
                    candidate_coordinates = byte_repr_of_chosen_neighbor
                else:
                    local = True

            j = j + 1
            T = T * 0.99

        if actual_no < best:
            best = actual_no

        #print(actual_no)

        t = t + 1

    # https://www.geeksforgeeks.org/how-to-check-the-execution-time-of-python-script/
    end = time.time()
    execution_time = (end - start)
    execution_times.append(round(execution_time, 3))
    solutions_values.append(round(best, 5))

    return best


def print_time_statistics():
    avg_time = sum(execution_times) / len(execution_times)
    print("Average execution time:", avg_time)


def print_function_values_statistics():
    avg_solution = sum(solutions_values) / len(solutions_values)
    best_solution = min(solutions_values)
    print("Average solution:", avg_solution)
    print("Best solution:", best_solution)


if __name__ == "__main__":

    improvement_methods = [first_improvement, best_improvement, worst_improvement]
    dimensions = [5,10,30]

    epsilon = 10 ** 3
    MAX = 10000
    RUNS = 30

    # 5000 for 5, 1000 for 10, 100 for 30

    # MICHALEWICZ
    # for dimension in dimensions:
    #     for improvement in improvement_methods:
    #         execution_times = []
    #         solutions_values = []
    #         for i in range(RUNS):
    #             simulated_annealing_iterative(improvement, dimension, michalewicz, MICHALEWICZ_LOWER, MICHALEWICZ_UPPER)
    #         print(f"Simulated Annealing : MICHALEWICZ in {dimension} dimensions with {improvement.__name__} in {MAX} iterations.")
    #         print_time_statistics();
    #         print_function_values_statistics();



    # # MICHALEWICZ
    # for dimension in dimensions:
    #     for improvement in improvement_methods:
    #         execution_times = []
    #         solutions_values = []
    #         for i in range(RUNS):
    #             hill_climbing(improvement, dimension, michalewicz, MICHALEWICZ_LOWER, MICHALEWICZ_UPPER)
    #         print(f"Hill Climbing : MICHALEWICZ in {dimension} dimensions with {improvement.__name__} in {MAX} iterations.")
    #         print_time_statistics();
    #         print_function_values_statistics();
    #
    # # RASTRIGIN
    # for dimension in dimensions:
    #     for improvement in improvement_methods:
    #         execution_times = []
    #         solutions_values = []
    #         for i in range(RUNS):
    #             hill_climbing(improvement, dimension, rastrigin, RASTRIGIN_LOWER, RASTRIGIN_UPPER)
    #         print(f"Hill Climbing : RASTRIGIN in {dimension} dimensions with {improvement.__name__} in {MAX} iterations.")
    #         print_time_statistics();
    #         print_function_values_statistics();
    #
    # # DE_JONG_1
    # for dimension in dimensions:
    #     for improvement in improvement_methods:
    #         execution_times = []
    #         solutions_values = []
    #         for i in range(RUNS):
    #             hill_climbing(improvement, dimension, de_jong_1, DE_JONG_1_LOWER, DE_JONG_1_UPPER)
    #         print(f"Hill Climbing : DE_JONG_1 in {dimension} dimensions with {improvement.__name__} in {MAX} iterations.")
    #         print_time_statistics();
    #         print_function_values_statistics();

    # SCHWEFEL
    # for dimension in dimensions:
    #     for improvement in improvement_methods:
    #         execution_times = []
    #         solutions_values = []
    #         for i in range(RUNS):
    #             hill_climbing(improvement, dimension, schwefel, SCHWEFEL_LOWER, SCHWEFEL_UPPER)
    #         print(f"Hill Climbing : SCHWEFEL in {dimension} dimensions with {improvement.__name__} in {MAX} iterations.")
    #         print_time_statistics();
    #         print_function_values_statistics();
