import random
import numpy as np
import time
import math
import matplotlib.pyplot as plt

generations_list = []
function_values_list = []
execution_times = []
solutions_values = []

# Constants for the lower and upper bounds of the functions
RASTRIGIN_LOWER = -5.12
RASTRIGIN_UPPER = 5.12
RASTRIGIN_DOMAIN = RASTRIGIN_UPPER- RASTRIGIN_LOWER

MICHALEWICZ_LOWER = 0
MICHALEWICZ_UPPER = np.pi
MICHALEWICZ_DOMAIN = MICHALEWICZ_UPPER - MICHALEWICZ_LOWER

SCHWEFEL_LOWER = -500
SCHWEFEL_UPPER = 500
SCHWEFEL_DOMAIN = SCHWEFEL_UPPER - SCHWEFEL_LOWER

DE_JONG_1_LOWER = -5.12
DE_JONG_1_UPPER = 5.12
DE_JONG_1_DOMAIN = DE_JONG_1_UPPER - DE_JONG_1_LOWER

mutation_probability = 0.001

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


def gray_to_binary(gray_code):
    binary_code = [gray_code[0]]
    for i in range(1, len(gray_code)):
        binary_code.append(gray_code[i] ^ binary_code[i - 1])
    return binary_code


def evaluate(cromozome_coordinates, l, function, lower, upper):
    num_elements = len(cromozome_coordinates) // l
    evaluated_values = []

    for i in range(num_elements):
        start_idx = i * l
        end_idx = (i + 1) * l

        binary_substring = cromozome_coordinates[start_idx:end_idx]

        # Convert Gray code to binary
        binary_substring = gray_to_binary(binary_substring)

        binary_substring_str = [str(bit) for bit in binary_substring]
        candidate_as_int = int(''.join(binary_substring_str), 2)

        value = lower + (candidate_as_int / ((2 ** l) - 1)) * (upper - lower)
        evaluated_values.append(value)

    results = function(evaluated_values)

    return results


def selection(pt, E):
    T = 0
    p = [0] * (pop_size-k )
    q = [0] * (pop_size-k+ 1)

    emax = -np.inf  # Inițializați cu o valoare negativă infinită
    emin = np.inf  # Inițializați cu o valoare pozitivă infinită

    for value in E:
        if value > emax:
            emax = value
        if value < emin:
            emin = value

    Fitness = []

    pressure = 1

    for value in E:
        Fitness.append(((emax - value) / (emax - emin + epsilon) + 1) ** pressure)
        T += Fitness[-1]  # Adăugați ultimul element adăugat la Fitness la suma

    for i in range(pop_size-k):
        p[i] = Fitness[i] / T

    q[0] = 0

    for i in range(pop_size-k):
        q[i + 1] = q[i] + p[i]

    q[pop_size-k] = 1

    next_pop = []

    for i in range(pop_size-k):
        r = random.random()
        for j in range(pop_size-k):
            if q[j] < r <= q[j + 1]:
                next_pop.append(pt[j])
                break

    return next_pop


def mutate(pt, no_of_bits):
    for i in range(len(pt)):
        for j in range(len(pt[i])):
            r = random.random()
            if r <= mutation_probability:
                pt[i][j] = 1 - pt[i][j]


def cx(c1, c2):
    a = c1.copy()
    b = c2.copy()
    pos1 = 1+random.randint(0, len(c1) - 3)
    pos2 = 1+random.randint(0, len(c1) - 2)
    while pos1 == pos2:
        pos2 = 1+random.randint(0, len(c1) - 2)

    if pos1 > pos2:
        pos1, pos2 = pos2, pos1

    a[pos1:pos2], b[pos1:pos2] = b[pos1:pos2], a[pos1:pos2]
    return a, b


def crossover(pt, l, function, lower, upper):
    p = [(i, random.random()) for i in range(pop_size)]
    p.sort(key=lambda x: x[1])

    i = 0

    while i < len(p):
        if i + 1 == len(p) or p[i + 1][1] >= crossover_probability:
            break

        x = cx(pt[p[i][0]], pt[p[i + 1][0]])

        Max_Of_4 = []

        Max_Of_4.append((evaluate(x[0], l, function, lower, upper),x[0]))
        Max_Of_4.append((evaluate(x[1], l, function, lower, upper),x[1]))
        Max_Of_4.append((evaluate(pt[p[i][0]], l, function, lower, upper),pt[p[i][0]]))
        Max_Of_4.append((evaluate(pt[p[i + 1][0]], l, function, lower, upper),pt[p[i + 1][0]]))

        sorted_max = sorted(Max_Of_4, key=lambda x: x[0], reverse=False)

        max1, max2 = sorted_max[0], sorted_max[1]
        alt_numar_max1, alt_numar_max2 = max1[1], max2[1]

        pt[p[i][0]] = alt_numar_max1
        pt[p[i + 1][0]] = alt_numar_max2

        i += 2


def elitism(Pt, E, k):
    p = [(i, E[i]) for i in range(len(E))]
    p.sort(key=lambda x: x[1])
    el = [Pt[p[i][0]] for i in range(k)]
    return el

def ag(dimensions, function, function_domain, lower, upper):
    start = time.time()
    epsilon2 = 10 ** 5

    l = math.ceil(math.log2((upper - lower) * epsilon2))
    # print(f"has exactly {l} bits.")

    no_of_bits = dimensions * l
    # print(f"has in total {no_of_bits} bits.")

    t = 0
    pt = []
    minimum = np.inf
    E = []

    global mutation_probability
    mini2 = np.inf
    minimum2 = np.inf
    best_in_t_multiple_of_30 = 0

    for curent_pop_size in range(pop_size):
        # Generate a random binary array of length no_of_bits
        cromozome_coordinates = [random.choice([0, 1]) for _ in range(no_of_bits)]
        # print(cromozome_coordinates)

        pt.append(cromozome_coordinates)

        E.append(evaluate(cromozome_coordinates, l, function, lower, upper))

    while t < generations:
        t += 1

        if t % 20 == 0 and t> 1600:
            minimum2 = np.inf
            for e in E:
                if e < minimum2:
                    minimum2 = e

            if abs(best_in_t_multiple_of_30 - minimum2) < function_domain/102.4:
                #print("Crestem probabilitatea de mutatie", best_in_t_multiple_of_30, minimum2 )
                mutation_probability = 0.002
            else:
                mutation_probability = 0.001
            best_in_t_multiple_of_30 = minimum2

        if t % 20 == 0 and t <= 1500:
            minimum2 = np.inf
            for e in E:
                if e < minimum2:
                    minimum2 = e
            # print(E)
            if abs(best_in_t_multiple_of_30 - minimum2) < function_domain/10.24:
                #print("Crestem probabilitatea de mutatie", best_in_t_multiple_of_30, minimum2 )
                mutation_probability = 0.003
            else:
                mutation_probability = 0.001
            best_in_t_multiple_of_30 = minimum2

        if t % 100 == 0:
            mini2 = np.inf
            for e in E:
                if e < mini2:
                    mini2 = e
            #print("Probabilitatea de mutatie este de: ", mutation_probability)
            print(mini2, "la iteratia curenta: ", t)


        # Adaugarea valorilor la liste
        mini3 = np.inf
        for e in E:
            if e < mini3:
                mini3 = e

        generations_list.append(t)
        function_values_list.append(mini3)


        elit = elitism(pt, E, k)


        pt = selection(pt, E)
        mutate(pt, no_of_bits)

        pt.extend(elit)

        crossover(pt, l, function, lower, upper)

        # print(E[-k:])  # print the last k elements of E

        for curent_pop_size in range(pop_size):
            E[curent_pop_size] = evaluate(pt[curent_pop_size], l, function, lower, upper)

    for e in E:
        if e < minimum:
            minimum = e

    end = time.time()
    execution_time = (end - start)
    execution_times.append(round(execution_time, 3))
    solutions_values.append(round(minimum, 5))

    return minimum


def print_time_statistics():
    avg_time = sum(execution_times) / len(execution_times)
    print("Average execution time:", avg_time)


def print_function_values_statistics():
    avg_solution = sum(solutions_values) / len(solutions_values)
    best_solution = min(solutions_values)
    print("Average solution:", avg_solution)
    print("Best solution:", best_solution)


if __name__ == "__main__":

    dimensions = [30]
    mutation_probability = 0.001
    crossover_probability = 0.8
    epsilon = 0.000001
    Cons = 30
    pop_size = 200
    k = (int)(pop_size * 0.07)
    generations = 2000

    RUNS = 1

    # MICHALEWICZ
    for dimension in dimensions:
        execution_times = []
        solutions_values = []
        for i in range(RUNS):
            ag(dimension, michalewicz, MICHALEWICZ_DOMAIN, MICHALEWICZ_LOWER, MICHALEWICZ_UPPER)
        print(f"MICHALEWICZ: {dimension} dimensions, {pop_size} pop_size, {generations} generations, {crossover_probability} crossover_probability, {mutation_probability} mutation_probability, {epsilon} epsilon.")
        print_time_statistics();
        print_function_values_statistics();

    # DE_JONG_1
    # for dimension in dimensions:
    #     execution_times = []
    #     solutions_values = []
    #     for i in range(RUNS):
    #         ag(dimension, de_jong_1, DE_JONG_1_DOMAIN, DE_JONG_1_LOWER, DE_JONG_1_UPPER)
    #     print(f"DE_JONG_1: {dimension} dimensions, {pop_size} pop_size, {generations} generations, {crossover_probability} crossover_probability, {mutation_probability} mutation_probability, {epsilon} epsilon.")
    #     print_time_statistics();
    #     print_function_values_statistics();


    # # RASTRIGIN
    # for dimension in dimensions:
    #     execution_times = []
    #     solutions_values = []
    #     for i in range(RUNS):
    #         ag(dimension, rastrigin, RASTRIGIN_DOMAIN, RASTRIGIN_LOWER, RASTRIGIN_UPPER)
    #     print(f"RASTRIGIN: {dimension} dimensions, {pop_size} pop_size, {generations} generations, {crossover_probability} crossover_probability, {mutation_probability} mutation_probability, {epsilon} epsilon.")
    #     print_time_statistics();
    #     print_function_values_statistics();

    # SCHWEFEL
    # for dimension in dimensions:
    #     execution_times = []
    #     solutions_values = []
    #     for i in range(RUNS):
    #         ag(dimension, schwefel, SCHWEFEL_DOMAIN, SCHWEFEL_LOWER, SCHWEFEL_UPPER)
    #     print(f"SCHWEFEL: {dimension} dimensions, {pop_size} pop_size, {generations} generations, {crossover_probability} crossover_probability, {mutation_probability} mutation_probability, {epsilon} epsilon.")
    #     print_time_statistics();
    #     print_function_values_statistics();

    # SCHWEFEL GRAFIC
    # # Afișarea graficului la finalul execuției
    # plt.scatter(generations_list, function_values_list, label='Schwefel Function', s=1)
    #
    # # Generarea tabloului pentru generații cu pas 1 până la 100 și pas logaritmic de la 100 încolo
    # generations_ticks = np.concatenate([[10, 50], np.logspace(2, np.log10(max(generations_list)), 10)])
    # plt.xticks(generations_ticks, rotation=45)  # Rotim marcajele pentru a le face mai lizibile
    #
    # #plt.xscale('log')  # Setarea scalei logaritmice pentru axa OX
    #
    # #Generarea tabloului pentru valorile funcției Rastrigin cu pas 1 până la 100 și pas logaritmic de la 100 încolo
    # function_values_ticks = [-4000, -6000, -7000, -8000, -9000, -10000, -10500, -11000, -11500, -12000]
    # plt.yticks(function_values_ticks)
    #
    # #plt.yscale('log')  # Setarea scalei logaritmice pentru axa OY
    #
    # plt.xlabel('Generations')
    # plt.ylabel('Minimum Values')
    # plt.legend()
    # plt.show()

    # RASTRIGIN GRAFIC
    # # Afișarea graficului la finalul execuției
    # plt.scatter(generations_list, function_values_list, label='Rastrigin Function', s=1)
    #
    # # Generarea tabloului pentru generații cu pas logaritmic
    # generations_ticks = np.concatenate([[10, 50], np.logspace(2, np.log10(max(generations_list)), 10)])
    # plt.xticks(generations_ticks, rotation=45)
    #
    # # Generarea tabloului pentru valorile funcției Rastrigin
    # function_values_ticks = [0.0, 20.0, 50.0, 75.0, 100.0, 116.4, 135.5, 157.7, 183.5, 213.6, 248.6, 289.3, 336.8, 392.0]
    # plt.yticks(function_values_ticks)
    #
    # plt.xlabel('Generations')
    # plt.ylabel('Minimum Values')
    # plt.legend()
    # plt.show()

    # DE JONG GRAFIC
    # # Afișarea graficului la finalul execuției
    # plt.scatter(generations_list, function_values_list, label='De Jong 1 Function', s=1)
    #
    # # Generarea tabloului pentru generații cu pas logaritmic
    # generations_ticks = np.concatenate([[10, 50], np.logspace(2, np.log10(max(generations_list)), 10)])
    # plt.xticks(generations_ticks, rotation=45)
    #
    # # Generarea tabloului pentru valorile funcției Rastrigin
    # function_values_ticks = [0.0, 20.0, 50.0, 75.0, 100.0, 104.8, 109.8, 115.1, 120.6, 126.4, 132.4, 138.8, 145.4,
    #                          152.4]
    # plt.yticks(function_values_ticks)
    #
    #
    # plt.xlabel('Generations')
    # plt.ylabel('Minimum Values')
    # plt.legend()
    # plt.show()

    # MICHALEWICZ GRAFIC
    # Afișarea graficului la finalul execuției
    plt.scatter(generations_list, function_values_list, label='Michalewicz Function', s=1)

    # Generarea tabloului pentru generații cu pas logaritmic
    generations_ticks = np.concatenate([[10, 50], np.logspace(2, np.log10(max(generations_list)), 10)])
    plt.xticks(generations_ticks, rotation=45)

    # Generarea tabloului pentru valorile funcției Rastrigin
    function_values_ticks = [-28, -10, -5, 0]
    plt.yticks(function_values_ticks)

    plt.xlabel('Generations')
    plt.ylabel('Minimum Values')
    plt.legend()
    plt.show()