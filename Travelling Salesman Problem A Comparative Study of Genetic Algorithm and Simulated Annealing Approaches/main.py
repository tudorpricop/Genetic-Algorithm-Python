import random
import math
import time
import csv

class Coords:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def get_distance(p1, p2):
    x_diff = p1.x - p2.x
    y_diff = p1.y - p2.y
    return int(math.sqrt(x_diff * x_diff + y_diff * y_diff) + 0.5)

def do_2opt(perm, left, right):
    while left < right:
        perm[left], perm[right] = perm[right], perm[left]
        left += 1
        right -= 1

class Chromosome:
    def __init__(self, p_size=None, nodes=None, p_perm=None):
        self.perm = []
        self.fitness = 0
        self.cost = 0
        self.elitist = False

        if p_perm is not None:
            self.perm = p_perm
        elif p_size is not None:
            self.fitness = self.cost = 0
            self.elitist = False

            self.perm = list(range(p_size))
            random.shuffle(self.perm)

            if nodes is not None:
                self.init_with_nodes(p_size, nodes)

    def init_with_nodes(self, p_size, nodes):
        self.fitness = self.cost = 0
        self.elitist = False

        improved = True
        while improved:
            improved = False

            for i in range(p_size - 1):
                for j in range(i + 1, p_size):
                    p_i1, p_i2 = self.perm[i], self.perm[(i + 1) % p_size]
                    p_j1, p_j2 = self.perm[j], self.perm[(j + 1) % p_size]
                    delta_len = get_distance(nodes[p_i1], nodes[p_i2]) + get_distance(nodes[p_j1], nodes[p_j2]) \
                                - get_distance(nodes[p_i1], nodes[p_j1]) - get_distance(nodes[p_i2], nodes[p_j2])

                    if delta_len > 0:
                        improved = True
                        do_2opt(self.perm, i + 1, j)

def init_random_generator():
    random.seed(int(time.time()))

def rand01():
    return random.random()

def read_data(file_name):
    nodes = []
    with open(file_name, 'r') as fin:
        for line in fin:
            if line.strip() == "NODE_COORD_SECTION":
                break
        for line in fin:
            if line.strip() == "EOF":
                break
            index, x, y = map(float, line.split())
            nodes.append(Coords(x, y))
    return nodes

def accept_worse_sol(delta, temperature):
    p = math.exp(-delta / temperature)
    return rand01() < p

def evaluate_dist(chr, nodes):
    chr.cost = 0
    for i in range(len(nodes) - 1):
        chr.cost += get_distance(nodes[chr.perm[i]], nodes[chr.perm[i + 1]])
    chr.cost += get_distance(nodes[chr.perm[-1]], nodes[chr.perm[0]])

def apply_displacement(p):
    mutated_perm = []
    left = random.randint(0, len(p.perm) - 1)
    right = random.randint(0, len(p.perm) - 1)
    if left > right:
        left, right = right, left
    subtour_size = right - left + 1
    insert_time = random.randint(0, len(p.perm) - subtour_size)
    inserted = False
    for i in range(len(p.perm)):
        if insert_time == 0 and not inserted:
            inserted = True
            mutated_perm.extend(reversed(p.perm[left:right + 1]))
        if i < left or i > right:
            insert_time -= 1
            mutated_perm.append(p.perm[i])
    if not inserted:
        mutated_perm.extend(p.perm[left:right + 1])
    p.perm = mutated_perm

def get_rand_neighbour(chr, nodes):
    ans = Chromosome(p_perm=chr.perm[:])
    apply_displacement(ans)
    evaluate_dist(ans, nodes)
    return ans

def simulated_annealing(chr, nodes, max_iterations=1000, cooling=0.99, epsilon=5):
    count = 0
    temperature = 100
    min_cost = chr.cost

    while temperature > 1e-epsilon:
        for _ in range(max_iterations):
            next_chr = get_rand_neighbour(chr, nodes)
            delta = chr.cost - next_chr.cost
            if delta < 0:
                delta *= -1
            if next_chr.cost < chr.cost:
                chr = next_chr
            elif accept_worse_sol(delta, temperature):
                chr = next_chr
            min_cost = min(min_cost, chr.cost)
        temperature *= cooling
        count += 1

    return min_cost

def evaluate_fitness(population, nodes, selection_pressure):
    max_cost = -1
    min_cost = -1
    for p in population:
        evaluate_dist(p, nodes)
        if max_cost == -1:
            max_cost = min_cost = p.cost
        else:
            max_cost = max(max_cost, p.cost)
            min_cost = min(min_cost, p.cost)
    if max_cost == min_cost:
        min_cost -= 1e-6
    for p in population:
        p.fitness = ((max_cost - p.cost) / (max_cost - min_cost) + 1) ** selection_pressure

def generate_offspring(offspring, parent, left, right):
    found = [False] * len(parent)
    for i in range(left, right + 1):
        found[offspring[i]] = True
    i = (right + 1) % len(parent)
    j = i
    while i != left:
        while found[parent[j]]:
            j = (j + 1) % len(parent)
        offspring[i] = parent[j]
        found[parent[j]] = True
        i = (i + 1) % len(parent)

def crossover(population, crossover_prob):
    random.shuffle(population)
    chr_size = len(population[0].perm)
    pop_size = len(population)
    for i in range(1, pop_size - 1, 2):
        if rand01() <= crossover_prob:
            left = random.randint(0, chr_size - 1)
            right = random.randint(0, chr_size - 1)
            if left > right:
                left, right = right, left
            parent1 = population[i].perm[:]
            parent2 = population[i + 1].perm[:]
            offspring1 = parent1[:]
            offspring2 = parent2[:]
            generate_offspring(offspring1, parent2, left, right)
            generate_offspring(offspring2, parent1, left, right)
            population[i].perm = offspring1
            population[i + 1].perm = offspring2

def select(population, cand_count, selected_count, pop_size):
    new_population = []
    if pop_size % selected_count != 0 or selected_count > cand_count:
        print("[ERROR]: not a valid selected_count")
        exit(-1)
    for _ in range(pop_size // selected_count):
        random.shuffle(population)
        population.sort(key=lambda x: x.fitness, reverse=True)
        for j in range(selected_count):
            new_population.append(population[j])
    population[:] = new_population

def mutate(population, mutation_prob):
    for p in population:
        if rand01() <= mutation_prob and not p.elitist:
            apply_displacement(p)

def genetic_search(nodes, pop_size=100, max_generations=2000, crossover_prob=0.3, c_tourney_pressure=50, selection_pressure=1):
    population = []
    size = len(nodes)
    c_mutation_prob = 0.1
    heuristic_count = 2 * pop_size // 10
    mutation_prob = c_mutation_prob
    not_found_cnt = 0

    for _ in range(heuristic_count):
        population.append(Chromosome(p_size=size, nodes=nodes))
    for _ in range(pop_size - heuristic_count):
        population.append(Chromosome(p_size=size))

    evaluate_fitness(population, nodes, selection_pressure)
    ans = population[0]
    ans.elitist = True

    generation = 0
    while generation < max_generations:
        select(population, c_tourney_pressure, c_tourney_pressure // 2, pop_size)
        mutate(population, c_mutation_prob)
        crossover(population, crossover_prob)
        evaluate_fitness(population, nodes, selection_pressure)
        generation += 1
        fittest = min(population, key=lambda x: x.cost)
        if fittest.cost < ans.cost:
            ans = fittest
            ans.elitist = True
        else:
            not_found_cnt += 1
        ans = fittest

    return simulated_annealing(ans, nodes, max_iterations=1000)

def get_optimal(instance):
    instances = ["berlin52", "st70", "kroA100", "d198", "kroB200", "pcb442", "pr439", "u2319", "pcb3038", "usa13509", "d18512"]
    optimal_val = [7542, 675, 21282, 15780, 29437, 50778, 107217, 234256, 137694, 19982889, 645488]
    for inst, opt in zip(instances, optimal_val):
        if instance == inst:
            return opt
    return 0

def create_csv(instance, sample_size):
    out_file_name = f"results/{instance}.csv"
    with open(out_file_name, 'a', newline='') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["functie", "dim", "minim", "timp"])

        best_res = 0
        results = []
        avg_time = 0
        avg_result = 0
        std_dev = 0

        for i in range(sample_size):
            print(f"Iteration {i + 1}/{sample_size}")
            t_start = time.time()

            nodes = read_data(instance + ".tsp")
            result = genetic_search(nodes)

            t_finish = time.time()
            elapsed = t_finish - t_start

            if i == 0:
                best_res = result
            else:
                best_res = min(best_res, result)

            avg_time += elapsed
            avg_result += result
            results.append(result)

            writer.writerow([instance, result, elapsed])

        avg_time /= sample_size
        avg_result /= sample_size

        for res in results:
            std_dev += (res - avg_result) ** 2
        std_dev = math.sqrt(std_dev / (sample_size - 1))

        opt_val = get_optimal(instance)
        precision = (1 - (avg_result - opt_val) / opt_val) * 100

        writer.writerow(["instance", "best result", "avg result", "optimum", "precision", "std dev", "avg time"])
        writer.writerow([instance, best_res, avg_result, opt_val, precision, std_dev, avg_time])

def main():
    import sys
    if len(sys.argv) != 3:
        print(f"{sys.argv[0]} <nume_fisier.tsp> <sample_size>")
        return -1

    sample_size = int(sys.argv[2])
    init_random_generator()
    print(sys.argv[1], sample_size)
    create_csv(sys.argv[1], sample_size)

if __name__ == "__main__":
    main()
