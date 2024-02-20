import numpy as np
import time

execution_times = []
solutions_values = []

# Constants for the lower and upper bounds of the functions
RASTRIGIN_LOWER = -5.12
RASTRIGIN_UPPER = 5.12

MICHALEWICZ_LOWER = 0
MICHALEWICZ_UPPER = np.pi

SPHERE_LOWER = -5.12
SPHERE_UPPER = 5.12

SUM_SQUARES_LOWER = -10
SUM_SQUARES_UPPER = 10


def rastrigin(x):
    # Rastrigin function
    sum_val = sum([(xi ** 2 - 10 * np.cos(2 * np.pi * xi)) for xi in x])
    return 10 * len(x) + sum_val


def michalewicz(x):
    d = len(x)
    _sum = 0

    for ii in range(1, d + 1):
        xi = x[ii - 1]
        new = np.sin(xi) * (np.sin(ii * xi**2 / np.pi))**(20)
        _sum += new

    y = -_sum
    return y


def sphere(x):
    return np.sum(np.square(x))


def sum_squares(x):
    return np.sum([(i + 1) * xi**2 for i, xi in enumerate(x)])


def function_minimum(function, lower, upper, deterministic, dimensions, iterations=100000):

    step = 0.01

    minim = float('inf')  # Initialize the searched minimum with infinite
    minim_points = []

    start = time.time()

    if deterministic:

        # Initialize the list of coordinates for the current point with values of 0
        current_point = [0] * dimensions

        # helper function
        def iterate_coordinates(index):

            # If we have finished generating coordinates in an iteration
            if index == dimensions:

                nonlocal minim, minim_points

                # Display the pair of values with at least 5 decimal places
                formatted_values = [f"{val:.5f}" for val in current_point]
                print(f"Pair of values: {formatted_values}")

                # calculating rastrigin function for current_point
                current_value = function(current_point)

                if current_value < minim:
                    minim = current_value
                    minim_points = [current_point.copy()]  # Reinitialize the list with the current point
                elif current_value == minim:
                    minim_points.append(current_point.copy())  # Adding the current point to the list

            else:
                # Iterating over the interval, incrementing the current point's coordinate
                for x in range(int((upper - lower) / step) + 1):
                    current_point[index] = lower + x * step
                    # Setting the next dimension
                    iterate_coordinates(index + 1)

        # generate all combinations of coordinates by calling the helper function
        iterate_coordinates(0)

    else:
        start = time.time()

        # run the code several times
        for i in range(iterations):
            # take a random point from the domain
            current_point = np.random.uniform(lower, upper, dimensions)
            # print(f"Pair of values: {current_point}")

            # calculating rastrigin function for current_point
            current_value = function(current_point)

            # finding the minimum value
            if current_value < minim:
                minim = current_value
                minim_points = [current_point.copy()]
            elif current_value == minim:
                minim_points.append(current_point.copy())

    # https://www.geeksforgeeks.org/how-to-check-the-execution-time-of-python-script/
    end = time.time()
    execution_time = (end - start)
    execution_times.append(round(execution_time, 3))
    return minim, minim_points


def print_time_statistics():
    # Calculate the minimum, maximum, and average values
    min_time = min(execution_times)
    max_time = max(execution_times)
    avg_time = sum(execution_times) / len(execution_times)

    print("Minimum execution time:", min_time)
    print("Maximum execution time:", max_time)
    print("Average execution time:", avg_time)


def print_function_values_statistics():
    # Calculate the minimum, maximum, and average values
    best_solution = min(solutions_values)
    worst_solution = max(solutions_values)
    avg_solution =  sum(solutions_values) / len(solutions_values)

    print("Worst solution:", worst_solution)
    print("Best solution:", best_solution)
    print("Average solution:", avg_solution)


def print_result(result):

    solutions_values.append(round(result[0], 5))

    print(f"Minimum value of the function: {result[0]}")
    print("Points that achieve the minimum:")

    # Iterate through the points that achieve the minimum
    for point in result[1]:
        formatted_values = [f"{val:.3f}" for val in point]
        print(f"Pair of values: {formatted_values}")
    print()


if __name__ == "__main__":

    # function to run each case to get the statistics about time and values
    for i in range(100):
        print_result(function_minimum(rastrigin, RASTRIGIN_LOWER, RASTRIGIN_UPPER, False, 2, 100000))  # nedeteministic, 2 dimensions, 100000 iterations
    print_time_statistics();
    print_function_values_statistics();

    # print_result(function_minimum(rastrigin, RASTRIGIN_LOWER, RASTRIGIN_UPPER, False, 2, 100000))   # nedeteministic, 2 dimensions, 100000 iterations
    # print_result(function_minimum(rastrigin, RASTRIGIN_LOWER, RASTRIGIN_UPPER, False, 10, 100000))  # nedeteministic, 10 dimensions, 100000 iterations
    # print_result(function_minimum(rastrigin, RASTRIGIN_LOWER, RASTRIGIN_UPPER, True, 2))    # deteministic, 2 dimensions
    # print_result(function_minimum(rastrigin, RASTRIGIN_LOWER, RASTRIGIN_UPPER, True, 10))   # deteministic, 10 dimensions
    #
    # print_result(function_minimum(michalewicz, MICHALEWICZ_LOWER, MICHALEWICZ_UPPER, False, 2, 100000))     # nedeteministic, 2 dimensions, 100000 iterations
    # print_result(function_minimum(michalewicz, MICHALEWICZ_LOWER, MICHALEWICZ_UPPER, False, 10, 100000))    # nedeteministic, 10 dimensions, 100000 iterations
    # print_result(function_minimum(michalewicz, MICHALEWICZ_LOWER, MICHALEWICZ_UPPER, True, 2))      # deteministic, 2 dimensions
    # print_result(function_minimum(michalewicz, MICHALEWICZ_LOWER, MICHALEWICZ_UPPER, True, 10))     # deteministic, 10 dimensions
    #
    # print_result(function_minimum(sphere, SPHERE_LOWER, SPHERE_UPPER, False, 2, 100000))     # nedeteministic, 2 dimensions, 100000 iterations
    # print_result(function_minimum(sphere, SPHERE_LOWER, SPHERE_UPPER, False, 10, 100000))    # nedeteministic, 10 dimensions, 100000 iterations
    # print_result(function_minimum(sphere, SPHERE_LOWER, SPHERE_UPPER, True, 2))      # deteministic, 2 dimensions
    # print_result(function_minimum(sphere, SPHERE_LOWER, SPHERE_UPPER, True, 10))     # deteministic, 10 dimensions
    #
    # print_result(function_minimum(sum_squares, SUM_SQUARES_LOWER, SUM_SQUARES_UPPER, False, 2, 100000))  # nedeteministic, 2 dimensions, 100000 iterations
    # print_result(function_minimum(sum_squares, SUM_SQUARES_LOWER, SUM_SQUARES_UPPER, False, 10, 100000)) # nedeteministic, 10 dimensions, 100000 iterations
    # print_result(function_minimum(sum_squares, SUM_SQUARES_LOWER, SUM_SQUARES_UPPER, True, 2))   # deteministic, 2 dimensions
    # print_result(function_minimum(sum_squares, SUM_SQUARES_LOWER, SUM_SQUARES_UPPER, True, 10))  # deteministic, 10 dimensions

