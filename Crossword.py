import random
import sys
import os

# ---------- Parameters ---------- #
population_size = 2000
number_of_generations = 200
mutation_rate = 0.15
grid_size = 20

# ---------- Chromosome grid representation ---------- #

# Build a grid representation from a chromosome
def build_grid(chromosome, grid_size):
    grid = [['-' for _ in range(grid_size)] for _ in range(grid_size)]
    
    for x, y, direction, word, in chromosome:
        if direction == 'horizontal':
            # Check if the word fits horizontally
            if y + len(word) <= grid_size:
                for i in range(len(word)):
                    grid[x][y + i] = word[i]
        else:
            if x + len(word) <= grid_size:
                for j in range(len(word)):
                    grid[x + j][y] = word[j]
                    
    return grid

# Print the crossword grid
def print_grid(grid):
    for row in grid:
        print(''.join(row))
    print('\n')
    
# Print the grids of the population
def print_population(population):
    for i, chromosome in enumerate(population):
        print(f'Grid {i+1}:')
        print_grid(build_grid(chromosome, grid_size))
    
# ---------- Initialize the population ---------- #
# Initialize the population of chromosomes
def initialize_population(population_size, word_list, grid_size):
    population = []
    
    while len(population) < population_size:
        while True:
            chromosome = []
            # Place all words must be place in the chromosome
            all_words_placed = True
            # Iterate over all words
            for word in word_list:
                placed = place_word_chromosome(chromosome, word, grid_size)
                if not placed:
                    all_words_placed = False
                    break
            
            if all_words_placed:
                population.append(chromosome)
                break
            
    return population

# Placing a word in a chromosme
def place_word_chromosome(chromosome, word, grid_size, max_attempts = 100):
    attemps = 0
    while attemps < max_attempts:
        # Generate random coordinates and direction for the word
        x = random.randint(0, grid_size - 1)
        y = random.randint(0, grid_size - 1)
        direction = random.choice(['horizontal', 'vertical'])
        
        # Check if the word fits in the grid 
        new_x, new_y = (x, y + len(word) - 1) if direction == 'horizontal' else (x + len(word) - 1, y)
        if in_boundaries(new_x, new_y, grid_size):
            chromosome.append((x, y, direction, word))
            return True
        
        attemps += 1
    
    return False

# Check if the coordinates are within the boundary
def in_boundaries(x, y, grid_size):
    return 0 <= x and x < grid_size and 0 <= y and y < grid_size

# ---------- Crossover ---------- #

# The crossover mechanism
def crossover(parent_1, parent_2):
    # Pick two random indices within the chromosome
    crossover_pivots = sorted(random.sample(range(1, len(parent_1) - 1), 2))

    # Create new child
    child_1 = parent_1[:crossover_pivots[0]] + parent_2[crossover_pivots[0]:crossover_pivots[1]] + parent_1[crossover_pivots[1]:]

    return child_1

# Select two parents randomly from a population
def select_parents(population):
    parent_1 = random.choice(population)
    parent_2 = random.choice(population)
    
    return parent_1, parent_2

# ---------- Mutation ---------- #

# The mutation mechanism
def mutation(chromosome, mutation_rate, grid_size, max_attempts=100):
    # If the mutation chance is less than the mutation_rate
    if random.random() < mutation_rate:
        # Pick a random word to mutate
        word_index = random.randint(0, len(chromosome) - 1)
        word = chromosome[word_index]
        
        # Find a new valid replacement
        attempts = 0
        
        while attempts < max_attempts:
            # Generate random coordinates and direction
            x = random.randint(0, grid_size - 1)
            y = random.randint(0, grid_size - 1)
            direction = 'horizontal' if random.random() < 0.5 else 'vertical'
            
            # Check if the word fits in the grid
            new_x, new_y = (x, y + len(word[3]) - 1) if direction == 'horizontal' else (x + len(word[3]) - 1, y)
            if in_boundaries(new_x, new_y, grid_size):
                chromosome[word_index] = (x, y, direction, word[3])
                break
            
            attempts += 1
        
    return chromosome

# ---------- Fitness test ---------- #

# Calculate the fitness of a chromosome
def calculate_fitness(chromosome, grid_size, word_list):
    score = 0.0

    # Check adjacent words and intersections (valid and invalid)
    invalid_intersections = 0
    neighbours = 0
    # Iterate over the chromosome
    for idx, current_word in enumerate(chromosome):
        # Crate new chromosme that has the rest of the words
        new_chromosome = chromosome[:idx] + chromosome[idx+1:]
        # Create the grid of the new chromosome
        grid = build_grid(new_chromosome, grid_size)
        
        x, y, direction, word = current_word
        
        # Iterate over the letter of the chosen word
        for i in range(len(word)):
            xi, yi = (x, y + i) if direction == 'horizontal' else (x + i, y)
            # Check intersections
            if grid[xi][yi] != '-':
                if grid[xi][yi] != word[i]:
                    invalid_intersections += 1
                else: # Valid intersection
                    score += 1
            else:
                # Check adjacent words
                x_up, y_up = (xi - 1, yi) if direction == 'horizontal' else (xi, yi + 1)
                x_down, y_down = (xi + 1, yi) if direction == 'horizontal' else (xi, yi - 1)
                if in_boundaries(x_up, y_up, grid_size):
                    if grid[x_up][y_up] != '-':
                        neighbours += 1
                if in_boundaries(x_down, y_down, grid_size):
                    if grid[x_down][y_down] != '-':
                        neighbours += 1
    
    # Check neighbours at the begining and at the end of each word
    grid = build_grid(chromosome, grid_size) # Create the grid
    # Iterate over each word
    for idx, current_word in enumerate(chromosome):
        x, y, direction, text = current_word
        word_length = len(text)
        # Check the start and the end of the word
        if direction == 'horizontal':
            if y > 0 and grid[x][y - 1] != '-':
                score -= 8
            if y + word_length < grid_size and grid[x][y + word_length] != '-':
                score -= 8
        else:
            if x > 0 and grid[x - 1][y] != '-':
                score -= 8
            if x + word_length < grid_size and grid[x + word_length][y] != '-':
                score -= 8
    
    # Calculate the number of connected word sets
    visited = set()
    number_of_connected_components = 0
    
    for i in range(grid_size):
        for j in range(grid_size):
            if grid[i][j] != '-' and (i, j) not in visited:
                dfs(grid, i, j, visited)
                number_of_connected_components += 1
    
    # Award and penalize the chromosome
    score -= invalid_intersections * 3 + neighbours * 4
    if invalid_intersections == 0:
        score += 10
    if neighbours == 0:
        score += 10
    if number_of_connected_components == 1:
        score += 20
    else:
        score -= pow(2, number_of_connected_components)
    return score

# DFS function
def dfs(grid, x, y, visited):
    # Base condition
    # If the cell is already visited or out of boundaries
    if (x, y) in visited or not(0 <= x < len(grid)) or not(0 <= y < len(grid)) or grid[x][y] == '-':
        return

    # Mark the current cell as visited
    visited.add((x, y))
    
    # Visited the horizontal and vertical cells
    dfs(grid, x + 1, y, visited) # Upper cell
    dfs(grid, x - 1, y, visited) # Bottom cell
    dfs(grid, x, y + 1, visited) # Right cell
    dfs(grid, x, y - 1, visited) # Left cell
    

# ---------- Genetic algorithm ---------- #
def genetic_algorithm(word_list, population_size, number_of_generations, mutation_rate, grid_size, elitism=population_size//2):
    # Create the initial population
    population = initialize_population(population_size, word_list, grid_size)
    
    for generation in range(number_of_generations):
        print_progress(generation+1, number_of_generations)
        # Pick the best 50% for the current population
        population = sorted(population, key=lambda x: calculate_fitness(x, grid_size, word_list), reverse=True)[:elitism]
        # Create the new offsprings
        offspring = []
        for _ in range(elitism):
            parent_1, parent_2 = select_parents(population)
            child = crossover(parent_1, parent_2)
            child = mutation(child, mutation_rate, grid_size)
            offspring.append(child)
        
        # Add the new offspring to the elite of the current population
        population.extend(offspring)

    print()
    # Pick the best chromosome
    population = sorted(population, key=lambda x: calculate_fitness(x, grid_size, word_list), reverse=True)
    best_chromosome = population[0]
    
    return best_chromosome

# Function to print the progress of each test case
def print_progress(generation, total_generations):
    sys.stdout.write(f'\rGeneration {generation}/{total_generations}')
    sys.stdout.flush()
    
# ---------- Input and Output ---------- #
def read_files(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]

# Main function
def main(input_dir, output_dir, num_files):
    # Create the output directory
    os.makedirs(output_dir)
    
    # Iterate over all input files and perform the algorithm
    for i in range(1, number_of_files + 1):
        input_file = os.path.join(input_dir, f'input{i}.txt')
        output_file = os.path.join(output_dir, f'output{i}.txt')

        # Initialize the list of words
        word_list = read_files(input_file)
        # Perform the genetic algorithm
        best_chromosome = genetic_algorithm(word_list, population_size, number_of_generations, mutation_rate, grid_size)
        # Write the output to the output file
        output = ""
        for element in best_chromosome:
            x, y, direction, word = element
            direc = 0 if direction == 'horizontal' else 1
            output += str(x) + " " + str(y) + " " + str(direc) + "\n"
        with open(output_file, 'w') as file:
            file.write(output)
        print(f'Output for file {i} generated.')

input_dir = 'inputs'  # Input directory
output_dir = 'outputs'  # Output directory
number_of_files = 100

main(input_dir, output_dir, number_of_files)
