import random
from rdkit import Chem # Ensure Chem is imported
from rdkit.Chem import Crippen, Descriptors, Draw, QED
try:
    from rdkit.Chem.SA_Score import sascorer
except ImportError:
    from rdkit.Contrib.SA_Score import sascorer
import matplotlib.pyplot as plt

# Define a basic set of characters that can appear in SMILES
# This is a simplification; real SMILES alphabets are more complex
# and depend on the chemical space being explored.
# User can extend this list based on their needs.
SMILES_CHARS = list("CNOSPClBrIF#()[]=123456789+-.cno")
def generate_random_smiles(length=20):
    """
    Generates a random (likely invalid) SMILES string of a given length.
    """
    if not SMILES_CHARS:
        raise ValueError("SMILES_CHARS must be defined to generate random SMILES.")
    return "".join(random.choice(SMILES_CHARS) for _ in range(length))

def mutate(smiles_string, mutation_rate):
    """
    Mutates a SMILES string.
    With probability `mutation_rate`, performs one mutation operation:
    1. Character replacement
    2. Character insertion
    3. Character deletion
    """
    if random.random() < mutation_rate:
        mutated_list = list(smiles_string)
        
        # Handle empty string case for mutation
        if not mutated_list:
            if SMILES_CHARS:
                return random.choice(SMILES_CHARS) # Insert a random char if empty
            return "" # Cannot mutate if no chars and no SMILES_CHARS

        mutation_type = random.choice(["replace", "insert", "delete", "branch", "ring"])
        if mutation_type == "ring":
            # This mutation inserts a ring.
            # insert_pos_1 is the start of the segment to be enclosed by ring_char
            # insert_pos_2 is the end of the segment (exclusive for slice) and also the index of the char to append.
            # Requires len(mutated_list) >= 1.
            if len(mutated_list) >= 1:
                insert_pos_1 = random.randint(0, len(mutated_list) - 1)
                insert_pos_2 = random.randint(insert_pos_1, len(mutated_list) - 1)

                ring_char_to_use = '1' # Default
                for i in range(1, 10): # Digits 1-9 from SMILES_CHARS
                    if str(i) not in mutated_list: # Find a digit not currently in the SMILES string
                        ring_char_to_use = str(i)
                        break
                
                segment_to_enclose = mutated_list[insert_pos_1 : insert_pos_2]
                char_to_append = mutated_list[insert_pos_2]
                mutated_list = mutated_list[:insert_pos_1] + [ring_char_to_use] + segment_to_enclose + [ring_char_to_use] + [char_to_append] + mutated_list[insert_pos_2+1:]
            else:
                mutation_type = "insert"
        if mutation_type == "replace" and mutated_list:
            idx_to_replace = random.randrange(len(mutated_list))
            mutated_list[idx_to_replace] = random.choice(SMILES_CHARS)
        if mutation_type == "insert" and SMILES_CHARS:
            insert_pos = random.randint(0, len(mutated_list))
            insert_char = random.choice(SMILES_CHARS)
            if insert_char in ["[]()"]:
                mutation_type = "branch"
            elif insert_char in ["123456789"]:
                mutation_type = "ring"
            else:
                mutated_list.insert(insert_pos, insert_char)
        if mutation_type == "delete":
            if len(mutated_list) > 1:
                del_pos = random.randrange(len(mutated_list))
                del mutated_list[del_pos]
            elif len(mutated_list) == 1 and SMILES_CHARS: # If only one char, replace it
                 mutated_list[0] = random.choice(SMILES_CHARS)
        if mutation_type == "branch" and len(mutated_list) > 0:
            # This mutation wraps a segment and appends a character.
            # insert_pos_1 is the start of the segment to be wrapped by '('
            # insert_pos_2 is the end of the segment (exclusive for slice) and also the index of the char to append.
            # Requires len(mutated_list) >= 1.
            insert_pos_1 = random.randint(0, len(mutated_list) - 1)
            # insert_pos_2 must be a valid index for mutated_list[insert_pos_2]
            insert_pos_2 = random.randint(insert_pos_1, len(mutated_list) - 1)
            
            segment_to_wrap = mutated_list[insert_pos_1 : insert_pos_2] # Can be empty if insert_pos_1 == insert_pos_2
            char_to_append = mutated_list[insert_pos_2]
            
            mutated_list = mutated_list[:insert_pos_1] + ['('] + segment_to_wrap + [')'] + [char_to_append] + mutated_list[insert_pos_2+1:]
        return "".join(mutated_list)
    return smiles_string # No mutation occurred

def crossover(parent1, parent2, crossover_rate):
    """
    Performs crossover between two parent SMILES strings using a random cut point for each.
    """
    if random.random() < crossover_rate:
        if not parent1 or not parent2: # Cannot crossover with an empty string effectively
            return parent1, parent2

        # Determine cut points, allowing for full string if length is 1
        p1_cut = random.randint(0, len(parent1)) if len(parent1) > 0 else 0
        p2_cut = random.randint(0, len(parent2)) if len(parent2) > 0 else 0
            
        child1 = parent1[:p1_cut] + parent2[p2_cut:]
        child2 = parent2[:p2_cut] + parent1[p1_cut:]
        return child1, child2
    return parent1, parent2

def tournament_selection(population_with_fitness, tournament_size):
    """
    Selects an individual using tournament selection.
    Assumes fitness is the first element of tuples in population_with_fitness,
    and higher fitness is better.
    """
    if not population_with_fitness:
        return None 
    
    actual_tournament_size = min(tournament_size, len(population_with_fitness))
    if actual_tournament_size == 0: # Should not happen if population_with_fitness is not empty
        return None
        
    tournament = random.sample(population_with_fitness, actual_tournament_size)
    tournament.sort(key=lambda x: x[0], reverse=True) # Sort by fitness, descending
    return tournament[0][1] # Return the SMILES string of the winner

def simple_genetic_algorithm(
    initial_population,
    fitness_function,
    generations,
    population_size,
    crossover_rate,
    mutation_rate,
    tournament_size,
    elitism_count=1,
    random_seed=42,
):
    """
    A simple genetic algorithm for SMILES string optimization.

    Args:
        initial_population (list): A list of initial SMILES strings.
        fitness_function (callable): A function that takes a SMILES string
                                     and returns a fitness score (higher is better).
                                     It should handle invalid SMILES (e.g., return a very low score or None).
        generations (int): The number of generations to run.
        population_size (int): The target size of the population.
        crossover_rate (float): The probability of crossover for a pair of parents.
        mutation_rate (float): The probability an individual undergoes mutation.
        tournament_size (int): The number of individuals in a selection tournament.
        elitism_count (int): Number of best individuals to carry over to the next generation.
        random_seed (int): Seed for the random number generator for reproducibility.

    Returns:
        dict: A dictionary containing the results of the GA run:
              - "best_smiles" (str): The best SMILES string found.
              - "best_fitness" (float): The fitness score of the best SMILES string.
              - "all_fitness_scores" (list): A list of lists, where each inner list
                                             contains the fitness scores of all individuals
                                             in that generation.
              - "best_fitness_history" (list): A list of the best fitness score in each generation.
    """
    if random_seed is not None:
        random.seed(random_seed)

    # Population initialization
    population = list(initial_population) if initial_population else []

    # Ensure population is of population_size
    if len(population) < population_size:
        num_to_generate = population_size - len(population)
        if not SMILES_CHARS and num_to_generate > 0 and not population : # Cannot generate if no chars and no base pop
            raise ValueError("SMILES_CHARS must be defined to generate random initial population if initial_population is empty.")
            
        for _ in range(num_to_generate):
            if population and SMILES_CHARS: # Base new individuals on existing ones if possible
                # Mutate a random existing individual significantly
                population.append(mutate(random.choice(population), 1.0))
            elif SMILES_CHARS: # Generate purely random individuals
                population.append(generate_random_smiles(random.randint(10,30))) # Default length range
            else: 
                # This case implies population is not empty but SMILES_CHARS is.
                # We can only rely on mutating existing population members.
                population.append(mutate(random.choice(population), 1.0)) # Fallback to mutating existing
    elif len(population) > population_size:
        population = random.sample(population, population_size) # Sample if too large

    current_best_smiles = None
    current_best_fitness = -float('inf')
    
    all_generation_fitness_scores = []
    best_fitness_per_generation = []

    for gen in range(generations):
        pop_with_fitness = []
        for smiles in population:
            fitness = fitness_function(smiles)
            if fitness is not None: # Allow fitness function to reject individuals
                pop_with_fitness.append((fitness, smiles))

        if not pop_with_fitness:
            print(f"Generation {gen+1}: No individuals with valid fitness scores. Stopping.")
            break

        pop_with_fitness.sort(key=lambda x: x[0], reverse=True)
        
        # Store fitness data for the current generation
        all_generation_fitness_scores.append([item[0] for item in pop_with_fitness])
        best_fitness_per_generation.append(pop_with_fitness[0][0])

        if pop_with_fitness[0][0] > current_best_fitness:
            current_best_fitness = pop_with_fitness[0][0]
            current_best_smiles = pop_with_fitness[0][1]
            print(f"Generation {gen+1}: New best -> Fitness: {current_best_fitness:.4f}, SMILES: {current_best_smiles}")
        else:
            # Log current generation's best if no improvement overall
            print(f"Generation {gen+1}: Best this gen -> Fitness: {pop_with_fitness[0][0]:.4f}, SMILES: {pop_with_fitness[0][1]}")

        new_population = []

        # Elitism
        if elitism_count > 0:
            elites = [ind[1] for ind in pop_with_fitness[:min(elitism_count, len(pop_with_fitness))]]
            new_population.extend(elites)

        # Generate new individuals
        while len(new_population) < population_size:
            parent1_smiles = tournament_selection(pop_with_fitness, tournament_size)
            parent2_smiles = tournament_selection(pop_with_fitness, tournament_size)

            if parent1_smiles is None or parent2_smiles is None: # Should only happen if pop_with_fitness is empty
                print(f"Warning: Parent selection failed in generation {gen+1}. Populating with random individuals.")
                # Fallback: if selection fails (e.g. pop_with_fitness became unexpectedly small)
                # add random or mutated individuals. For simplicity, break if no parents.
                if not pop_with_fitness: break 
                parent1_smiles = random.choice(pop_with_fitness)[1]
                parent2_smiles = random.choice(pop_with_fitness)[1]

            child1_smiles, child2_smiles = crossover(parent1_smiles, parent2_smiles, crossover_rate)
            
            new_population.append(mutate(child1_smiles, mutation_rate))
            if len(new_population) < population_size:
                new_population.append(mutate(child2_smiles, mutation_rate))
        
        population = new_population[:population_size] # Ensure population size

        if not population:
            print(f"Generation {gen+1}: Population became empty. Stopping.")
            break

    return {
        "best_smiles": current_best_smiles,
        "best_fitness": current_best_fitness,
        "all_fitness_scores": all_generation_fitness_scores,
        "best_fitness_history": best_fitness_per_generation,
    }
     
# Plotting functions
def plot_ga_results(results):
    if not results:
        print("No results to plot.")
        return

    plt.style.use('seaborn-v0_8-whitegrid') # Apply a nicer style
    # Visualize the best molecule found
    print("\nVisualizing the best molecule...")
    if results['best_smiles']:
        mol = Chem.MolFromSmiles(results['best_smiles'])
        if mol:
            img = Draw.MolToImage(mol, size=(300, 300)) # Generate image
            plt.figure(figsize=(4, 4)) # Create a new figure for the molecule
            plt.imshow(img)
            plt.axis('off') 
            plt.title(f"Best Molecule Structure\nSMILES: {results['best_smiles']}\nFitness: {results['best_fitness']:.2f}")
            plt.show() # Display the molecule plot
        else:
            print(f"Could not generate RDKit molecule object for SMILES: {results['best_smiles']}")
    else:
        print("No best SMILES found to visualize.")
    if results['best_fitness_history']:
        plt.figure(figsize=(12, 6))
        generations_ran = len(results['best_fitness_history'])

        # Plot 1: Best and Average Fitness per Generation
        plt.subplot(1, 2, 1)
        plt.plot(results['best_fitness_history'], marker='o', linestyle='-', color='dodgerblue', label='Best Fitness')
        
        if results['all_fitness_scores']:
            avg_fitness_per_generation = [
                sum(scores) / len(scores) for scores in results['all_fitness_scores'] if scores
            ]
            if avg_fitness_per_generation:
                plt.plot(avg_fitness_per_generation, marker='x', linestyle='--', color='green', label='Average Fitness')
        
        plt.title('Fitness Progression Over Generations')
        plt.xlabel(f'Generation (0-{generations_ran-1})')
        plt.ylabel('Fitness Score')
        plt.legend()
        plt.grid(True)

        # Plot 2: Fitness Distribution per Generation (Box Plot)
        if results['all_fitness_scores']:
            plt.subplot(1, 2, 2)
            valid_generation_scores = [gen_scores for gen_scores in results['all_fitness_scores'] if gen_scores]
            if valid_generation_scores:
                bp = plt.boxplot(valid_generation_scores, patch_artist=True, medianprops={'color':'black'})
                for patch in bp['boxes']:
                    patch.set_facecolor('lightblue')
                plt.title('Fitness Distribution per Generation')
                plt.xlabel(f'Generation (0-{generations_ran-1})')
                plt.ylabel('Fitness Score')
                step = max(1, generations_ran // 10 if generations_ran > 10 else 1)
                plt.xticks(ticks=range(1, generations_ran + 1, step), labels=[str(i-1) for i in range(1, generations_ran + 1, step)])
            else:
                plt.text(0.5, 0.5, 'No fitness distributions to plot.', ha='center', va='center')
        
        plt.tight_layout()
        plt.show()
    else:
        print("No fitness history to plot.")
        
# Some example fitness functions
def rdkit_logp_fitness(smiles):
    """
    Calculates the LogP for a given SMILES string using RDKit.
    Returns a very low score for invalid SMILES.
    """
    if not smiles:
        return -1000.0 # Penalize empty SMILES

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return -1000.0  # Invalid SMILES, return a very low fitness
    try:
        logp = Crippen.MolLogP(mol)
        return logp
    except Exception: # Catch any other RDKit errors during calculation
        return -1000.0
    
    
    
