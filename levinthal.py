import random
import string
import time

TARGET_STRUCTURE = "ACDEFGHIKLMNPQRSTVWY" 
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"      
MAX_ATTEMPTS = 10000           
TARGET_LENGTH = len(TARGET_STRUCTURE)

def generate_random_structure(length):
    return ''.join(random.choice(AMINO_ACIDS) for _ in range(length))

def levinthal_simulation():
    attempts = 0
    start_time = time.time()
    
    while attempts < MAX_ATTEMPTS:
        random_structure = generate_random_structure(TARGET_LENGTH)
        attempts += 1
        
        # Check if the random structure matches the target
        if random_structure == TARGET_STRUCTURE:
            elapsed_time = time.time() - start_time
            print(f"Attempts: {attempts} Time: {elapsed_time:.2f} seconds.")
            return
        
    print(f"NO {MAX_ATTEMPTS} attempts.")

# Run the Levinthal paradox simulation
levinthal_simulation()
