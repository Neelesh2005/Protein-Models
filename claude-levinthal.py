import random
import string
import time
from typing import Tuple, Optional

class LevinthalSimulation:
    # Standard amino acids in single letter code
    AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
    
    def __init__(self, target_structure: str, max_attempts: int = 10000):
        """
        Initialize the Levinthal paradox simulation.
        
        Args:
            target_structure (str): The target protein structure to find
            max_attempts (int): Maximum number of attempts before giving up
        """
        self.target_structure = target_structure.upper()
        self.target_length = len(target_structure)
        self.max_attempts = max_attempts
        self._validate_input()
    
    def _validate_input(self) -> None:
        """Validate the input parameters."""
        if not self.target_structure:
            raise ValueError("Target structure cannot be empty")
        
        if not all(aa in self.AMINO_ACIDS for aa in self.target_structure):
            invalid_aas = set(self.target_structure) - set(self.AMINO_ACIDS)
            raise ValueError(f"Invalid amino acids in target structure: {invalid_aas}")
    
    def generate_random_structure(self) -> str:
        """Generate a random protein structure of the target length."""
        return ''.join(random.choice(self.AMINO_ACIDS) 
                      for _ in range(self.target_length))
    
    def calculate_similarity(self, structure: str) -> float:
        """
        Calculate similarity between a structure and the target.
        Returns percentage of matching positions.
        """
        matches = sum(1 for a, b in zip(structure, self.target_structure) if a == b)
        return (matches / self.target_length) * 100
    
    def run(self) -> Tuple[bool, dict]:
        """
        Run the Levinthal paradox simulation.
        
        Returns:
            Tuple[bool, dict]: Success status and statistics dictionary
        """
        attempts = 0
        start_time = time.time()
        best_similarity = 0
        best_structure = ""
        
        while attempts < self.max_attempts:
            random_structure = self.generate_random_structure()
            attempts += 1
            
            current_similarity = self.calculate_similarity(random_structure)
            
            # Update best found structure
            if current_similarity > best_similarity:
                best_similarity = current_similarity
                best_structure = random_structure
            
            # Check if we found the target
            if random_structure == self.target_structure:
                elapsed_time = time.time() - start_time
                stats = {
                    "success": True,
                    "attempts": attempts,
                    "time": elapsed_time,
                    "best_similarity": 100.0,
                    "best_structure": random_structure
                }
                return True, stats
        
        elapsed_time = time.time() - start_time
        stats = {
            "success": False,
            "attempts": attempts,
            "time": elapsed_time,
            "best_similarity": best_similarity,
            "best_structure": best_structure
        }
        return False, stats

def run_levinthal_simulation(target_structure: str, max_attempts: int = 10000) -> None:
    """
    Run the Levinthal simulation and print results.
    
    Args:
        target_structure (str): Target protein structure to find
        max_attempts (int): Maximum number of attempts
    """
    try:
        simulation = LevinthalSimulation(target_structure, max_attempts)
        success, stats = simulation.run()
        
        print("\nLevinthal Paradox Simulation Results:")
        print("-" * 40)
        print(f"Target Structure: {target_structure}")
        print(f"Attempts Made: {stats['attempts']:,}")
        print(f"Time Elapsed: {stats['time']:.2f} seconds")
        print(f"Best Similarity: {stats['best_similarity']:.2f}%")
        print(f"Best Structure Found: {stats['best_structure']}")
        
        if success:
            print("\nSuccess! Target structure found!")
        else:
            print(f"\nFailed to find target structure in {max_attempts:,} attempts.")
            
    except ValueError as e:
        print(f"Error: {e}")

# Example usage
if __name__ == "__main__":
    TARGET = "ACDEFGHIKLMNPQRSTVWY"
    run_levinthal_simulation(TARGET)