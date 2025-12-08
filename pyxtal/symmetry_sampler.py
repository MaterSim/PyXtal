from typing import Counter
from pyxtal.symmetry import Group
import random

class Sampler:
    """
    A class to sample space group and Wyckoff position combinations
    for a given chemical composition and maximum atom count.
    """

    def __init__(self, composition=[1, 2], max_atoms=100, max_wp=20, verbose=False):
        """
        Initializes the Sampler with a given composition and atom limit.

        Args:
            composition (list): A list of chemical element symbols, e.g., [1, 2].
            max_atoms (int): The maximum number of atoms allowed in the unit cell.
            max_wp (int): The maximum number of Wyckoff positions to consider.
        """
        self.composition = composition
        self.max_atoms = max_atoms
        self.max_wp = max_wp
        self.spg_range = range(1, 231)
        self.verbose = verbose
        self.max_comp_per_spg = 50
        if self.verbose:
            print(f"Sampler initialized for {self.composition} with <= {self.max_atoms} atoms")

    def generate_combinations(self):
        """
        Generates all valid (space group, Wyckoff positions) combinations
        for the given composition and atom count constraints.
        """
        self.base_symmetries = {}
        count = 0
        for spg in self.spg_range:
            # Get all valid Wyckoff combinations for the given constraints
            g = Group(spg)
            Z_max = self.max_atoms // sum(self.composition)
            for z in range(1, Z_max + 1):
                comp = [n * z for n in self.composition]
                wps, _, ids = g.list_wyckoff_combinations(comp, numWp=(len(comp), self.max_wp))
                # sort ids by the length of each sublist to get the sorted index
                indices = sorted(range(len(ids)), key=lambda i: sum(len(x) for x in ids[i]))
                wps = [wps[i] for i in indices]
                ids = [ids[i] for i in indices]
                if ids:
                    for i, id in enumerate(ids):
                        if self.verbose:
                            print(f"Space Group {spg}, {wps[i]}")
                        frozen_id = tuple(tuple(sorted(x)) for x in id)
                        key = (spg, frozen_id)
                        self.base_symmetries[key] = 1
                        count += 1
                        if count == self.max_comp_per_spg:
                            print("Too many combinations, moving to next space group.")
                            break
                if count == self.max_comp_per_spg:
                    count = 0
                    break

        if self.verbose:
            print(f"Generated {len(self.base_symmetries)} valid (spg, wps) combinations.")


    def augment_data(self, train_data, weight=10):
        """
        Augments the data by calling get_subgroup_composition

        Args:
            train_data (list): The DataFrame containing (spg, wps) combinations.
            weight (int): The weight to assign to the new combinations.

        Returns:
            list: A list of generated pyxtal structure objects.
        """
        for data in train_data:
            self.add_composition(data, multiplier=weight)
            spg, wps = data
            g = Group(spg)
            sub_comps = g.get_subgroup_composition(wps,
                                                   max_atoms=self.max_atoms,
                                                   max_wps=self.max_wp)

            for sub_comp in sub_comps:
                # sum_comp: (125, [[2], [0, 1]])
                self.add_composition(sub_comp, multiplier=weight//2)

    def add_composition(self, data, multiplier=1):
        """
        Adds a new (spg, wps) combination to the base_symmetries.

        Args:
            data (tuple): A tuple containing (spg, wps).
            multiplier (int): The weight to assign to the new combination.
        """
        spg, wps = data
        frozen_id = tuple(tuple(sorted(x)) for x in wps)
        key = (spg, frozen_id)
        if key in self.base_symmetries.keys():
            self.base_symmetries[key] += multiplier
            if self.verbose:
                print(f"Updated existing composition: {data} with multiplier {multiplier}.")
        else:
            self.base_symmetries[key] = multiplier
            if self.verbose:
                print(f"Added new composition: {data} with multiplier {multiplier}.")

    def get_top_combinations(self, top_n=10):
        """
        Retrieves the top N (spg, wps) combinations based on their weights.

        Args:
            top_n (int): The number of top combinations to retrieve.

        Returns:
            list: A list of tuples containing the top (spg, wps) combinations.
        """
        sorted_combinations = sorted(self.base_symmetries.items(), key=lambda item: item[1], reverse=True)
        return sorted_combinations[:top_n]

    def sample(self, N=1):
        """
        Samples a (spg, wps) combination based on their weights.

        Args:
            N (int): The number of samples to draw.

        Returns:
            list: A list of sampled (spg, wps) combinations.
        """
        # Use random.choices for efficient weighted sampling
        population = list(self.base_symmetries.keys())
        weights = list(self.base_symmetries.values())

        if not population:
            return []

        return random.choices(population, weights=weights, k=N)

if __name__ == '__main__':
    # --- Example Usage ---

    # 1. Initialize the sampler for a composition (e.g., SiO2) with a max atom count
    sio2_sampler = Sampler(composition=[1, 2], max_atoms=120, max_wp=15, verbose=True)

    # 2. Generate all possible (spg, wps) combinations
    sio2_sampler.generate_combinations()

    # 3. (Optional) Augment data using subgroup compositions
    data = [(152, [[1], [0]]), (181, [[8], [1]])]
    sio2_sampler.augment_data(data, weight=1000)

    print("Top combinations based on weights:")
    top_combinations = sio2_sampler.get_top_combinations(top_n=20)
    for combo, weight in top_combinations:
        print(f"Combination: {combo}, Weight: {weight}")


    samples = sio2_sampler.sample(5000)

    # Count the frequency of each sampled combination
    sample_counts = Counter(samples)

    print("\nTop 20 most frequent sampled combinations:")
    sorted_samples = sample_counts.most_common(20)
    for combo, count in sorted_samples:
        print(f"Combination: {combo}, Count: {count}")
