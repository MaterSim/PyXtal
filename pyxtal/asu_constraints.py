"""
Module for handling Asymmetric Unit (ASU) constraints with complex conditions.

This module provides tools to:
1. Define ASU box constraints (x_min ≤ x ≤ x_max, etc.)
2. Define additional inequality constraints (e.g., x ≤ y, y ≤ 1/2 - x)
3. Check if coordinates satisfy ASU conditions
4. Generate random points within the ASU
5. Project coordinates into the ASU

Examples of ASU conditions:
- Space Group 1 (P1): 0 ≤ x ≤ 1, 0 ≤ y ≤ 1, 0 ≤ z ≤ 1
- Space Group 75 (P4): 0 ≤ x ≤ 1/2, 0 ≤ y ≤ 1/2, 0 ≤ z ≤ 1, with x ≤ y
- Space Group 146 (R3): x ≤ (1+y)/2, y ≤ min(1-x, (1+x)/2), z ≤ 1
"""

import numpy as np
from fractions import Fraction
from typing import List, Callable, Optional


class ASUCondition:
    """
    Represents a single inequality condition for ASU.
    
    Examples:
        x ≤ y  -->  lambda coords: coords[0] <= coords[1]
        y ≤ 1/2 - x  -->  lambda coords: coords[1] <= 0.5 - coords[0]
        x ≤ (1+y)/2  -->  lambda coords: coords[0] <= (1 + coords[1])/2
    """
    
    def __init__(self, condition: Callable, description: str):
        """
        Initialize an ASU condition.
        
        Args:
            condition: A callable that takes coords [x, y, z] and returns True/False
            description: A human-readable description of the condition
        """
        self.condition = condition
        self.description = description
    
    def check(self, coords: np.ndarray) -> bool:
        """
        Check if coordinates satisfy this condition.
        
        Args:
            coords: Array of shape (..., 3) with [x, y, z] coordinates
            
        Returns:
            Boolean or array of booleans indicating if condition is satisfied
        """
        return self.condition(coords)
    
    def __repr__(self):
        return f"ASUCondition({self.description})"


class ASU:
    """
    Represents the Asymmetric Unit for a space group with box constraints
    and additional inequality conditions.
    """
    
    def __init__(self, 
                 x_min: float, x_max: float,
                 y_min: float, y_max: float,
                 z_min: float, z_max: float,
                 conditions: Optional[List[ASUCondition]] = None):
        """
        Initialize ASU with box bounds and optional additional conditions.
        
        Args:
            x_min, x_max: x coordinate bounds
            y_min, y_max: y coordinate bounds
            z_min, z_max: z coordinate bounds
            conditions: List of additional ASUCondition objects
        """
        self.bounds = np.array([x_min, x_max, y_min, y_max, z_min, z_max])
        self.conditions = conditions or []
    
    @property
    def x_min(self): return self.bounds[0]
    
    @property
    def x_max(self): return self.bounds[1]
    
    @property
    def y_min(self): return self.bounds[2]
    
    @property
    def y_max(self): return self.bounds[3]
    
    @property
    def z_min(self): return self.bounds[4]
    
    @property
    def z_max(self): return self.bounds[5]
    
    def check_box_bounds(self, coords: np.ndarray) -> np.ndarray:
        """
        Check if coordinates are within the box bounds.
        
        Args:
            coords: Array of shape (..., 3) with [x, y, z] coordinates
            
        Returns:
            Boolean array indicating if each point is within bounds
        """
        coords = np.asarray(coords)
        if coords.ndim == 1:
            coords = coords.reshape(1, -1)
        
        in_bounds = (
            (coords[..., 0] >= self.x_min) & (coords[..., 0] <= self.x_max) &
            (coords[..., 1] >= self.y_min) & (coords[..., 1] <= self.y_max) &
            (coords[..., 2] >= self.z_min) & (coords[..., 2] <= self.z_max)
        )
        return in_bounds
    
    def check_conditions(self, coords: np.ndarray) -> np.ndarray:
        """
        Check if coordinates satisfy all additional conditions.
        
        Args:
            coords: Array of shape (..., 3) with [x, y, z] coordinates
            
        Returns:
            Boolean array indicating if each point satisfies all conditions
        """
        if not self.conditions:
            return np.ones(coords.shape[:-1], dtype=bool)
        
        coords = np.asarray(coords)
        result = np.ones(coords.shape[:-1], dtype=bool)
        
        for condition in self.conditions:
            result &= condition.check(coords)
        
        return result
    
    def is_valid(self, coords: np.ndarray) -> np.ndarray:
        """
        Check if coordinates are within the ASU (both box bounds and conditions).
        
        Args:
            coords: Array of shape (..., 3) with [x, y, z] coordinates
            
        Returns:
            Boolean array indicating if each point is in the ASU
        """
        return self.check_box_bounds(coords) & self.check_conditions(coords)
    
    def generate_random_points(self, n: int, max_attempts: int = 10000) -> np.ndarray:
        """
        Generate random points uniformly within the ASU.
        
        Args:
            n: Number of points to generate
            max_attempts: Maximum number of attempts per point
            
        Returns:
            Array of shape (n, 3) with valid ASU coordinates
        """
        points = []
        attempts = 0
        
        while len(points) < n and attempts < max_attempts:
            # Generate random points in the box
            candidates = np.random.rand(n - len(points), 3)
            candidates[:, 0] = candidates[:, 0] * (self.x_max - self.x_min) + self.x_min
            candidates[:, 1] = candidates[:, 1] * (self.y_max - self.y_min) + self.y_min
            candidates[:, 2] = candidates[:, 2] * (self.z_max - self.z_min) + self.z_min
            
            # Check which ones satisfy all conditions
            valid = self.check_conditions(candidates)
            points.extend(candidates[valid])
            
            attempts += 1
        
        if len(points) < n:
            raise ValueError(f"Could not generate {n} valid points after {max_attempts} attempts")
        
        return np.array(points[:n])
    
    def project_to_asu(self, coords: np.ndarray, method: str = 'clamp') -> np.ndarray:
        """
        Project coordinates to be within the ASU.
        
        Args:
            coords: Array of shape (..., 3) with [x, y, z] coordinates
            method: Projection method ('clamp' or 'modulo')
            
        Returns:
            Projected coordinates within the ASU
        """
        coords = np.asarray(coords, dtype=float).copy()
        original_shape = coords.shape
        coords = coords.reshape(-1, 3)
        
        if method == 'clamp':
            # Clamp to box bounds
            coords[:, 0] = np.clip(coords[:, 0], self.x_min, self.x_max)
            coords[:, 1] = np.clip(coords[:, 1], self.y_min, self.y_max)
            coords[:, 2] = np.clip(coords[:, 2], self.z_min, self.z_max)
            
            # Handle additional conditions (simple approach: iterative adjustment)
            for _ in range(10):  # Max iterations
                if np.all(self.check_conditions(coords)):
                    break
                
                for condition in self.conditions:
                    invalid = ~condition.check(coords)
                    if np.any(invalid):
                        # Simple heuristic: move slightly toward center of ASU
                        center = np.array([
                            (self.x_min + self.x_max) / 2,
                            (self.y_min + self.y_max) / 2,
                            (self.z_min + self.z_max) / 2
                        ])
                        coords[invalid] = 0.9 * coords[invalid] + 0.1 * center
        
        elif method == 'modulo':
            # Apply periodic boundary conditions first
            coords[:, 0] = coords[:, 0] % 1.0
            coords[:, 1] = coords[:, 1] % 1.0
            coords[:, 2] = coords[:, 2] % 1.0
            
            # Then clamp to ASU box
            coords[:, 0] = np.clip(coords[:, 0], self.x_min, self.x_max)
            coords[:, 1] = np.clip(coords[:, 1], self.y_min, self.y_max)
            coords[:, 2] = np.clip(coords[:, 2], self.z_min, self.z_max)
        
        return coords.reshape(original_shape)
    
    def __repr__(self):
        bounds_str = (f"x: [{self.x_min}, {self.x_max}], "
                      f"y: [{self.y_min}, {self.y_max}], "
                      f"z: [{self.z_min}, {self.z_max}]")
        if self.conditions:
            cond_str = ", conditions: " + ", ".join(str(c) for c in self.conditions)
        else:
            cond_str = ""
        return f"ASU({bounds_str}{cond_str})"


# ============================================================================
# Example ASU definitions with additional conditions
# ============================================================================

def create_asu_for_space_group(space_group: int) -> ASU:
    """
    Create ASU object for a given space group number.
    
    This is a demonstration with a few examples. For complete implementation,
    you would need to add all 230 space groups.
    
    Args:
        space_group: Space group number (1-230)
        
    Returns:
        ASU object with appropriate bounds and conditions
    """
    # Basic ASU bounds for common space groups (partial list for demonstration)
    # Full data can be loaded from pyxtal/database/asymmetric_unit.txt
    ASU_BOUNDS_DICT = {
        1: [0, 1, 0, 1, 0, 1],
        2: [0, 0.5, 0, 1, 0, 1],
        75: [0, 0.5, 0, 0.5, 0, 1],
        143: [0, 2/3, 0, 2/3, 0, 1],
        146: [0, 2/3, 0, 2/3, 0, 1/3],
        195: [0, 0.5, 0, 0.5, 0, 0.5],
    }
    
    bounds = ASU_BOUNDS_DICT.get(space_group, [0, 1, 0, 1, 0, 1])
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    
    # Define additional conditions for specific space groups
    conditions = []
    
    # Space Group 75 (P4): Additional condition x ≤ y
    if space_group == 75:
        conditions.append(ASUCondition(
            lambda c: c[..., 0] <= c[..., 1],
            "x ≤ y"
        ))
    
    # Space Group 143 (P3): Additional conditions
    elif space_group == 143:
        conditions.extend([
            ASUCondition(
                lambda c: c[..., 0] <= c[..., 1],
                "x ≤ y"
            ),
            ASUCondition(
                lambda c: c[..., 1] <= 0.5,
                "y ≤ 1/2"
            )
        ])
    
    # Space Group 146 (R3): More complex conditions
    elif space_group == 146:
        conditions.extend([
            ASUCondition(
                lambda c: c[..., 0] <= (1 + c[..., 1]) / 2,
                "x ≤ (1+y)/2"
            ),
            ASUCondition(
                lambda c: c[..., 1] <= np.minimum(1 - c[..., 0], (1 + c[..., 0]) / 2),
                "y ≤ min(1-x, (1+x)/2)"
            )
        ])
    
    # Space Group 195 (P23): Additional condition x ≤ y and y ≤ z
    elif space_group == 195:
        conditions.extend([
            ASUCondition(
                lambda c: c[..., 0] <= c[..., 1],
                "x ≤ y"
            ),
            ASUCondition(
                lambda c: c[..., 1] <= c[..., 2],
                "y ≤ z"
            )
        ])
    
    # Add more space groups as needed...
    
    return ASU(x_min, x_max, y_min, y_max, z_min, z_max, conditions)


# ============================================================================
# Utility functions
# ============================================================================

def parse_fraction(frac_str: str) -> float:
    """Convert fraction string like '1/2' to float."""
    frac_str = frac_str.strip()
    if '/' in frac_str:
        return float(Fraction(frac_str))
    return float(frac_str)


def test_asu_implementation():
    """Test the ASU implementation with examples."""
    print("Testing ASU implementation\n" + "=" * 60)
    
    # Test 1: Simple box constraint (Space Group 1)
    print("\n1. Space Group 1 (P1) - Simple box")
    asu1 = create_asu_for_space_group(1)
    print(asu1)
    
    test_points = np.array([
        [0.5, 0.5, 0.5],  # Valid
        [1.5, 0.5, 0.5],  # Invalid (x > 1)
        [0.0, 0.0, 0.0],  # Valid
    ])
    valid = asu1.is_valid(test_points)
    print(f"Test points validity: {valid}")
    
    # Test 2: Box with condition (Space Group 75)
    print("\n2. Space Group 75 (P4) - Box with x ≤ y")
    asu75 = create_asu_for_space_group(75)
    print(asu75)
    
    test_points = np.array([
        [0.3, 0.4, 0.5],  # Valid (x < y)
        [0.4, 0.3, 0.5],  # Invalid (x > y)
        [0.25, 0.25, 0.5],  # Valid (x = y)
    ])
    valid = asu75.is_valid(test_points)
    print(f"Test points validity: {valid}")
    
    # Test 3: Generate random points
    print("\n3. Generate random points in ASU 75")
    random_points = asu75.generate_random_points(5)
    print(f"Generated {len(random_points)} random points:")
    for i, pt in enumerate(random_points):
        print(f"  Point {i+1}: x={pt[0]:.4f}, y={pt[1]:.4f}, z={pt[2]:.4f}, valid={asu75.is_valid(pt)}")
    
    # Test 4: Project to ASU
    print("\n4. Project points to ASU")
    outside_points = np.array([
        [0.6, 0.3, 0.5],  # Outside condition (x > y)
        [1.5, 0.5, 0.5],  # Outside box
    ])
    projected = asu75.project_to_asu(outside_points)
    print(f"Original points:\n{outside_points}")
    print(f"Projected points:\n{projected}")
    print(f"Projected points valid: {asu75.is_valid(projected)}")
    

if __name__ == "__main__":
    test_asu_implementation()
