"""
Spectrum Python for THz-TDS (TeraHertz Time Domain Spectroscopy)  
"""
import numpy as np
import scipy.constants as constant
from numpy.linalg import inv

# Constants
light_speed = constant.speed_of_light * 1e-6  # µm/ps


class PropagationTerms:
    """
    # PropagationTerms
    This class is used to calculate any terms related to propagation through the media
    It consider an homogeneus media with for constant n(ω)
    ## Note:
    The Default units are lenght μm time ps frequency THz
    """

    def __init__(self, n1, n2, w: float, d: float):
        """
        The functions used for calculating the terms of the transmission trough the material
        ### Parameters:
        * `n1` = Material's complex refractive index
        * `n2` = External's medium complex refractive index
        * `w` = Angular Frequency measurement
        * `d` = Material thickness
        """
        self.n1 = n1
        self.n2 = n2
        self.w = w
        self.d = d

    def reflection_interface(self):
        """Calculates the reflection coefficients from medium 1 to 2 and from medium 2 to 1 assuming normal incidence
        ### Output:
        * `(R12,R21)`
        """
        x = self.n1
        y = self.n2
        R12 = (x - y) / (x + y)
        R21 = (y - x) / (x + y)

        return R12, R21

    def transmission_interface(self):
        """
        Calculate the transmission coefficient from medium 1 to 2 and from 2 to 1
        """

        T12 = 2 * self.n1 / (self.n1 + self.n2)
        T21 = 2 * self.n2 / (self.n1 + self.n2)

        return T12, T21

    def propagation(self):
        """
        Propagation term
        """
        exp_arg = self.n1 * self.w * self.d / light_speed
        P = np.exp(1j * exp_arg)
        return P

    def fabry_perot(self, num_of_terms="all"):
        """
        Calculate the Frabry-Perot term of the transmission function
        * `num_of_terms = Fabry-Perot number of terms`
        """
        R12, R21 = self.reflection_interface()
        P = self.propagation()

        if num_of_terms == "all":
            FP = 1 / (1 - R12 * R21 * (P**2))
        else:
            FP = 0
            for i in range(num_of_terms, step=-1):
                FP += (R12 * R21 * (P**2)) ** i

        return FP

    def total_propagation_term(self):
        """
        Calculate the total propagation coefficient
        """

        FP = self.fabry_perot()
        T12, T21 = self.transmission_interface()
        P = self.propagation()

        return FP * P * (T12 * T21)

    def argument(self):
        """
        The total argument of propagation trough the material
        """
        P_arg = np.real(self.n1 * self.w * self.d / light_speed)
        T12, T21 = self.transmission_interface()
        # T_arg = np.angle(T12 * T21)  # + np.angle(T21)
        # FP_arg = np.angle(self.fabry_perot())
        return P_arg  # + FP_arg  # + T_arg

    def __str__(self) -> str:
        return f"R12, R21 = {self.reflection_interface()}\nT12,T21 = {self.transmission_interface()}"


def transmission(n1, n2, w, d):
    """
    ### Input
    * `n1` = Sample complex refractive index
    * `n2` = External medium complex refractive index
    * `w` = Angular frequency
    * `d` = Thickness
    ### Output
    * `t`, `arg` = The complex tranmission and the transmission argument
    """

    # Sample data
    sample = PropagationTerms(n1, n2, w, d)
    sample_propagation = sample.total_propagation_term()
    sample_argument = sample.argument()

    # Reference data
    reference = PropagationTerms(n2, n2, w, d)
    reference_propagation = reference.total_propagation_term()
    reference_argument = reference.argument()

    # Transmission
    t = sample_propagation / reference_propagation
    arg = sample_argument - reference_argument

    return t, arg


def error(Tc, Tc_arg, Tm, Tm_arg):
    """
    Error function returns a scalar value
    ### Input
    * `Tc` = calculated transmission
    * `Tc_arg` = argument of calculated transmission
    * `Tm` = measured transmission
    * `Tm_arg` = argument od measured transmission

    ### Output
    * `error_value` = error value
    """
    dp = np.log(np.abs(Tc)) - np.log(np.abs(Tm))
    dq = Tc_arg - Tm_arg

    error_value = dp**2 + dq**2

    return error_value


def hessian_matrix_2d(func, step):
    """
    Calculate the numerical hessian matrix of a 2D scalar function whithin a dimension 3x3
    The matrix is has a regular step distance
    ### Input
    * `func` = 3x3 matrix (the center value is the one calculated)
    * `step` = step size
    ### Examples

    ```python
    >>> A = [[1,2,3],[1,2,3],[1,2,3]]
    >>> H = hessian_matrix_2d(A,0.1)
    >>> print(H)
    ```
    """
    # First order derivatives
    funcy, funcx = np.array(np.gradient(func)) / step

    # Second order derivatives
    funcxy, funcxx = 2 * np.array(np.gradient(funcx)) / step
    funcyy, funcyx = 2 * np.array(np.gradient(funcy)) / step

    # Hessian matrix calculated at the center
    hessian = [[funcxx[1][1], funcxy[1][1]], [funcyx[1][1], funcyy[1][1]]]

    return hessian


def calculate_transmission_in_region(n1, n2, w, d, step):
    """
    This function calculates the matrix of the transmission T within the region of the complex refracive index `n1` = `n + 1j*k` in (n-step, n+step),(k-step, k+step)
    ### Input
    * `n1` = sample complex refractive index
    * `n2` = external medium complex refractive index
    * `w` = angular frequency
    * `d` = thickness
    * `step` = step size for the matrix

    ### Output
    * `T_matrix` = transmission matrix
    ### Examples
    """
    # Makes a 3x3 matrix of tuples (n,k) for computation within region

    N1 = n1 * np.ones((3, 3))
    N2 = n2 * np.ones((3, 3))

    # Make a B "step-matrix" 3X3 of tuples (n,k) for computating the new N1
    B = np.meshgrid([-step, 0, step], [-step, 0, step])
    B = B[0] + 1j * B[1]

    # Define the new N1 matrix
    N1 = np.complex128(N1) + B

    # Use N1 (sample) and N2 (reference) to calculate the transmission matrix
    T_matrix = transmission(N1, N2, w, d)

    # Return transmission matrix
    return T_matrix


def find_root(
    r, Tm, Tm_arg, w, d, max_num_of_iterations=10, step=0.0001, error_stop_value=0.05
):
    """
    Find the value for the complex refractive index n_complex = [n,k] using the Newton-Rhapson method
    ### Input
    * `r` = initial guess refractive index
    * `Tm` = The measured transmission
    * `Tm_arg` = measured transmission argument/phase
    * `w` = angular frequency ω
    * `d` = material's thickness
    * `max_num_of_iterations` = number of iterations
    * `step` = step size for the calculated transmission matrix

    It stops if the number of iterations is superior to `max_num_of_iterations`
    """

    # add new criteria for stopping the iteration from what we should get
    i = 0
    initial_value = r
    if not (r):
        r = 2 + 3j
    while i < max_num_of_iterations:
        # Calculate the theoretical transmission within the small region with the phase cycle
        Tc, Tc_arg = calculate_transmission_in_region(r, 1, w, d, step)

        # Calculate the error from the calculated and measured transmission in the small region (3x3 matrix)
        error_val = error(Tc, Tc_arg, Tm, Tm_arg)

        # Calculate the error gradient from the `error_val` (3x3 matrix)
        error_grad = np.array(np.gradient(error_val)) / step
        # transform into a 2x1 matrix
        error_grad = [error_grad[1][1][1], error_grad[0][1][1]]

        # Calculate the Hessian matrix of the error function in r = [n,k] region
        H = hessian_matrix_2d(error_val, step)

        # Inversion of the Hessian matrix
        H_inv = inv(H)

        # Calculate the next approximation for [n,k]
        r = [np.real(r), np.imag(r)]
        r -= np.dot(H_inv, error_grad)

        r = r[0] + 1j * r[1]

        if error_val[1][1] <= error_stop_value:
            return r

        i += 1

    return None


class Medium:
    """
    Class that charactirezes the material
    """

    def __init__(self, material_name, thickness) -> None:
        self.material = material_name
        self.thickness = thickness

    def refractive_index(
        self,
        measured_transmission: list,
        measured_transmission_phase: list,
        angular_frequency: list,
        intial_refractive_index: complex,
        max_number_of_iterations: int,
    ):
        """
        This function will calculate the
        ### Input
        * `measured_transmission` = The list of the measured transmission as a function of the angular frequency
        * `angular_frequency` = The angular frequency w
        * `intial_refractive_index` = The guess initial refrective index to start the Newton-Raphson method
        * `max_number_of_iterations` = This gives how many iterations will be done with Newton-Raphson method

        ### Output
        * `n`: list of complex refractive index as funcion of the angular frequency ω
        """

        # The list of refractive index
        n = []

        # Set the first refractive index to use as starting value
        n_value = intial_refractive_index

        for w, T, phase in zip(
            angular_frequency, measured_transmission, measured_transmission_phase
        ):
            # Find the refractive index for the specific angular frequency w and for the measured transmission T(w)
            n_value = find_root(
                n_value,
                T,
                phase,
                w,
                self.thickness,
                max_number_of_iterations,
            )
            n.append(n_value)

        return n

    def optical_conductivity(self):
        pass

    def dielectric_constant(self):
        pass

    def __str__(self):
        name = f"Material: {self.material}"
        thickness = f"Thickness: {self.thickness}"
        return name + thickness
