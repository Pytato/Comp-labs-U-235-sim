import numpy as np
import neutrons


class TopologyError(Exception):
    pass


class Uranium:
    """Generates Uranium Object, can be 1D or 3D based on 'cubey'.

    arg: l: float: length of 1D strip, side lengths of U-235 object in 'm',
        radius of sphere.

    kwarg: shape="line": string: defined which of the preset shapes this class
        represents for fissions.
    kwarg: a=1.7e-2: float: mean free path between neutron hitting nucleus
        and bouncing off it in 'm'.
    kwarg: b=21.0e-2: float: mean free path between neutron hitting nucleus
        and causing a fission event in 'm'.
    """

    def __enter__(self):
        return self

    def __init__(self, l, shape="line", a=1.7e-2, b=21.0e-2):
        # Declaring these parameters as class wide
        length = float(l)
        self.possible_shapes = ["line", "cube", "sphere"]
        self.x = length
        self.shape = shape.lower()

        if self.shape not in self.possible_shapes:
            raise TopologyError("kwarg shape as defined: {0} not in {1}".format(
                shape.lower(), self.possible_shapes))

        if self.shape == "cube":
            self.y = length
            self.z = length
        elif self.shape == "line":
            self.y = 0.0
            self.z = 0.0
        elif self.shape == "sphere":
            self.radius = length

        self.a = a
        self.b = b

    def this_neutron_is_important(self, decay_pos):
        """Generates the path a neutron will take within the U-235 object
        represented by this class and returns whether or not it exists within
        the U-235 after travelling it's average distance.

        #NeutronDecaysMatter
        """

        if self.shape == "line":
            # Figure out whether neutrons are moving left or right
            n_distance_travelled = 0
            while n_distance_travelled == 0:
                n_distance_travelled = np.random.randint(-1, 2)
            # Calculate the actual distance travelled for neutron
            mean_free_path = np.sqrt(2 * self.a * self.b)
            n_distance_travelled *= mean_free_path
            n_final_pos = decay_pos + n_distance_travelled

            return 0 <= n_final_pos[0] <= self.x

        # Looks like we have got something 3D on our hands
        # Generate a spherical coordinate vector for the direction
        # the product neutron we're considering will move.
        phi = 2.0 * np.pi * np.random.random()
        theta = np.arccos(2.0 * np.random.random() - 1.0)
        radius = neutrons.diffusion() * np.sqrt(2 * self.a * self.b)

        if self.shape == "cube":
            # Find the final position of the neutron in cartesian
            n_final_pos = np.add(
                np.array(decay_pos),
                np.array([radius * np.sin(theta) * np.cos(phi),
                          radius * np.sin(theta) * np.sin(theta),
                          radius * np.cos(theta)]))

            upper_bound = np.array([self.x, self.y, self.z])

            # Test whether it sits within all three bounds of the cube
            for i in range(3):
                if not 0 < n_final_pos[i] < upper_bound[i]:
                    return False

            return True

        if self.shape == "sphere":
            # Calculate the final position in cartesian coords
            n_final_pos = np.add(
                decay_pos,
                np.array([radius * np.sin(theta) * np.cos(phi),
                          radius * np.sin(theta) * np.sin(phi),
                          radius * np.cos(theta)]))

            # Return if the sqrt of the sum of squares of each
            # axis is less than the radius of the sphere.
            return np.linalg.norm(n_final_pos) < self.radius

    def fission(self, n_decays):
        """Complete a fission run on the U-235 object

        arg: int: n_decays: number of initial decays occurring within the
        U-235 length.

        returns:
            tuple(int, float):
                (secondary_fission_count, secondary_count/n_decays).
        """

        secondary_fission_count = 0

        # Now we loop for every single decay and figure out if its
        # neutrons will do anything interesting.

        if self.shape == "line" or self.shape == "cube":
            for i in range(n_decays):
                decay_pos = np.array([self.x * np.random.random(),
                                      self.y * np.random.random(),
                                      self.z * np.random.random()])

                # Generating a value for number of neutrons produced.
                n_neutrons_produced = neutrons.neutrons()

                # For each neutron produced, calculate travel distance
                for j in range(n_neutrons_produced):
                    if self.this_neutron_is_important(decay_pos):
                        secondary_fission_count += 1

        if self.shape == "sphere":
            for i in range(n_decays):
                # Find out how many neutrons are produced by this event.
                n_neutrons_produced = neutrons.neutrons()
                '''
                decay_pos = np.multiply(np.random.rand(3), 2 * self.radius)

                while np.linalg.norm(decay_pos) > self.radius:
                    decay_pos = np.multiply(np.random.rand(3), 2 * self.radius)

                '''
                # Generate a spherical coord position vector to represent the
                # initial neutron event, in phi, radius and theta
                phi = 2.0 * np.pi * np.random.random()
                theta = np.arccos(2 * np.random.random() - 1)
                radius = self.radius * np.cbrt(np.random.random())
                
                # Calculate the cartesian coordinates for the decay event
                decay_pos = np.array([radius * np.sin(theta) * np.cos(phi),
                                      radius * np.sin(theta) * np.sin(phi),
                                      radius * np.cos(theta)])

                # For each neutron produced, check if it remains within the sphere.
                for j in range(n_neutrons_produced):
                    if self.this_neutron_is_important(decay_pos):
                        secondary_fission_count += 1

        return secondary_fission_count, (secondary_fission_count/float(n_decays))

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


if __name__ == "__main__":
    with Uranium(10.0e-2) as uranium_c:
        print("Line test:", uranium_c.fission(10000))

    with Uranium(10.0e-2, shape="cube") as uranium_c:
        print("Cube test:", uranium_c.fission(10000))

    with Uranium(10.0e-2, shape="sphere") as uranium_c:
        print("Sphere test:", uranium_c.fission(10000))

'''
local_average_input = np.array([])
np.append(local_average_input, uranium_c.fission(n_decays=100))
'''
