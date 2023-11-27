import numpy as np

class MassSpringDamper:
    def __init__(self, mass, spring, damper, sim_time=30, sim_step=0.05, x0=1, v0=0):
        self.mass = mass
        self.spring = spring
        self.damper = damper
        self.sim_time = sim_time
        self.sim_step = sim_step
        self.x0 = x0
        self.v0 = v0

    def forward_euler_method(self):
        num_steps = int(self.sim_time / self.sim_step)
        time = np.linspace(0, self.sim_time, num_steps + 1)
        positions = np.zeros(num_steps + 1)
        velocities = np.zeros(num_steps + 1)

        positions[0] = self.x0
        velocities[0] = self.v0

        for i in range(1, num_steps + 1):
            # Forward Euler integration
            acceleration = - (self.damper * velocities[i - 1] + self.spring * positions[i - 1]) / self.mass
            velocities[i] = velocities[i - 1] + self.sim_step * acceleration
            positions[i] = positions[i - 1] + self.sim_step * velocities[i]

        return time, positions

    def backward_euler_method(self):
        num_steps = int(self.sim_time / self.sim_step)
        time = np.linspace(0, self.sim_time, num_steps + 1)
        positions = np.zeros(num_steps + 1)
        velocities = np.zeros(num_steps + 1)

        positions[0] = self.x0
        velocities[0] = self.v0

        for i in range(1, num_steps + 1):
            # Backward Euler integration using fsolve
            def implicit_eqn(x):
                return x - positions[i - 1] - self.sim_step * (self.damper * x + self.spring * x) / self.mass

            positions[i] = fsolve(implicit_eqn, positions[i - 1])[0]
            velocities[i] = (positions[i] - positions[i - 1]) / self.sim_step

        return time, positions
    
    def plot_simulation(self):
        return

    def animate_simulation(self):
        return

if __name__ == "__main__":
    system = MassSpringDamper(2,5,3)
