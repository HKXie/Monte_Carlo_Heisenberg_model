# Monte_Carlo_Heisenberg_model
Based on Heisenberg model, do Monte Carlo simulation calculation for MnBi2Te4 to produce magnetization and susceptibility.  
MnBi2Te4 is a topological material with A-type antiferromagnetic (AFM) structure, the interaction between layers are positive(antiferromagnetic), and it's negetive in the layer(ferromagnetic),as the figure below shows: 
![image](https://github.com/HKXie/Monte_Carlo_Heisenberg_model/blob/master/images/Magnetic%20structure_1.png)  
The J' is ferromagnetic interaction, and the J is antiferromagnetic. D is Anisotropic energy. H is external magnetic field along z axis.  
The Hamiltonian for the spin system in MnBi2Te4 can be written as:  
![image](https://github.com/HKXie/Monte_Carlo_Heisenberg_model/blob/master/images/Heisenberg_model.png)  
The spin system follows a Boltzmann distribution, expressed as:  
![image](https://github.com/HKXie/Monte_Carlo_Heisenberg_model/blob/master/images/Boltzmann%20distribution.png)  
The beta is -1/kT, where k is Boltzmann constant, T is the temperature.  
We use Markov chain Monte Carlo (MCMC) and Metropolis-Hastings algorithm to simulate the magnetization at different temperature or magnetic field.


