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
We use Markov chain Monte Carlo (MCMC) to simulate the magnetization at different temperature or magnetic field.  
We can't sample from the Boltzmann distribution directly for the whole spin system because the partial function it's difficult to calculate and the state space for system is infinite.  
We use the Metropolis-Hastings algorithm to generate the Markov chain for the spin system. The state transition function can be written as:  
*p*(x,x’) = *q*(x,x’)*a*(x,x’)  
Where the *q*(x,x’) is proposal distribution and the *a*(x,x’) is acceptance distribution.  
The a(x,x’) can be expressed as  
*a*(x,x’) = min{1, p(x’)q(x’,x)/p(x)q(x,x’)}  
We can choose the proposal distribution p(x,x’) to be the Gaussian distribution form.  
q(x,x’)=N(x,sigma)  
q(x',x)=N(x',sigma')  
the x and x' is the spin direction at t and t+1 respectively, if we control the sigma and sigma' to be the same, the a(x,x') can be simplified as  
*a*(x,x’) = min{1, p(x’)/p(x)}  
because the proposal distribution q(x’,x) and q(x',x) are symmetric

