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
*p*(s,s’) = *q*(s,s’)*a*(s,s’)  
Where the *q*(s,s’) is proposal distribution and the *a*(s,s’) is acceptance distribution.  
The a(s,s’) can be expressed as  
*a*(s,s’) = min{1, p(s’)q(s’,s)/p(s)q(s,s’)}  
We can choose the proposal distribution p(s,s’) to be the Gaussian distribution form.  
q(s,s’)=N(s,sigma)  
q(s',s)=N(s',sigma')  
the s and s' is the spin direction at t and t+1 respectively, if we control the sigma and sigma' to be the same, the a(s,s') can be simplified as  
*a*(s,s’) = min{1, p(s’)/p(s)}  
because the proposal distribution q(s’,s) and q(s',s) are symmetric.  
Although the p(s) can not be calculted directly, the ratio of the state  p(s’)/p(s) can be expressed as  
p(s’)/p(s) = exp(-H(s)/kT)/exp(-H(s')/kT)  
which is easy to calculate

