import numpy as np

class MC():
    def __init__(self, Layer_nums,L,J_intraSS, J_intraQQ ,J_interSS, J_interSQ, J_interQQ,K,B,T):
        self.Layer_nums=4#Layer_nums
        self.L=L
        # self.J_intra=J_intra
        # self.J_inter=J_inter

        self.J_intraSS=J_intraSS
        self.J_intraQQ=J_intraQQ

        self.J_interSS=J_interSS
        self.J_interQQ=J_interQQ
        self.J_interSQ=J_interSQ


        self.K=K
        self.B=B
        self.T=T
    
    # #system energy
    # def H(self,data_spin):
    #     for ee in range(self.Layer_nums):
    #         The_energy = 0
            
    #         for ii in range(self.L):
    #             for jj in range(self.L):

    #                 I=2*(ii%2)-1
    #                 E=2*(ee%2)-1
    #         #intralayer interaction
    #                 The_energy = The_energy + np.dot(data_spin[ee,ii,jj,:],    data_spin[ee,ii,jj-1,:]) * (-self.J_intra)
    #                 The_energy = The_energy + np.dot(data_spin[ee,ii,jj,:],    data_spin[ee,ii,(jj+1)%self.L,:]) * (-self.J_intra)
    #                 The_energy = The_energy + np.dot(data_spin[ee,ii,jj,:],    data_spin[ee,ii-1,jj,:]) * (-self.J_intra)
    #                 The_energy = The_energy + np.dot(data_spin[ee,ii,jj,:],    data_spin[ee,(ii+1)%self.L,jj,:]) * (-self.J_intra)
    #                 The_energy = The_energy + np.dot(data_spin[ee,ii,jj,:],    data_spin[ee,ii-1,(jj+I)%self.L,:]) * (-self.J_intra)
    #                 The_energy = The_energy + np.dot(data_spin[ee,ii,jj,:],    data_spin[ee,(ii+1)%self.L,(jj+I)%self.L,:]) * (-self.J_intra)
                    
    #         #interlayer interaction
    #                 The_energy = The_energy + np.dot(data_spin[ee,ii,jj,:],     data_spin[ee-1,ii,jj,:]) * (-self.J_inter)
    #                 The_energy = The_energy + np.dot(data_spin[ee,ii,jj,:],     data_spin[ee-1,(ii+E)%self.L,jj,:]) * (-self.J_inter)
    #                 The_energy = The_energy + np.dot(data_spin[ee,ii,jj,:],     data_spin[ee-1,(ii+E)%self.L,(jj+I)%self.L,:]) * (-self.J_inter)
    #                 The_energy = The_energy + np.dot(data_spin[ee,ii,jj,:],     data_spin[(ee+1)%self.Layer_nums,ii,jj,:]) * (-self.J_inter)
    #                 The_energy = The_energy + np.dot(data_spin[ee,ii,jj,:],     data_spin[(ee+1)%self.Layer_nums,(ii+E)%self.L,jj,:]) * (-self.J_inter)
    #                 The_energy = The_energy + np.dot(data_spin[ee,ii,jj,:],     data_spin[(ee+1)%self.Layer_nums,(ii+E)%self.L,(jj-I)%self.L,:]) * (-self.J_inter)

    #     # Anisotropy K and external magnetic field B
    #     The_energy = The_energy/2.0
    #     for ee in range(self.Layer_nums):
    #         for ii in range(self.L):
    #             for jj in range(self.L):
    #                 The_energy = The_energy - self.K*data_spin[ee,ii,jj,2]**2-self.B*data_spin[ee,ii,jj,2]

    #     return The_energy

    #The nearest energy for temporary_S at (ee,ii,jj) lattice site 
    def Nearest_energy(self,rnd_S,ee,ii,jj,data_spin): 
        # The_energy = 0
        
        # I=2*(ii%2)-1
        # E=2*(ee%2)-1

        # #intralayer
        # The_energy = The_energy + np.dot(temporary_S,     data_spin[ee,ii,jj-1,:])*(-self.J_intra)
        # The_energy = The_energy + np.dot(temporary_S,     data_spin[ee,ii,(jj+1)%self.L,:])*(-self.J_intra)
        # The_energy = The_energy + np.dot(temporary_S,     data_spin[ee,ii-1,jj,:])*(-self.J_intra)
        # The_energy = The_energy + np.dot(temporary_S,     data_spin[ee,(ii+1)%self.L,jj,:])*(-self.J_intra)
        # The_energy = The_energy + np.dot(temporary_S,     data_spin[ee,ii-1,(jj+I)%self.L,:])*(-self.J_intra)
        # The_energy = The_energy + np.dot(temporary_S,     data_spin[ee,(ii+1)%self.L,(jj+I)%self.L,:])*(-self.J_intra)

        # #interlayer
        # The_energy = The_energy + np.dot(temporary_S,     data_spin[ee-1,ii,jj,:])*(-self.J_inter)
        # The_energy = The_energy + np.dot(temporary_S,     data_spin[ee-1,(ii+E)%self.L,jj,:])*(-self.J_inter)
        # The_energy = The_energy + np.dot(temporary_S,     data_spin[ee-1,(ii+E)%self.L,(jj+I)%self.L,:])*(-self.J_inter)
        # The_energy = The_energy + np.dot(temporary_S,     data_spin[(ee+1)%self.Layer_nums,ii,jj,:])*(-self.J_inter)
        # The_energy = The_energy + np.dot(temporary_S,     data_spin[(ee+1)%self.Layer_nums,(ii+E)%self.L,jj,:])*(-self.J_inter)
        # The_energy = The_energy + np.dot(temporary_S,     data_spin[(ee+1)%self.Layer_nums,(ii+E)%self.L,(jj-I)%self.L,:])*(-self.J_inter)

        # The_energy = The_energy - self.K*temporary_S[2]**2-self.B*temporary_S[2]
	
        # return The_energy

        The_energy = 0
        
        I=2*(ii%2)-1
        E=2*(ee%2)-1

        if ee==2 or ee==3: #QL层

            #intralayer QQ
            The_energy = The_energy + rnd_S@data_spin[ee,ii,jj-1,:]*(-self.J_intraQQ)
            The_energy = The_energy + rnd_S@data_spin[ee,ii,(jj+1)%self.L,:]*(-self.J_intraQQ)
            The_energy = The_energy + rnd_S@data_spin[ee,ii-1,jj,:]*(-self.J_intraQQ)
            The_energy = The_energy + rnd_S@data_spin[ee,(ii+1)%self.L,jj,:]*(-self.J_intraQQ)
            The_energy = The_energy + rnd_S@data_spin[ee,ii-1,(jj+I)%self.L,:]*(-self.J_intraQQ)
            The_energy = The_energy + rnd_S@data_spin[ee,(ii+1)%self.L,(jj+I)%self.L,:]*(-self.J_intraQQ)

            #interlayer SQ or QQ
            The_energy = The_energy + rnd_S@data_spin[ee-1,ii,jj,:]*(-self.J_interSQ if ee==2 else -self.J_interQQ)
            The_energy = The_energy + rnd_S@data_spin[ee-1,(ii+E)%self.L,jj,:]*(-self.J_interSQ if ee==2 else -self.J_interQQ)
            The_energy = The_energy + rnd_S@data_spin[ee-1,(ii+E)%self.L,(jj+I)%self.L,:]*(-self.J_interSQ if ee==2 else -self.J_interQQ)
            The_energy = The_energy + rnd_S@data_spin[(ee+1),ii,jj,:]*(-self.J_interQQ if ee==2 else -self.J_interSQ)#%self.Layer_nums
            The_energy = The_energy + rnd_S@data_spin[(ee+1),(ii+E)%self.L,jj,:]*(-self.J_interQQ if ee==2 else -self.J_interSQ)#%self.Layer_nums
            The_energy = The_energy + rnd_S@data_spin[(ee+1),(ii+E)%self.L,(jj-I)%self.L,:]*(-self.J_interQQ if ee==2 else -self.J_interSQ)	#%self.Layer_nums

            The_energy = The_energy - self.K*rnd_S[2]**2-self.B*rnd_S[2]

        if ee==1 or ee==4: #SL层

            #intralayer SS
            The_energy = The_energy + np.dot(rnd_S,     data_spin[ee,ii,jj-1,:])*(-self.J_intraSS)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[ee,ii,(jj+1)%self.L,:])*(-self.J_intraSS)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[ee,ii-1,jj,:])*(-self.J_intraSS)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[ee,(ii+1)%self.L,jj,:])*(-self.J_intraSS)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[ee,ii-1,(jj+I)%self.L,:])*(-self.J_intraSS)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[ee,(ii+1)%self.L,(jj+I)%self.L,:])*(-self.J_intraSS)

            #interlayer  SQ
            The_energy = The_energy + np.dot(rnd_S,     data_spin[ee-1,ii,jj,:])*(-self.J_interSQ)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[ee-1,(ii+E)%self.L,jj,:])*(-self.J_interSQ)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[ee-1,(ii+E)%self.L,(jj+I)%self.L,:])*(-self.J_interSQ)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[(ee+1)%self.Layer_nums,ii,jj,:])*(-self.J_interSQ)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[(ee+1)%self.Layer_nums,(ii+E)%self.L,jj,:])*(-self.J_interSQ)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[(ee+1)%self.Layer_nums,(ii+E)%self.L,(jj-I)%self.L,:])*(-self.J_interSQ)


            #interlayer  SS
            The_energy = The_energy + np.dot(rnd_S,     data_spin[ee-3,ii,jj,:])*(-self.J_interSS)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[ee-3,(ii+E)%self.L,jj,:])*(-self.J_interSS)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[ee-3,(ii+E)%self.L,(jj+I)%self.L,:])*(-self.J_interSS)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[(ee+3)%self.Layer_nums,ii,jj,:])*(-self.J_interSS)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[(ee+3)%self.Layer_nums,(ii+E)%self.L,jj,:])*(-self.J_interSS)
            The_energy = The_energy + np.dot(rnd_S,     data_spin[(ee+3)%self.Layer_nums,(ii+E)%self.L,(jj-I)%self.L,:])*(-self.J_interSS)

            The_energy = The_energy - self.K*rnd_S[2]**2-self.B*rnd_S[2]
	
	
        return The_energy

    def Metropolis(self,data_spin,acccept_rate,sigma):

        if sigma <= 60 and acccept_rate < 1.0 :
            factor = 0.5/(1-acccept_rate)
            sigma = factor*sigma
        else:
            sigma = 60

        # spin_result = data_spin
        # create a counter to record acccept_rate
        counter = 0
        all = 0
        # rejected_counter=0
        for ee in range(1, self.Layer_nums+1):
            for ii in range(self.L):
                for jj in range(self.L):

                    if data_spin[ee,ii,jj,0]==0 and data_spin[ee,ii,jj,1]==0 and data_spin[ee,ii,jj,2]==0:
                        pass
                    
                    else:
                    
                        #generate the temporary spin with random direction
                        temporary_S = np.zeros(3)
                        gama = np.random.randn(3)

                        temporary_S[0] = gama[0]*sigma+data_spin[ee,ii,jj,0]
                        temporary_S[1] = gama[1]*sigma+data_spin[ee,ii,jj,1]
                        temporary_S[2] = gama[2]*sigma+data_spin[ee,ii,jj,2]
                        temporary_S = temporary_S/np.sqrt(np.dot(temporary_S,temporary_S))

                        data_S = data_spin[ee,ii,jj,:]
                        delta_E = self.Nearest_energy(temporary_S,ee,ii,jj,data_spin)-self.Nearest_energy(data_S,ee,ii,jj,data_spin)

                        if delta_E < 0:
                            data_spin[ee,ii,jj,:] = temporary_S
                            counter += 1
                            all+=1

                        elif np.random.rand()<np.exp(-delta_E/self.T):
                            data_spin[ee,ii,jj,:] = temporary_S
                            counter += 1
                            all+=1
                        else:
                            all += 1

        acccept_rate = counter/float(all) #1-counter/float(self.Layer_nums*self.L**2)#self.Layer_nums*self.L**2
        
        return data_spin, acccept_rate, sigma 


    def Summation(self,data_spin):
        M_x = 0
        M_y = 0
        M_z = 0
        for ee in range(1, self.Layer_nums+1):
            for ii in range(self.L):
                for jj in range(self.L):
                    
                    M_x += data_spin[ee,ii,jj,0]
                    M_y += data_spin[ee,ii,jj,1]
                    M_z += data_spin[ee,ii,jj,2]

        M_x = M_x/(self.L**2)#Layer_nums*
        M_y = M_y/(self.L**2)#Layer_nums*
        M_z = M_z/(self.L**2)#Layer_nums*
        M_total=[M_x,M_y,M_z]
        
        return np.sqrt(np.dot(M_total,M_total))*np.sign(M_z)#[M_x,M_y,M_z]

    def Summation_Z(self,data_spin):
        
        M_z = 0
        for ee in range(1, self.Layer_nums+1):
            for ii in range(self.L):
                for jj in range(self.L):
    
                    M_z += data_spin[ee,ii,jj,2]
        
        M_z = M_z/(self.L**2)#Layer_nums*

        return M_z