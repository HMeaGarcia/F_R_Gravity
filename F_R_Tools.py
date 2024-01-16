import numpy as np






class ModelFR(object):

    """
    Base class for the RHS of the f(R) equations assuming 
    Spherically Simmetric Spacetimes

    RHS Equations:

    There are two equivalent formulations of the RHS of the SSS f(R)
    equations. 

    ----------
    
     first_test: Corresponds to the Eqs. 4.30 of the Thesis.
     second_test

     second_test: Uses an aditional equation for the n(r) variable.


    ----------
    alpha : float
        Dimensionless parameter .
    status : string
        Current status of the solver: 'running', 'finished' or 'failed'.
    t_bound : float
        Boundary time.

    """
    


    def __init__(self, F , alpha, beta, rho=None):
        
        if not callable(F):
            raise TypeError('F is %s, not a function' % type(F))
        
        self.F = F
        self.rho = rho 
        self.alpha = alpha
        self.beta = beta
    
    def first_test(self, r, z):
        alpha, beta, rho = self.alpha, self.beta, self.rho
        R,Rprime, m, n, p = z
        f, f_R, f_RR, f_RRR  = self.F(R)
        
        """Returns a vector with the righ-hand side of the system of equations
        for the model selected. 
        """
        
        
        Tbar = 3*p - rho
        T_rr = p
        T_tt = -rho
    
        A = (3 * f_RR )**(-1)
        B = m * (alpha * Tbar + beta*(2 * f - R * f_R) ) - 3 * f_RRR * Rprime**2
    
        mprime1 = m/(r*(2 * f_R + r * Rprime* f_RR))
        mprime2 = 2 * f_R * (1 - m) - 2* alpha*m * r**2*T_tt + m*r**2/3*( beta*(R* f_R + f) + 2*alpha*Tbar)
        mprime3 = r*Rprime*f_RR/f_R
        mprime4 = m*r**2/3*( beta*(2 * R * f_R - f) + alpha*Tbar) - alpha*m*r**2*(T_tt + T_rr) + 2*(1 - m )* f_R + 2*r* Rprime * f_RR
        MPRIME = mprime1*(mprime2 +  mprime3*mprime4)
    
        nprime1 = n/(r*(2* f_R + r* Rprime * f_RR))
        nprime2 = m*r**2*( beta*(f - R * f_R) + 2*alpha*T_rr) + 2*f_R*(m-1)- 4*r*Rprime * f_RR
        NPRIME = nprime1 * nprime2 
        
        C = (MPRIME/(2*m) - NPRIME/(2*n) - 2/r)*Rprime
    
        H1 = Rprime 
        H2 = A*B+ C
        H3 = MPRIME
        H4 = NPRIME
        H5 = -(p + rho)* NPRIME/(2*n)
        H = [H1,H2,H3,H4,H5]
    
        return H
    
    def second_test(self, r, z):

        alpha, beta, rho = self.alpha, self.beta, self.rho
        R, Rprime, m, n, NPRIME, p = z
        f, f_R, f_RR, f_RRR  = self.F(R)

        """Returns a vector with the righ-hand side of the system of equations
        for the second set of equations, including the extra one
        """
        Tbar = 3*p - rho
        T_rr = p
        T_tt = -rho
        T_tetateta= p 
    
        A = (3 * f_RR)**(-1)
        B = m*(alpha/2*Tbar + beta*(2* f - R * f_R ) - 3*f_RRR*Rprime**2)
    
        mprime1 = m/(r*(2* f_R + r*Rprime*f_RR ))
        mprime2 = 2*f_R*(1 - m) - alpha*m*r**2*T_tt + m*r**2/3*(beta*(R*f_R + f ) + alpha*Tbar)
        mprime3 = r*Rprime*f_RR/f_R
        mprime4 = m*r**2/3*(beta*(2*R*f_R - f) + alpha/2*Tbar) - alpha/2*m*r**2*(T_tt + T_rr ) + 2*(1 - m )*f_R + 2*r*Rprime*f_RR
        MPRIME = mprime1*(mprime2 +  mprime3*mprime4)
   
        n_pprime1 = 2*n*m/f_R*(T_tetateta*alpha/2 -1/6*(beta*(R*f_R + f ) + alpha*Tbar) + Rprime*f_RR/(r*m))
        n_pprime2 = n/(2*r)*(2*(MPRIME/m - NPRIME/n) + NPRIME*r/n*(MPRIME/m + NPRIME/n))
        N_PPRIME =  n_pprime1 +  n_pprime2
    
        C = (MPRIME/(2*m) -NPRIME/(2*n) - 2/r)*Rprime
    
        H1 = Rprime 
        H2 = A*B+ C
        H3 = MPRIME
        H4 = NPRIME
        H5 = N_PPRIME
        H6 = -(p + rho)* NPRIME/(2*n)
    
        f = [H1,H2,H3,H4,H5,H6]
        return  np.array(f)
    
    def beta_coef(self,R): 
        """
        Returns the coefficient equivalent to 
        the derivative of the potential
        """
        f, f_R, f_RR, f_RRR  = self.F(R)
        return (2 *f - R * f_R)
    
    def beta_coef_wtrace(self,R, p=None): 
        """
        Returns the coefficient plus the trace of the EMT
        """
        alpha, beta, rho = self.alpha, self.beta, self.rho
        f, f_R, f_RR, f_RRR  = self.F(R)
        Tbar = 3*p - rho
        return (alpha*Tbar + 2 *f - R * f_R)
    

    def effe_mass(self, R, p = None):
        """
        Returns the coefficient plus the trace of the EMT and 
        the extra term in the Eqs
        """
        alpha, beta, rho = self.alpha, self.beta, self.rho
        f, f_R, f_RR, f_RRR  = self.F(R)
        Tbar = 3*p - rho 
        return (f_R - R*f_RR)/(3*f_RR) - f_RRR/f_RR*((2 *f - R * f_R)/3*f_RR)



class STBS(ModelFR):
    """
    Class for the Starobinsky f(R) Model as defined in:

         Starobinsky, A. A. “Disappearing Cosmological Constant in f(R) Gravity”. 
         JETP Letters, 86, (2007), 157–163

    PARAMETERS
    ----------
    q : integer
        Dimensionless parameter.
    llambda : float
        Dimensionless parameter.
    R_tilde : float
        Dimensionless parameter  
    """
    
    def __init__(self, R_tilde,llambda, q):
        self.R_tilde = R_tilde
        self.llambda =  llambda
        self.q = q
    
    def __call__(self, R): 
        llambda =  self.llambda
        R_S = self.R_tilde
        q = self.q
        f = R + llambda*R_S* ((1  +R**2/R_S**2)**(-q) - 1)
        f_R = -2*R*llambda*q*(R**2/R_S**2 + 1)**(-q)/(R_S*(R**2/R_S**2 + 1)) + 1
        f_RR = 2*llambda*q*(R**2/R_S**2 + 1)**(-q)*(2*R**2*q/(R_S**2*(R**2/R_S**2 + 1)) + 2*R**2/(R_S**2*(R**2/R_S**2 + 1)) - 1)/(R_S*(R**2/R_S**2 + 1))
        f_RRR =  4*R*llambda*q*(R**2/R_S**2 + 1)**(-q)*(-2*R**2*q**2/(R_S**2*(R**2/R_S**2 + 1)) - 6*R**2*q/(R_S**2*(R**2/R_S**2 + 1)) - 4*R**2/(R_S**2*(R**2/R_S**2 + 1)) + 3*q + 3)/(R_S**3*(R**2/R_S**2 + 1)**2)
        return [f, f_R, f_RR, f_RRR]
    def potential(self,R):
        beta = self.R_tilde
        q = self.q
        llambda =  self.llambda
        
        if q==1: 
            V1 = 1/3*(R/2* (R - 4*llambda - 2*llambda*(1 + R**2)**(-1) ) + 3*llambda*np.arctan(R) )
        elif q ==2:
            V1 = 1/6*(R**2 - llambda*R*R_S* (4*R**4 + 5*R**2*R_S**2 + 3*R_S**4)/(R**2 + R_S**2)**2 ) + llambda *R_S**2/2*np.arctan(R/R_S)
        return V1

class MJW(ModelFR):
    """
    Class for the Miranda f(R) Model as defined in:
    
        Miranda, V. et al. “Viable Singularity-Free f(R) Gravity Without a 
        Cosmological Constant”. Physical Review Letters, 102, (2009), 221101
        
    
    PARAMETERS
    ----------
    
    R_tilde : float
        Dimensionless parameter.
    alpha_M : float
        Dimensionless parameter
    """

    def __init__(self, R_tilde, alpha_M):
        self.R_tilde= R_tilde
        self.alpha_M =  alpha_M

    def __call__(self,R):
        R_tilde, alpha_M = self.R_tilde,self.alpha_M 
        f = R - alpha_M*R_tilde*np.log(1 + R_tilde**(-1)*R)
        f_R = 1 - alpha_M*(1 + R_tilde**(-1)*R)**(-1)
        f_RR = alpha_M/R_tilde*(1 + R_tilde**(-1)*R)**(-2)
        f_RRR = 2*alpha_M/R_tilde**2*(1 + R_tilde**(-1)*R)**(-3)
        return [f, f_R, f_RR, f_RRR]
    
    def potential(self,R):
        R_tilde= self.R_tilde
        alpha_M = self. alpha_M
        V = 1/3*((1 + R_tilde*R)*(R_tilde*R+ 6* alpha_M-1 ) - 2* alpha_M*(3 + 2*R*R_tilde)*np.log(1 + R_tilde*R))
        return V

class HUSAW(ModelFR):

    """
    Class for the Hu & Sawicki f(R) Model as defined in:

         Hu, W. and Sawicki, I. “Models of f(R) Cosmic Acceleration 
         That Evade Solar-System Tests”. Physical Review D, 76, (2007), 064004.


    PARAMETERS
    ----------
    
    m2 : float
        Dimensionless parameter.
    c_1 : float
        Dimensionless parameter.
    c_2 : float
        Dimensionless parameter.
    n : integer
        Dimensionless parameter  
    """

    def __init__(self,m2,c_1,c_2,n):
        self.m2= m2
        self.c_1 =  c_1
        self.c_2 =  c_2
        self.n = n
    
    def __call__(self, R):
        m2, c_1, c_2, n = self.m2, self.c_1, self.c_2, self.n
        
        f = R - c_1*m2*(R/m2)**n/(c_2*(R/m2)**n + 1)
        f_R = 1 + c_1*c_2*m2*n*(R/m2)**(2*n)/(R*(c_2*(R/m2)**n + 1)**2) - c_1*m2*n*(R/m2)**n/(R*(c_2*(R/m2)**n + 1))
        f_RR = c_1*m2*n*(R/m2)**n*(-2*c_2**2*n*(R/m2)**(2*n)/(c_2*(R/m2)**n + 1)**2 + 3*c_2*n*(R/m2)**n/(c_2*(R/m2)**n + 1) - c_2*(R/m2)**n/(c_2*(R/m2)**n + 1) - n + 1)/(R**2*(c_2*(R/m2)**n + 1))
        f_RRR = c_1*m2*n*(R/m2)**n*(6*c_2**3*n**2*(R/m2)**(3*n)/(c_2*(R/m2)**n + 1)**3 - 12*c_2**2*n**2*(R/m2)**(2*n)/(c_2*(R/m2)**n + 1)**2 + 6*c_2**2*n*(R/m2)**(2*n)/(c_2*(R/m2)**n + 1)**2 + 7*c_2*n**2*(R/m2)**n/(c_2*(R/m2)**n + 1) - 9*c_2*n*(R/m2)**n/(c_2*(R/m2)**n + 1) + 2*c_2*(R/m2)**n/(c_2*(R/m2)**n + 1) - n**2 + 3*n - 2)/(R**3*(c_2*(R/m2)**n + 1))
        return [f, f_R, f_RR, f_RRR]
    
    def potential(self, R):
        R_tilde= self.R_tilde
        alpha_M = self. alpha_M
        V = 1/3*((1 + R_tilde*R)*(R_tilde*R+ 6* alpha_M-1 ) - 2* alpha_M*(3 + 2*R*R_tilde)*np.log(1 + R_tilde*R))
        return V
    
