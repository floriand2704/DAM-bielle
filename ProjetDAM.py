import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


D = 0.076 #[m]
C = 0.07 #[m]
L = 0.0122 #[m]
mpiston = 0.353#[kg]
mbielle = 0.549 #[kg]
tau = 10 #[-]
Mair_carb_essence =  14.5#[kg_air/kg_fuel] pour un moteur essence ou 
Mair_carb_diesel = 26 #[kg_air/kg_fuel] pour un moteur diesel  @ #[kg_air/kg_fuel]
R = C/2
beta = L/R

def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    #VOTRE CODE
    thetaRad = theta*(np.pi/180) 
    thetaCRad = -1*thetaC*(np.pi/180) 
    deltaThetaCRad = deltaThetaC*(np.pi/180)  


    V_output = vol ( thetaRad)
    Q_output = chaleur( thetaRad, thetaCRad, deltaThetaCRad)
    F_pied_output = force_p_o (rpm,s, thetaRad, thetaCRad, deltaThetaCRad)
    F_tete_output = force_t_o (rpm,s, thetaRad, thetaCRad, deltaThetaCRad)
    p_output = pression (s, thetaRad, thetaCRad, deltaThetaCRad)
    t = epaisseur (rpm,s, thetaRad, thetaCRad, deltaThetaCRad)

    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)

# normallement correct
def vol (theta):
    Vmax = ((np.pi*D**2)/2)*R
    Vmin= 0
    Vc = Vmax - Vmin
    
    V_output = (Vc/2) * ( 1- np.cos(theta)+ beta - np.sqrt(beta**2 - np.sin(theta)**2)) + 1/(tau-1)*Vc
    return V_output
# manque Q_tot
def chaleur( theta, thetaC, deltaThetaC):


    # je sais pas comment trouver Q_tot
    Q_tot = 0
    # Q que quand combustion
    Q_output = Q_tot *0.5*(1- np.cos(np.pi*(theta - thetaC)/deltaThetaC))

    return Q_output

# on considère que la bielle a une forme de section en "I"
# normallement correct 
def force_p_o (rpm,s, theta, thetaC, deltaThetaC):

    P_flattened = pression(rpm,s, theta, thetaC, deltaThetaC)
    w = rpm *2 *np.pi/60
    F_pied_output = (np.pi * (D**2)*P_flattened/4) - mpiston*R*(w**2)*np.cos(theta)                  #force sur le pied de la bielle en fonction de l'angle de rotation [N]
   

    return F_pied_output
# normallement correct 
def force_t_o (rpm,s, theta, thetaC, deltaThetaC):

    P_flattened = pression(rpm,s, theta, thetaC, deltaThetaC)
    w = rpm *2 *np.pi/60
    F_tete_output = -(np.pi * (D**2)*P_flattened/4) + (mpiston + mbielle)*R*(w**2)*np.cos(theta)     #force sur la tête de la bielle en fonction de l'angle de rotation [N]
    
    return F_tete_output

# il manque Qtot
def pression (s, theta, thetaC, deltaThetaC):

    gamma = 1,3

    def EDO (p,theta):
        Vmax = ((np.pi*D**2)/2)*R
        Vmin= 0
        Vc = Vmax - Vmin
        dV_dthet = (Vc/2) * (np.sin(theta) + (np.sin(theta*2)/np.sqrt((beta**2) - np.sin(theta)**2)))
        # apport à lieu que lors de la phase de combustion ( je sais pas quand c'est)
        # je sais pas comment trouver Q_tot
        dQ_dthet = Q_tot * 0.5 * np.pi * ((theta-thetaC)/deltaThetaC) * np.sin(np.pi * (theta - thetaC)/deltaThetaC)
        return -gamma * (p/vol(theta))*(dV_dthet)+((gamma - 1)/vol(theta))*dQ_dthet
    p_admission = s * 10**5
    pres = odeint(EDO, p_admission, theta)
    p_flattend = np.ravel (pres)

    return p_flattend
#vérif si max ou min
def epaisseur (rpm,s, theta, thetaC, deltaThetaC):

    E = 200 * 10**9
   
    I_x = 419/12
    I_y = 131/12
    
    K_x = 1
    K_y = 0.5
    
    F_euler_x = ((np.pi **2) * E *I_x)/(K_x*L)**2
    F_euler_y = ((np.pi **2) * E *I_y)/(K_y*L)**2
    
    sigam_c = 450*10**6 #[MPa]
    
    # F max dans le sens de la compression
    F_p_o_max = np.max(force_p_o(rpm,s, theta, thetaC, deltaThetaC))
    F_t_o_max = np.min(force_t_o(rpm,s, theta, thetaC, deltaThetaC))
    F_critique = np.max(-F_t_o_max, F_p_o_max)
   
    t_x = np.roots(F_critique, 0, -1/(11*sigam_c), 0, F_euler_x)
    t_y = np.roots(F_critique, 0, -1/(11*sigam_c), 0, F_euler_y)
    t_x1 = np.array([])
    t_y1 = np.array([])
    for i in range (t_x.size):
        if (t_x[i].imag == 0):
            np.append(t_x1, t_x[i])
        if (t_y[i].imag == 0):
            np.append(t_y1, t_y[i])    
    # ici je sais pas si c'est max ou min !!!!
    t_xmax = np.max(t_x1) if t_x1.size > 0 else 0
    t_ymax = np.max(t_y1) if t_y1.size > 0 else 0

    t = np.max (t_xmax, t_ymax)
    
    return t
