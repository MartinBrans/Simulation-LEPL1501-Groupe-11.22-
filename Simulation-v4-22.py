# Simulation d'une charge accrochée à la grue flottante

import math
import matplotlib.pyplot as plt
import numpy as np
import Modelisation_physique_v3 as mod

### Paramètres du système

D = 0.6              # Coefficient d'amortissement
I = 0.12            # Moment d'inertie

### Paramètres de la simulation

step = 0.001         # pas (dt) [s]
end = 10.0           # durée [s]
theta_0 = 0.0        # position initiale [m]
omega_0 = 0.0        # vitesse initiale [m/s]
Ca_0 = - mod.G[3][0] * (mod.poids_grue + mod.poids_charge)          # couple initial [N m]
temps = 3            # temps que met la charge à arriver à sa position finale [s]

t = np.arange(0, end, step)
theta = np.empty_like(t)          
omega = np.empty_like(t)
a = np.empty_like(t)
Gx = np.empty_like(t)
Gy = np.empty_like(t)
Ca = np.empty_like(t)

E_G = np.empty_like(t)
E_C = np.empty_like(t)
E_K = np.empty_like(t)
E_A = np.empty_like(t)
E_Total = np.empty_like(t)

def simulation():
    """
    pre: -
    post: exécute une simulation jusqu'à t = end par pas de dt = step.
          Remplit les listes theta, omega et a des positions, vitesses et accélérations.
    """
    # conditions initiales
    
    theta[0] = theta_0
    omega[0] = omega_0
    Ca[0] = Ca_0
    Gx[0] = mod.G[3][0]
    Gy[0] = mod.G[3][1]
    
    for i in range(len(t)-1):
        
        dt = step
        
        """Eqt de la droite reliant la position initiale a la position finale de G dans la situation 3 en fct du temps :
            Gx[i+1] = 0.055903333333333333*(i+1)*dt + 0.09579
            Gy[i+1] = 0.05294*(i+1)*dt + 0.18272
        """
        
        if i * dt < temps:
            Gx[i+1] = 0.055903333333333333*(i+1)*dt + 0.09579
            Gy[i+1] = 0.05294*(i+1)*dt + 0.18272
        else :
            Gx[i+1] = mod.G[4][0]
            Gy[i+1] = mod.G[4][1]
        
        Ca[i+1] = - Gx[i+1] * (mod.poids_grue + mod.poids_charge)
        
        C = mod.somme_des_couples(theta[i], Ca[i])[0] - D*omega[i]
        
        theta[i+1] = theta[i] + omega[i] * dt
        omega[i+1] = omega[i] + a[i] * dt
        a[i+1] = C / I
        
    ### Calculs d'énergie
        
        E_G[i+1] = (mod.poids_grue + mod.poids_charge + mod.poids_base)*(mod.rotation(mod.Gtotal,math.atan(mod.Gtotal[1]/mod.Gtotal[0])+theta[i])[1]-mod.Gtotal[1])
        #E_C[i+1] = -(mod.poids_ensemble + mod.poids_charge)*(y_C-y_C[0])
        E_K[i+1] = I*omega[i]*omega[i]/2
        E_A[i+1] = -mod.somme_des_couples(theta[i],Ca[i])[1] * theta[i]
        E_Total[i+1] = E_K[i+1] + E_G[i+1] + E_A[i+1] # + E_C[i+1]
        
def distance_couple() :
    
    plt.figure(1)
    plt.subplot(1,2,1)
    plt.plot(t, Gx, label = "Position horizontale du centre de gravite")
    plt.xlabel("Temps [s]")
    plt.ylabel("distance [m]")
    plt.legend()
    plt.subplot(1,2,2)
    plt.plot(t, Ca, label = "Couple applique")
    plt.xlabel("Temps [s]")
    plt.ylabel("Couple [N m]")
    plt.legend()
    plt.show()
                
def graphiques():
    
    plt.figure(1)
    plt.subplot(3,1,1)
    plt.plot(t,theta*180/math.pi, label="Angle")
    plt.xlabel("Temps [s]")
    plt.ylabel("Theta [rad]")
    plt.legend(loc="upper right")
    plt.subplot(3,1,2)
    plt.plot(t,omega, label="Vitesse angulaire")
    plt.xlabel("Temps [s]")
    plt.ylabel("Omega [rad/s]")
    plt.legend()
    plt.subplot(3,1,3)
    plt.plot(t,a, label="Acceleration angulaire")
    plt.xlabel("Temps [s]")
    plt.ylabel("a [rad/s**2]")
    plt.legend()
    plt.show()
    
def graphique_theta() :
    
    plt.figure(1)
    plt.plot(t, theta, label = "Theta")
    plt.plot([0,end], [mod.theta_sub, mod.theta_sub], '--r', label='Submersion')
    plt.plot([0,end], [-mod.theta_sub, -mod.theta_sub], '--r')
    plt.plot([0,end], [mod.theta_soul, mod.theta_soul], '--m', label='Soulevement')
    plt.plot([0,end], [-mod.theta_soul, -mod.theta_soul], '--m')
    plt.legend()
    plt.show()
    
def diagramme_de_phase() :
    
    plt.figure(1)
    plt.title("Diagramme de phase")
    plt.plot(theta, omega)
    plt.xlabel("Theta [rad]")
    plt.ylabel("Omega [rad/s]")
    plt.show()
    
def graphique_energies() :
    
    plt.figure(1)
    plt.title("Energies en fonction du temps")
    plt.plot(t, E_G,label="En. pot. gravifique")
    plt.plot(t, E_K,label="En. cinetique")
    plt.plot(t, E_Total,label="En. total du systeme")
    plt.xlabel("Temps [s]")
    plt.ylabel("Energie [J]")
    plt.legend(loc="upper right")
    plt.show()
    
simulation()
distance_couple()
graphiques()
graphique_theta()
diagramme_de_phase()
graphique_energies()