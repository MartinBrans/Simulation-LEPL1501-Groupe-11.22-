# Simulation d'une charge accrochée à la grue flottante et qui se déplace pendant 3 secondes vers l'avant

import math
import matplotlib.pyplot as plt
import numpy as np
import Modelisation_physique_v3 as mod

### Paramètres du système

D = 0.5             # Coefficient d'amortissement, cette valeur a été estimée, mais ne peut être confirmée que expérimentalement
I = 0.12            # Moment d'inertie, calculé sur Fusion 360


### Paramètres de la simulation

step = 0.001         # pas (dt) [s]
end = 10.0           # durée [s]
theta_0 = 0.0        # position initiale [m]
omega_0 = 0.0        # vitesse initiale [m/s]
Ca_0 = - mod.G[3][0] * (mod.poids_grue + mod.poids_charge)          # couple initial [N m]
temps = 3            # temps que met la charge à arriver à sa position finale [s]

t = np.arange(0, end, step) # Temps
theta = np.empty_like(t)    # Angle d'inclinaison      
omega = np.empty_like(t)    # Vitesse angulaire
a = np.empty_like(t)        # Accélération angulaire
Gx = np.empty_like(t)       # Position du centre de gravité du bras + de la charge (composantes horizontale)
Gy = np.empty_like(t)       # Position du centre de gravité du bras + de la charge (composantes verticale)
Ca = np.empty_like(t)       # Couple d'affaissement

E_G = np.empty_like(t)      # Energie potentielle gravifique
E_C = np.empty_like(t)      # Energie potentielle de la poussée d'archimède
E_K = np.empty_like(t)      # Energie cinétique
E_A = np.empty_like(t)      # Engergie potentielle du couple d'affaissement
E_Total = np.empty_like(t)  # Energie totale du système

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
        

        
        if i * dt < temps:
            # Hypothèse : le centre de gravité du bras + de la charge se déplace en ligne droite
            Gx[i+1] = ((mod.G[4][0]-mod.G[3][0])/temps)*(i+1)*dt + mod.G[3][0] # x = mt + p avec m = delta x/delta t
            Gy[i+1] = ((mod.G[4][1]-mod.G[3][1])/temps)*(i+1)*dt + mod.G[3][1] # y = mt + p avec m = delta y/delta t
        else :
            Gx[i+1] = mod.G[4][0] # Lorsque la masse a fini son mouvement, sa position est constante
            Gy[i+1] = mod.G[4][1]
        
        Ca[i+1] = - Gx[i+1] * (mod.poids_grue + mod.poids_charge) # Couple d'affaissement généré par le bras et la masse
        
        C = mod.somme_des_couples(theta[i],Ca[i])[0] - D*omega[i] # Somme des couples = Cr + Cg + Ca + Cd
        
        theta[i+1] = theta[i] + omega[i] * dt # On utilise des approximations linéaires à chaque pas pour calculer la position du point suivant
        omega[i+1] = omega[i] + a[i] * dt
        a[i+1] = C / I
        
    ### Calculs d'énergie
        
        E_G[i+1] = mod.poids_base*(mod.rotation(mod.Gbase,math.atan(mod.Gbase[1]/mod.Gbase[0])+theta[i])[1]-mod.Gbase[1])
        #E_C[i+1] = -(mod.poids_ensemble + mod.poids_charge)*(y_C-y_C[0]) <- cette formule n'a pas été intégrée, pour le moment
        #                                                                    car aucune fonction ne renvoie la composante verticale de la position du centre de poussée
        E_K[i+1] = I*omega[i]*omega[i]/2
        E_A[i+1] = Ca[i] * theta[i]
        E_Total[i+1] = E_K[i+1] + E_G[i+1] + E_A[i+1] # + E_C[i+1]
        
def distance_couple() :
    """ Cette fonction sert à s'assurer que la position de la charge varie bien au cours du temps"""
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
    plt.plot(t,theta, label="Angle")
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
    
def graphique_detaille_theta() :
    
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
    
def graphiques_energies() :
    
    plt.figure(1)
    plt.title("Energies en fonction du temps")
    plt.plot(t, E_G,label="En. pot. gravifique")
    plt.plot(t, E_K,label="En. cinetique")
    plt.plot(t, E_A,label="En. pot. couple d'affaissement")
    plt.plot(t, E_Total,label="En. total du systeme")
    plt.xlabel("Temps [s]")
    plt.ylabel("Energie [J]")
    plt.legend(loc="upper right")
    plt.show()

if __name__ == "__main__" :
    simulation()
    # distance_couple()
    graphiques()
    graphique_detaille_theta()
    diagramme_de_phase()
    graphiques_energies()
