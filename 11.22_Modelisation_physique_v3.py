# Modelisation physique de la grue flottante
# Groupe 11-22

import math as m

############
# Constantes
############

g = 9.81

############################################
# Dimensions  et poids de la base de la grue
############################################

long_base = 0.6                    # [m]
larg_base = 0.6                    # [m]
h_base = 0.06                      # [m]
m_vol_base = 240                   # [kg/m**3] Cette valeur est inférieure à celle du bois, car la base n'est en réalité pas pleine
                                   # Afin de simplifier les calculs, elle a donc été ajustée pour correspondre à son poids dans la vraie vie

poids_base = long_base*larg_base*h_base*m_vol_base*g
m_base = (poids_base / g)

##############################################
# Dimensions et poids du bras et du grappin
##############################################

long_bras = [0.3, 0.2575, 0.2125, 0.1475]     # [m]   
m_vol_bras = 415                              # [kg/m**3] Valeur réelle de la densité du bois, calculée pour correspondre aux valeurs dans la vraie vie
larg_bras = [0.038, 0.028, 0.028, 0.023]      # [m]               
prof_bras = [0.038, 0.028, 0.021, 0.023]      # [m]        

poids_bras = []
for i in range(len(long_bras)):
    poids_bras.append((long_bras[i]*larg_bras[i]*prof_bras[i])*m_vol_bras*g)
    
m_bras = []
for i in range(len(poids_bras)):
    m_bras.append(poids_bras[i] / g) # On a ainsi les valeurs des poids des bras individuels, ainsi qu'une valeur qui regroupe l'ensemble des bras

poids_grappin = 0.086 * g                     # [N]

poids_grue = 0
for i in range(len(poids_bras)):
    poids_grue += poids_bras[i]

poids_grue += poids_grappin + 0.4 * 9.81 # <- poids auxiliaire du bras, càd des seringues, du système de rotation, etc. Estimé via Fusion 360
m_grue = poids_grue / g


###################
# Charge a soulever
###################

m_charge = 0.2                  # [kg]
poids_charge = g * m_charge
d_charge = 0.7                  # [m]
h_charge = 0.623                # [m]

position_charge = (d_charge, h_charge)

########################################################
# Calcul du centre de gravité dans différentes positions
########################################################

# Ces valeurs ont été mesurées via Fusion 360, G est le centre de gravité du bras, du grappin et de la charge

G = ((0.06996,0.17589,"1er et/ou 2eme fut ramasse"),(0.20851,0.16151,"1er fut place"),(0.20328,0.2341,"2eme fut place"),(0.09579,0.18272,"3eme fut ramasse"),(0.2635, 0.34154, "3eme fut place"))
Gbase = (-0.0022,0.0053) 
Gtotal = ((G[4][0]*(m_grue+m_charge)+Gbase[0]*m_base)/(m_grue+m_charge+m_base),(G[4][1]*(m_grue+m_charge)+Gbase[1]*m_base)/(m_grue+m_charge+m_base))
# Il s'avère en réalité qu'on n'a simulé que la dernière situation dans notre programme
# Gtotal est donc le centre de gravité de tout le prototype lorsqu'il tient le 3eme fut à la hauteur maximale

################################################
# Calcul de la hauteur de la ligne de flottaison
################################################

h_c = (poids_grue + poids_base + poids_charge)/((4 * 0.095) * long_base * 1000 * g)


#########################################
# Calcul de l'angle d'inclinaison maximal
#########################################

theta_sub =  - m.atan((h_base-h_c)/(long_base/2)) # Angle maximal avant la submersion de la base
theta_soul =  - m.atan((h_c)/(long_base/2)) # Angle maximal avant le soulèvement du fond de la base
global theta_max
theta_max = max(theta_sub,theta_soul) # L'angle à ne pas dépasser est le plus petit (en valeur absolue) des deux angles critiques


####################################################################################
# Calcul de la position du centre de poussée après une perturbation d'un angle theta
####################################################################################

def centre_de_poussee(theta):
    
    """Retourne la position en x du centre de poussee lorsque la base est desequilibree d'un angle theta"""
    
    # Note : on fait ici l'hypothèse que, bien que le volume immergé n'est pas un trapèze, la position de son centre de gravité est la même 
        
    Base = h_c + abs((long_base/2) * m.tan(theta))
    base = h_c - abs((long_base/2) * m.tan(theta))
    global vol_trapeze
    vol_trapeze = h_c * 0.095 * 4 * long_base # Comme précisé plus haut, la base n'est en réalité pas pleine, le calcul du vrai volume immergé est donc celui-ci 
    
    # Position en x du centre de poussee
    
    global X_c
    X_c = (long_base/2) - (long_base/3) * ((Base+2*base)/(Base+base)) # Formule du centre de gravité du trapèze décrite dans la modélisation physique
    
    return X_c


########################################################################
# Calcul des nouvelles positions d'un couple après une rotation de theta
########################################################################

def rotation(position, theta) :
    
    """Retourne la position en x et en y d'une position de depart apres une rotation d'un angle theta
    'position' est un tuple (x,y) et theta l'angle entre la position post-rotation, l'origine de la rotation et l'horizontale"""
    
    x,y = position
    r = m.sqrt(x**2+y**2)
    
    # on calcule l'angle que forme le point de base avec l'horizontale
    if x == 0 : # Si x vaut 0, il faut éviter qu'il y ait une ZeroDivisionError
        angle = m.pi/2
    else :
        angle = m.atan(y/x) 
    return (r*m.cos(angle + theta),r*m.sin(angle + theta))


############################################################
# Calcul de la somme des couples après une rotation de theta
############################################################

def somme_des_couples(theta,Ca):
    
    "Retourne la somme totale des couples, puis Ca, Cg et Cr, le tout est un tuple"
    
    centre_de_poussee(theta) # Cette fonction globalise certaines valeurs, il faut donc l'éxécuter pour mettre à jour les valeurs du centre de poussée
    
    Cg = - rotation(Gbase, theta)[0] * poids_base # Couple généré par la base
    Cr = (centre_de_poussee(theta) - rotation(Gtotal, theta)[0]) * vol_trapeze * 1000 * g # Couple de redressement généré par le flotteur
    # Note : On effectue ici une rotation de Gtotal et de Gbase, car ils ont changé légèrement de position lors de l'inclinaison de la grue
    
    somme = (Ca + Cg + Cr) 
    return (somme, Ca, Cg, Cr)

def modele_1():
    """ Recherche dichotomique de l'angle d'équilibre, càd lorsque Ca + Cg + Cr == 0 """
    first = 0
    last = theta_max
    while first >= last :
        theta = (first + last)/2
        
        Ca = - (poids_grue+poids_charge) * rotation((G[4][0],G[4][1]), theta)[0] # A nouveau, seul la derniere situation nous intéresse
        Cg = somme_des_couples(theta, Ca)[2]
        Cr = somme_des_couples(theta, Ca)[3]
        
        if abs(somme_des_couples(theta, Ca)[0]) < 0.0000001 :
            break
        elif abs(Cr) < abs(Ca) :
            first = theta - (1/10000000)
        else:
            last = theta + (1/10000000)
    return theta * 180 / m.pi
  
def modele_2():
    """ Recherche dichotomique de l'angle d'équilibre, càd lorsque le centre de gravité est au dessus du centre de poussée """
    first = 0
    last = theta_max
    while first >= last :
        theta = (first + last) / 2
        
        if abs(rotation(Gtotal, theta)[0] - centre_de_poussee(theta)) < 0.000001 :
            break 
        elif rotation(Gtotal, theta)[0] < centre_de_poussee(theta) :
            last = theta + (1/1000000)
        else :
            first = theta - (1/1000000)
    return (theta / m.pi) * 180 
                      
if __name__ == "__main__" :           
    print("Modèle 1 : theta = ", modele_1())
    print("Modèle 2 : theta = ", modele_2())
