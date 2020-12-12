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
m_vol_base = 240                   # [kg/m**3]

poids_base = long_base*larg_base*h_base*m_vol_base*g
m_base = (poids_base / g)

##############################################
# Dimensions et poids de la grue et du grappin
##############################################

long_bras = [0.3, 0.2575, 0.2125, 0.1475]     # [m]   
m_vol_bras = 415                              # [kg/m**3] 
larg_bras = [0.038, 0.028, 0.028, 0.023]      # [m]               
prof_bras = [0.038, 0.028, 0.021, 0.023]      # [m]        

poids_bras = []
for i in range(len(long_bras)):
    poids_bras.append((long_bras[i]*larg_bras[i]*prof_bras[i])*m_vol_bras*g)
    
m_bras = []
for i in range(len(poids_bras)):
    m_bras.append(poids_bras[i] / g)

poids_grappin = 0.086 * g                     # [N]

poids_grue = 0
for i in range(len(poids_bras)):
    poids_grue += poids_bras[i]

poids_grue += poids_grappin + 0.4 * 9.81 # <- poids du support du bras
m_grue = poids_grue / g


###################
# Charge a soulever
###################

m_charge = 0.2                  # [kg]
poids_charge = g * m_charge
d_charge = 0.7                  # [m]
h_charge = 0.623                # [m]
charges = [0.2,0.14,0.07]
hCharges = [0.623,0.25,0.45]

position_charge = (d_charge, h_charge)

'''
##################################################################
# Liste des angles des différents bras par rapport à l'horizontale
##################################################################

liste_des_inclinaisons = [(13 * m.pi)/30, (11 * m.pi)/60, (4 * m.pi)/45, m.pi/180]


########################################################
# Calcul du centre de gravite de chaque bras sans charge
########################################################

def G_x_bras(liste_des_inclinaisons) :                                                 # Par rapport a l'axe parallèle a la base
    
    "Retourne une liste des positions en x du centre de gravite de chaque bras"
    
    global G_x_par_bras
    G_x_par_bras = []
    
    G_x_par_bras.append((m.cos(liste_des_inclinaisons[0])*((long_bras[0])/2)))         # Valeur pour le premier bras
    for i in range(1, len(liste_des_inclinaisons)) :
        G_x_par_bras.append((G_x_par_bras[-1])+(m.cos(liste_des_inclinaisons[i-1])*((long_bras[i-1])/2))+(m.cos(liste_des_inclinaisons[i])*((long_bras[i])/2)))
    
    return G_x_par_bras

def G_y_bras(liste_des_inclinaisons) :
    
    "Retourne une liste des positions en y du centre de gravite de chaque bras"
    
    global G_y_par_bras
    G_y_par_bras = []
    
    G_y_par_bras.append((m.sin(liste_des_inclinaisons[0])*((long_bras[0])/2) + (h_base / 2)))
    for i in range(1, len(liste_des_inclinaisons)) :
        G_y_par_bras.append((G_y_par_bras[-1])+(m.sin(liste_des_inclinaisons[i-1])*((long_bras[i-1])/2))+(m.sin(liste_des_inclinaisons[i])*((long_bras[i])/2)))
        
    return G_y_par_bras


########################################################################
# Calcul de la position du centre de gravité de tous les bras de la grue
########################################################################

def G_x_grue():
    
    "Retourne la position en x du centre de gravite des bras de la grue (sans masse et sans contrepoids)"
    
    G_x_bras(liste_des_inclinaisons)        # on exécute cette fonction pour avoir accès à la liste G_x_par_bras
    
    somme_xi_mi = 0
    G_x_bras(liste_des_inclinaisons)
    for i in range(len(G_x_par_bras)):
        somme_xi_mi += (G_x_par_bras[i])*(m_bras[i])
        
    somme_mi = sum(m_bras) + m_base
    
    G_x_grue = somme_xi_mi / somme_mi
    
    return G_x_grue

def G_y_grue():
    
    "Retourne la position en y du centre de gravite des bras de la grue (sans masse et sans contrepoids)"
    
    G_y_bras(liste_des_inclinaisons)
    
    somme_yi_mi = 0
    for i in range(len(G_y_par_bras)):
        somme_yi_mi += (G_y_par_bras[i])*(m_bras[i])

    somme_mi = sum(m_bras) + m_base

    G_y_grue = somme_yi_mi / somme_mi
    
    return G_y_grue


G_grue = (G_x_grue(), G_y_grue())


#######################################################################
# Calcul du centre de gravité extrême (grue chargée à 40 cm de la base)
#######################################################################

def G_x_charge():
    
    "Retourne la position en x du centre de gravite de la grue dans le premier cas"

    G_x_bras(liste_des_inclinaisons)     # Mesures sur Fusion 360, version du 15/11/2020
    G_x_charge = ((G_x_grue() * m_grue) + (d_charge * m_charge))/(m_grue + m_charge + m_base)
    return G_x_charge

def G_y_charge():
    
    "Retourne la position en x du centre de gravite de la grue dans le premier cas"
    
    G_y_bras(liste_des_inclinaisons)
    G_y_charge = ((G_y_grue() * m_grue) + (h_charge * m_charge))/(m_grue + m_charge + m_base)
    return G_y_charge


G_chargé = (G_x_charge(), G_y_charge())
'''
G = ((0.06996,0.17589,"1er et/ou 2eme fut ramasse"),(0.20851,0.16151,"1er fut place"),(0.20328,0.2341,"2eme fut place"),(0.09579,0.18272,"3eme fut ramasse"),(0.2635, 0.34154, "3eme fut place"))
Gbase = (-0.0022,0.0053)
Gtotal = ((G[4][0]*(m_grue+m_charge)+Gbase[0]*m_base)/(m_grue+m_charge+m_base),(G[4][1]*(m_grue+m_charge)+Gbase[1]*m_base)/(m_grue+m_charge+m_base))

###############
# Calcul de h_c
###############

h_c = (poids_grue + poids_base + poids_charge)/((4 * 0.095) * long_base * 1000 * g)           # Volume immergé


#########################################
# Calcul de l'angle d'inclinaison maximal
#########################################

theta_sub =  m.atan((h_base-h_c)/(long_base/2))
theta_soul =  m.atan((h_c)/(long_base/2))
global theta_max
theta_max = min(theta_sub,theta_soul)


####################################################################################
# Calcul de la position du centre de poussée après une perturbation d'un angle theta
####################################################################################

def centre_de_poussee(theta):
    
    "Retourne la position en x du centre de poussee lorsque la base est desequilibree d'un angle theta"
    
    # Dimensions du trapeze lorsque la base est desequilibrée
        
    Base = h_c + abs((long_base/2) * m.tan(theta))
    base = h_c - abs((long_base/2) * m.tan(theta))
    global vol_trapeze
    #vol_trapeze = (((Base + base) * long_base)/2) * larg_base
    vol_trapeze = h_c * 0.095 * 4
    
    # Position en x du centre de poussee
    
    global X_c
    X_c = (long_base/2) - (long_base/3) * ((Base+2*base)/(Base+base))
    
    return X_c


########################################################################
# Calcul des nouvelles positions d'un couple après une rotation de theta
########################################################################

def rotation(position, theta) :
    
    """Retourne la position en x et en y d'une position de depart apres une rotation d'un angle theta
    'position' est un tuple (x,y) et theta l'angle entre la position post-rotation, l'origine de la rotation et l'horizontale"""
    
    x,y = position
    r = m.sqrt(x**2+y**2)
    
    # on calcule l'angle que forme le point de base avec la verticale
    if x == 0 :
        angle = m.pi/2
    else :
        angle = m.atan(y/x)
    return (r*m.cos(angle + theta),r*m.sin(angle + theta))


############################################################
# Calcul de la somme des couples après une rotation de theta
############################################################

def somme_des_couples(theta, Ca):
    
    "Retourne la somme totale des couples, puis Ca, Cg et Cr, le tout est un tuple"
    
    centre_de_poussee(theta)
    
    # Ca = - poids_grue * rotation((G[4][0],G[4][1]), theta)[0]
    Cg = - rotation(Gtotal, theta)[0] * (poids_grue +poids_charge + poids_base)
    Cr = (centre_de_poussee(theta) - rotation((Gtotal), theta)[0]) * vol_trapeze * 1000 * g
    
    somme = (Cr + Cg)
    return (somme, Ca, Cg, Cr)

def modele_1():
    first = 0
    last = theta_max
    while first <= last :
        theta = (first + last)/2
        
        Ca = - (poids_grue+poids_charge) * rotation((G[4][0],G[4][1]), theta)[0]
        
        Ca = somme_des_couples(theta, Ca)[1]
        Cg = somme_des_couples(theta, Ca)[2]
        Cr = somme_des_couples(theta, Ca)[3]
        
        if abs(somme_des_couples(theta, Ca)[0]) < 0.001 :
            return (theta / m.pi) * 180
        elif abs(Cr) < abs(Cg) :
            last = theta - (1/1000)
            print(theta)
        else:
            first = theta + (1/1000)
            print(theta)
                
def modele_2():
    first = 0
    last = theta_max
    while first <= last :
        theta = (first + last) / 2
        
        if abs(rotation(Gtotal, theta)[0] - centre_de_poussee(theta)) < 0.000001 :
            return (theta / m.pi) * 180 
        elif rotation(Gtotal, theta)[0] < centre_de_poussee(theta) :
            last = theta - (1/1000000)
            #print(theta)
        else :
            first = theta + (1/1000000)
            #print(theta)
                
                      
if __name__ == "__main__" :           
    print(modele_1())
    print(modele_2())
    
    
#def somme_des_couples(theta, Ca):
    
""""Retourne la somme totale des couples, puis Ca, Cg et Cr, le tout est un tuple"
    
    centre_de_poussee(theta)
    
    Ca = - poids_charge * rotation(position_charge, theta)[0]
    Cg = - rotation(G_grue, theta)[0] * (poids_grue)
    Cr = (X_c - rotation(G_chargé, theta)[0]) * vol_trapeze * 1000 * g
    
    somme = (Cr + Cg + Ca)
    return (somme, Ca, Cg, Cr)"""
