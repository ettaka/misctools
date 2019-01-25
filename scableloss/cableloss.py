import math
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

mu0 = 4e-7 * math.pi
"""Jc Nb3Sn"""

NB3SNC_matparam = {}
NB3SNC_matparam['C'] = 4.3e10
NB3SNC_matparam['Bc20'] = 27.012
NB3SNC_matparam['Tc0'] = 18

HFMD07B_strandparam = {}
HFMD07B_strandparam['a'] = 0.7e-3/2.
HFMD07B_strandparam['RRR'] = 100.
HFMD07B_strandparam['lf'] = 14e-3
HFMD07B_strandparam['cu_to_noncu'] = 1.106
#Bc2(T) = Bc20*(1-(T/Tc0)**2)*(1-0.31*(T/Tc0)**2.)*(1-1.77*math.log(T/Tc0))
#Jc(B,T) = C/sqrt(B)*(1-B/Bc2(T))**2.*(1-(T/Tc0)**2)**2

def get_Bc2(T,matparam):
    Bc20 = matparam['Bc20']
    Tc0 = matparam['Tc0']
    return Bc20*(1-(T/Tc0)**2)*(1-0.31*(T/Tc0)**2.*(1-1.77*np.log(T/Tc0)))

def get_Jc(B,T,matparam):
    Tc0 = matparam['Tc0']
    C = matparam['C']
    return C/np.sqrt(B)*(1-B/get_Bc2(T,matparam))**2.*(1-(T/Tc0)**2)**2

def test_summers():
    matparam = {}
    matparam['C'] = 4.3e10
    matparam['Bc20'] = 27.012
    matparam['Tc0'] = 18
    T = 4.2
    B = 12
    print "Bc2:", get_Bc2(T, matparam)
    print "Jc:",  get_Jc(B, T, matparam)

""" Hysterisis loss [Wilson]"""
def get_round_filament_loss_per_cycle_par(Bm, Jc_theta, a):
    beta = Bm/(2*mu0*Jc_theta*a)
    loss = Bm**2./(2.*mu0)
    if beta <= 1:
        loss *= (2*beta/3-beta**2/3)
    elif beta > 1:
        loss *= (2/(3*beta)-1/(3*beta**2))
    return loss

def get_round_filement_loss_par(T_cycle, Bm, Jc_theta, a):
    return get_round_filament_loss_per_cycle_par(Bm, Jc_theta, a)/Tcycle

def loss_factor_integrand(e_r):
    return e_r-np.arcsin(np.sqrt(1-e_r**2.))/np.sqrt(1-e_r**2.)

def beta_vs_e_m(e_m):
    return 1-e_m*np.arcsin(np.sqrt(1-e_m**2))/np.sqrt(1-e_m**2)

def get_beta_root_function(beta):

    def beta_root_function(e_m):
        return beta_vs_e_m(e_m)-beta

    return beta_root_function

def solve_e_m(beta):
    beta_root_function = get_beta_root_function(beta)

    e_m_sol = optimize.root(beta_root_function, 0)

    return e_m_sol['x'][0]

def round_filament_loss_factor(beta):
    #print "beta: ", beta
    e_m = solve_e_m(beta)
    #print "e_m:", e_m
    e_r_space = np.linspace(1,e_m,1e7)[1:]
    loss_integration = np.trapz(loss_factor_integrand(e_r_space), x=e_r_space)
    #print "loss_integration:",loss_integration
    loss_factor = 8/(3*beta**2) * loss_integration - 4/(3*beta)*(1-e_m**2)
    #print "loss_factor: ", loss_factor
    return loss_factor

def test_round_filament_loss_factor():
    T = 4.2
    B = 12
    beta = 1
    print "beta:", beta
    print "round filament loss factor:", round_filament_loss_factor(beta)

def get_round_filament_loss_per_cycle_perp(Bm, Jc, a):
    print "filament perpendicular loss:"
    beta = Bm*math.pi/(4*mu0*Jc*a)
    print "Bm:", Bm
    print "Bp:", Bm/beta
    print "beta:", beta
    loss = Bm**2./(2.*mu0)
    if beta <= 1:
        loss *= round_filament_loss_factor(beta)
    elif beta > 1:
        loss *= (4/(3*beta)-0.710/(beta**2.))
    return loss

def test_loss_per_cycle():
    matparam = NB3SNC_matparam
    T = 4.2
    B = 0.1
    Bc2 = get_Bc2(T, matparam)
    Jc = get_Jc(B, T, matparam)
    a = 55e-6/2.
    print "Bc2:", Bc2
    print "Jc:", Jc 
    print "filament half width:", a

    loss_per_cycle = get_round_filament_loss_per_cycle_perp(B, Jc, a)
    print "loss_per_cycle", loss_per_cycle

def get_round_filement_loss_perp(T_cycle, Bm, Jc_theta, a):
    return get_round_filament_loss_per_cycle_perp(Bm, Jc_theta, a)/Tcycle

""" Penetration loss filamentary composites [Wilson p.182] """

def composite_penetration_loss_per_cycle_perp(f, Bm, Jc, matparam, strandparam):
    RRR = strandparam['RRR']
    f_eff = get_f_eff_11T(False)
    rho0, rho1 = get_copper_rhos(RRR)
    rho_eff = get_rho_eff(f_eff, rho0, rho1, Bm)
    lf = strandparam['lf']
    a = strandparam['a']
    cu_to_noncu = strandparam['cu_to_noncu']
    omega = 2. * math.pi * f
    tau = get_tau_if(rho_eff, lf)
    Bm = 0.1
    cu_to_noncu = 1.106
    lambd = 1/(cu_to_noncu+1.)
    
    Bm_eff = Bm*omega*tau/np.sqrt(omega**2.*tau**2.+1)
    beta_eff = math.pi * Bm*omega*tau/(4*mu0*lambd*Jc*a*np.sqrt(omega**2.*tau**2+1))
    loss = Bm_eff**2./(2.*mu0)
    if beta_eff <= 1:
        loss *= round_filament_loss_factor(beta_eff)
    elif beta_eff > 1:
        loss *= (4/(3*beta_eff)-0.710/(beta_eff**2.))
    return loss

def composite_penetration_loss(f, Bm, Jc, matparam, strandparam):
    loss_per_cycle = composite_penetration_loss_per_cycle_perp(f, Bm, Jc, matparam, strandparam)
    q = loss_per_cycle * f
    return q
    
def test_composite_pen_loss_per_cycle():
    strandparam = HFMD07B_strandparam
    matparam = NB3SNC_matparam
    T = 1.9
    B = 0.1
    Bc2 = get_Bc2(T, matparam)
    Jc = get_Jc(B, T, matparam)
    f = 50.
    print "Bc2:", Bc2
    print "Jc:", Jc 
    print "Strand half width:", strandparam['a']

    loss_per_cycle = composite_penetration_loss_per_cycle_perp(f, B, Jc, matparam, strandparam)
    print "loss_per_cycle", loss_per_cycle

def test_composite_pen_loss():
    strandparam = HFMD07B_strandparam
    matparam = NB3SNC_matparam
    T = 1.9
    B = 0.01
    Bc2 = get_Bc2(T, matparam)
    Jc = get_Jc(B, T, matparam)
    f = 50.
    print "Bc2:", Bc2
    print "Jc:", Jc 
    print "Strand half width:", strandparam['a']

    q = composite_penetration_loss(f, B, Jc, matparam, strandparam)
    print "composite_penetration_loss", q


""" Inter-filament losses
parameter description (11T example)
lf is filament twist pitch (14 mm)
"""
def get_f_eff_11T(iph_res=False):
	if iph_res: return 0.4
	else: return 2.7
def get_f_eff(f_sc, iph_res=False):
	if iph_res: return (1.-f_sc)/(1.+f_sc)
	else: return (1.+f_sc)/(1.-f_sc)
def get_copper_rhos(RRR_in):
	RRR = [20 ,40 ,60 ,80 ,100,120,140,160,180,200]
	rho0 = [7.60E-10,3.70E-10,2.40E-10,1.80E-10,1.40E-10,1.20E-10,1.00E-10,8.70E-11,7.70E-11,6.80E-11]
	rho1 = [3.50E-11,3.80E-11,4.00E-11,4.10E-11,4.10E-11,4.10E-11,4.20E-11,4.20E-11,4.20E-11,4.20E-11]
	rho0_out = np.interp(RRR_in, RRR, rho0)
	rho1_out = np.interp(RRR_in, RRR, rho1)
	return rho0_out, rho1_out
def get_rho_eff(f_eff, rho0, rho1, B):
	return f_eff*(rho0 + rho1*B)
def get_tau_if(rho_eff, lf):
	return mu0/2./rho_eff * (lf/(2.*math.pi))**2.

def get_P_if(tau_if, dBdt):
	return 2.*tau_if/mu0 * dBdt**2.

def get_P_if_winkler(tau_if, B, f):
        omega = 2 * math.pi * f
	return 2.*math.pi*tau_if * omega * B**2 * f/(mu0*(1+(omega*tau_if)**2.))

########################## 

# cross-check winkler p.36
# reference 3.23 mW/cm^3 
# print get_P_if_winkler(26e-3, 10.3e-3, 50)*1e-2**3*1e3


""" Inter-Strand losses
parameter description (11T example)
b is cable thickness (1.25 mm)
c is cable width (14.7 mm)
N is the number of strands (40)
p is the cable twist pitch (100 mm)
Rc is the cable crossover resistance (typically between 10-300uohms)
Ra is the cable adjacent resistance (typically Rc/100) 
"""
def get_tau_is_a_par(b, p, c, Ra):
	return mu0*b*p/(16.*c*Ra)
def get_tau_is_a_perp(b, p, c, Ra):
	return mu0*c*p/(6.*b*Ra)
def get_tau_is_c_perp(N, b, p, c, Rc):
	return mu0*N*(N-1.)*c*p/(240.*b*Rc)
def get_P_is(tau, dBdt):
	return 2.*tau*dBdt**2./mu0

def get_tot_area(b, c, nof_cables):
	b = 1.25e-3
	c = 14.7e-3
	return b*c*nof_cables

########################## 

# 11T piece computation

def compute_11T_piece_loss(B0, Rc=3e-5, Ra=3e-7):
        angle_bface_field = 106./180.*math.pi
	RRR = 100
	lf = 14e-3
	#B0 = 10e-3
	f = 50.
	B_rms = B0/math.sqrt(2.)                      # sinusoidal signal
	dBdt_rms = math.sqrt(2.) * math.pi * f * B0
	dBdt_rms_par = abs(math.sqrt(2.) * math.pi * f * B0 * math.sin(angle_bface_field))
	dBdt_rms_perp = abs(math.sqrt(2.) * math.pi * f * B0 * math.cos(angle_bface_field))
        #print "dBdt", dBdt_rms * math.sqrt(2)
        #print "dBdt_rms", dBdt_rms
        #print "dBdt_rms_par", dBdt_rms_par
        #print "dBdt_rms_perp", dBdt_rms_perp

	b = 1.25e-3
	p = 100e-3
	c = 14.7e-3
	#Rc = 30e-6
	#Ra = Rc/100.
	N = 40
	# P_is
	tau_is_a_par = get_tau_is_a_par(b, p, c, Ra)
	tau_is_a_perp = get_tau_is_a_perp(b, p, c, Ra)
	tau_is_c_perp = get_tau_is_c_perp(N, b, p, c, Rc)
	P_is_a_par = get_P_is(tau_is_a_par, dBdt_rms_par) 
	P_is_a_perp = get_P_is(tau_is_a_perp, dBdt_rms_perp)
	P_is_c_perp = get_P_is(tau_is_c_perp, dBdt_rms_perp)

	# P_if
	f_eff = get_f_eff_11T(False)
	rho0, rho1 = get_copper_rhos(RRR)
	rho_eff = get_rho_eff(f_eff, rho0, rho1, B_rms)
	tau_if = get_tau_if(rho_eff, lf)
        print "tau_if:", tau_if
	P_if = get_P_if(tau_if, dBdt_rms)
	P_if_winkler = get_P_if_winkler(tau_if, B0, f)

	#nof_cables = 9+8+3+2+16+18
	nof_cables = 9+16
	tot_piece_area = get_tot_area(b, c, nof_cables)
        piece_len = 0.14
	#tot_quadrant_area = 5.93e1 * 1e-4 # cm^2 to m

	#print "total piece area 11T:", tot_piece_area
	#print "piece len 11T:", piece_len
	print "Piece volume 11T:", piece_len * tot_piece_area * 1e3**2
	#print "P_if [W] 11T piece: ", P_if * tot_piece_area * piece_len
	#print "P_is_a_par [W] 11T piece: ", P_is_a_par * tot_piece_area * piece_len
	#print "P_is_a_perp [W] 11T piece: ", P_is_a_perp * tot_piece_area * piece_len
	#print "P_is_c_perp [W] 11T piece: ", P_is_c_perp * tot_piece_area * piece_len

        # P_penetration
        strandparam = HFMD07B_strandparam
        matparam = NB3SNC_matparam
        T = 1.9
        B = 0.01
        Bc2 = get_Bc2(T, matparam)
        Jc = get_Jc(B, T, matparam)

#        print "Bc2:", Bc2
#        print "Jc:", Jc 
#        print "Strand half width:", strandparam['a']
#
        P_pen = composite_penetration_loss(f, B0, Jc, matparam, strandparam)
    
        Pif_piece = P_if * tot_piece_area * piece_len
        Pif_piece_winkler = P_if_winkler * tot_piece_area * piece_len
        
        Pis_piece = P_is_a_par * tot_piece_area * piece_len + P_is_a_perp * tot_piece_area * piece_len + P_is_c_perp * tot_piece_area * piece_len
        Ppen_piece = P_pen * tot_piece_area * piece_len
        return dBdt_rms, Pif_piece, Pif_piece_winkler, Pis_piece, Ppen_piece
        
def piece_11T_losses_test():
    Bmax_list = [10.3e-3,
                 9.8e-3,
                 9.1e-3,
                 8.1e-3,
                 6.9e-3,
                 5.5e-3,
                 4.8e-3,
                 4.0e-3,
                 17.9e-3,
                 17.0e-3,
                 15.8e-3,
                 14.1e-3,
                 12.0e-3,
                 9.6e-3,
                 8.3e-3,
                 7.0e-3]

    dBdt_rms_list = []
    Pif_list = []
    Pif_winkler_list = []
    Pis_list = []
    Ppen_list = []
    
    print "#dBdt_rms(T/s) Ppen(W) Pif(W) Pif_winkler(W) Pis(W) Ptot(W)"
    for Bmax in Bmax_list:
        dBdt_rms, Pif, Pif_winkler, Pis, Ppen = compute_11T_piece_loss(Bmax)
        Ptot = Pif + Pis + Ppen 
        print dBdt_rms, Ppen, Pif, Pif_winkler, Pis, Ptot
        dBdt_rms_list.append(dBdt_rms)
        Pif_list.append(Pif)
        Pif_winkler_list.append(Pif_winkler)
        Pis_list.append(Pis)
        Ppen_list.append(Ppen)

    plt.plot(Bmax_list, Pif_list)
    plt.show()

if __name__ == '__main__':
    piece_11T_losses_test()
    #test_round_filament_loss_factor()
    #test_summers()
    #test_loss_per_cycle()
    #test_composite_pen_loss_per_cycle()
    #test_composite_pen_loss()

