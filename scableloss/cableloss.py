import math
import numpy as np
from scipy import optimize

mu0 = 4e-7 * math.pi
""" Hysterisis loss [Wilson]"""
def get_round_filament_loss_per_cycle_par(Bm, Jc_theta, a):
    beta = Bm/(2*mu0*Jc_theta*a)
    loss = Bm**2./(2.*mu_0)
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

    return e_m_sol['x']

def round_filament_loss_factor(beta):
    e_m=0
    e_r_space = np.linspace(1,e_m,1e7)[1:]
    loss_integration = np.trapz(loss_factor_integrand(e_r_space), x=e_r_space)
    return 8/(3*beta**2) * loss_integration - 4/(3*beta)*(1-e_m**2)

def test_round_filament_loss_factor():
    print solve_e_m(0.2)
    print round_filament_loss_factor(1)

def get_round_filament_loss_per_cycle_perp(Bm, Jc_theta, a):
    beta = Bm*math.pi/(4*mu0*Jc*a)
    loss = Bm**2./(2.*mu_0)
    if beta <= 1:
        loss *= round_filament_loss_factor(Bm,beta)
    elif beta > 1:
        loss *= (4/(3*beta)-0.710/(beta**2.))
    return loss

def get_round_filement_loss_perp(T_cycle, Bm, Jc_theta, a):
    return get_round_filament_loss_per_cycle_par(Bm, Jc_theta, a)/Tcycle


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
	return mu0/2./rho_eff * (lf/2.*math.pi)**2.
def get_P_if(tau_if, dBdt):
	return 2.*tau_if/mu0 * dBdt

########################## 

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
	return 2.*tau*dBdt/mu0

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
	P_if = get_P_if(tau_if, dBdt_rms)

	#nof_cables = 9+8+3+2+16+18
	nof_cables = 9+16
	tot_piece_area = get_tot_area(b, c, nof_cables)
        piece_len = 0.14
	#tot_quadrant_area = 5.93e1 * 1e-4 # cm^2 to m

	#print "total piece area 11T:", tot_piece_area
	#print "piece len 11T:", piece_len
	#print "Piece volume 11T:", piece_len * tot_piece_area
	#print "P_if [W] 11T piece: ", P_if * tot_piece_area * piece_len
	#print "P_is_a_par [W] 11T piece: ", P_is_a_par * tot_piece_area * piece_len
	#print "P_is_a_perp [W] 11T piece: ", P_is_a_perp * tot_piece_area * piece_len
	#print "P_is_c_perp [W] 11T piece: ", P_is_c_perp * tot_piece_area * piece_len

        Pif_piece = P_if * tot_piece_area * piece_len
        Pis_piece = P_is_a_par * tot_piece_area * piece_len + P_is_a_perp * tot_piece_area * piece_len + P_is_c_perp * tot_piece_area * piece_len
        return dBdt_rms, Pif_piece, Pis_piece
        
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

    print "#dBdt_rms(T/s) Pif(W) Pis(W) Ptot(W)"
    for Bmax in Bmax_list:
        dBdt_rms, Pif, Pis = compute_11T_piece_loss(Bmax)
        Ptot = Pif + Pis
        print dBdt_rms, Pif, Pis, Ptot


if __name__ == '__main__':
    test_round_filament_loss_factor()

