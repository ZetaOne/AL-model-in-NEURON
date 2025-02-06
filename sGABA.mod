COMMENT
-----------------------------------------------------------------------------

  Kinetic model of slow GABA receptors for PN-LN connections based on PSST from neurorishika
  =================================


  Original model by : 
  Alain Destexhe, Salk Institute and Laval University, 1995
  ===================================

-----------------------------------------------------------------------------
ENDCOMMENT




NEURON {
	POINT_PROCESS nr_sGABA
	RANGE r1,r2,r3,r4, beta, e
	NONSPECIFIC_CURRENT i
	POINTER vpre
	RANGE gsyn
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	r1	= 0.50	(/ms mM)	: forward binding rate to receptor
	r2	= 0.0013 (/ms)		: backward (unbinding) rate of receptor
	r3	= 0.01 (/ms)		: rate of G-protein production
	r4	= 0.033 (/ms)		: rate of G-protein decay
	KD	= 100			: dissociation constant of K+ channel
	ek= -95	(mV)
	gsyn=0.015 (uS)
}

ASSIGNED {
	vpre (mV)
	v (mV)
	i (nA)
}

STATE {

	R
	G 
	
}

INITIAL {


    
	R = 0.0
	G = 0.0

	                                              
}

BREAKPOINT {
	SOLVE state METHOD cnexp

	i = gsyn*(G*G*G*G)*(v - ek)/(G*G*G*G + KD)
}

DERIVATIVE state {
	R' = r1*(1 - R)*T(vpre) - r2*R
	G' = r3*R - r4*G
	
}


NET_RECEIVE(weight (uS)) {
	gsyn = gsyn + weight
	
}


FUNCTION T(vpre) {  
 	T =  1/(1 + exp(-(vpre + 20)/1.5))
}

