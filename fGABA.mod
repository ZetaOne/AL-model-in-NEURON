COMMENT
-----------------------------------------------------------------------------

  Kinetic model of fast GABA receptors for PN-LN connections based on PSST from neurorishika
  =================================


  Original model by : 
  Alain Destexhe, Salk Institute and Laval University, 1995
  ===================================

-----------------------------------------------------------------------------
ENDCOMMENT




NEURON {
	POINT_PROCESS nr_fGABA
	RANGE alpha, beta, e
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
	alpha = 10 (1/ms) <1e-9,1e9>
	beta = 0.16 (1/ms) <1e-9,1e9>
	e= -70	(mV)
	gsyn=0.8 (uS)
}

ASSIGNED {
	vpre (mV)
	v (mV)
	i (nA)
}

STATE {
	OP 
	
}

INITIAL {


    OP = 0.0

	                                              
}

BREAKPOINT {
	SOLVE state METHOD cnexp

	i = gsyn*OP*(v - e)
}

DERIVATIVE state {
	OP' = alpha * (1 - OP)*T(vpre) - beta* OP
}


NET_RECEIVE(weight (uS)) {
	OP = OP + weight
	
}


FUNCTION T(vpre) {  
 	T=  1/(1 + exp(-(vpre + 20)/1.5))
}

