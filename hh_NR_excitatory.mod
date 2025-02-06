TITLE neurorishika_excitatory.mod   primary neuron sodium, potassium, and leak channels
 
COMMENT

 This file is based on the original hh.mod file (see original comment
 below). It was modified to match the model that was used in the
 simulations of Wang and Buzsaki (1996, J. Neurosci. 16).

 ***************************************************************************
 This is the original Hodgkin-Huxley treatment for the set of sodium, 
  potassium, and leakage channels found in the squid giant axon membrane.
  ("A quantitative description of membrane current and its application 
  conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544 (1952).)
 Membrane voltage is in absolute mV and has been reversed in polarity
  from the original HH convention and shifted to reflect a resting potential
  of -65 mV.
 Remember to set celsius=6.3 (or whatever) in your HOC file.
 See squid.hoc for an example of a simulation using this model.
 SW Jaslove  6 March, 1992
 ***************************************************************************

 changes:

  - m is substituted by it"s steady state value: m_inf - see 'BREAKPOINT'
  {as a result mtau is not needed, 'minf' is removed from
  GLOBAL declaration and 'm' is included in the RANGE var list
  otherwise it will be handled as a GLOBAL var and will not be
  evaluated separately for the 'sections'; for 'h' an 'n' this 
  is not a problem}

  - for h and n alpha and beta values are multiplied by 5 
  (see factor "Phi" in the W&B model)

  - USEION removed as we don't want to deal with ions and set eNa and
  eK directly. Rev potentials 'egna' and 'egk' are in the PARAMETERS
  list
    
  - temp: set to 6.3 Celsius, alpha and beta values are set/manipulated
  directly to simulate characteristic firing pattern

  I. Vida, Nov. 2000

  ***************************************************************************
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

? interface

NEURON {
        SUFFIX hh_nrPN
        NONSPECIFIC_CURRENT ina,ik,il,ikl,ia

        RANGE gnabar,gna,egna, gkbar,gk,egk, gl,el,ea,gabar
	GLOBAL hinf, ninf, htau, ntau ,mtau ,minf , mainf , hainf , matau , hatau

}
 
PARAMETER {
        gnabar = 7.15 (mho/cm2)	<0,1e9>
	egna	= 50 (mV)	
        gkbar = 1.43 (mho/cm2)	<0,1e9>
	egk	= -95 (mV)	
        gl = 0.021 (mho/cm2)	<0,1e9>
        el = -55 (mV)
        gkl = 0.00572 (mho/cm2)	<0,1e9>
        ekl = -95 (mV)		
        gabar = 1.43 (mho/cm2)	<0,1e9>
        ea = -95 (mV)	        	
}
 
STATE {
        m h n ma ha
}
 
ASSIGNED {
        v (mV)
	celsius (degC)

	gna (mho/cm2)
        ina (mA/cm2)
	gk (mho/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        ikl (mA/cm2)
        ga (mho/cm2)
        ia (mA/cm2)
        minf hinf ninf hainf mainf
	htau (ms) ntau (ms) mtau (ms) matau (ms) hatau (ms)
}
 
LOCAL mexp, hexp, nexp , haexp , maexp       
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
	
        gna = gnabar*(m^3)*h
	ina = gna*(v - egna)
        gk = gkbar*(n^4)
        ga = gabar*(ma^4)*ha
	ik = gk*(v - egk)      
        il = gl*(v - el)
        ikl = gkl*(v - ekl)
        ia  = ga*(v - ea)
        
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
	ha = hainf
	ma = mainf        
}

? states
DERIVATIVE states {  
        rates(v)
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
        m' = (minf -m)/mtau
        ha' = (hainf-ha)/hatau
        ma' = (mainf -ma)/matau
}
 
LOCAL q10 , q11


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
		      
        LOCAL  alpha, beta, sum ,vd
        TABLE minf, hinf, htau, ninf, ntau, mtau, minf, mainf, hainf DEPEND celsius FROM -100 TO 100 WITH 200


UNITSOFF
        q10 = 3^((22.0 - 36)/10.0)
        q11 = 3^((36.0 - 23.5)/10.0)
        

               :"m" sodium activation system
        vd = v + 50
        alpha = 0.32 * vtrap((13 - vd),4)
        beta =  0.28 * vtrap1((vd - 40),5)
        sum = alpha + beta
        minf = alpha/sum
        mtau = 1/(q10*sum)

                :"h" sodium inactivation system
        vd = v + 50
        alpha = 0.128*exp((17 -vd)/18.0)
        beta = 4.0/(exp((40 -vd)/5.0)+ 1.0 ) 
        sum = alpha + beta
	htau = 1/(q10*sum)
        hinf = alpha/sum

                :"n" potassium activation system
        vd = v + 50
        alpha = 0.02*vtrap((15- vd),5) 
        beta = 0.5*exp((10 -vd)/40.0)
	sum = alpha + beta
        ntau = 1.0/(q10*sum)
        ninf = alpha/sum


                :"ha" potassium transient current
        vd = 0
        alpha = 0
        beta = 0 
        sum = exp((v + 46.05)/5.0) + exp(-(v + 238.4)/37.45) 
        if (v < -63.0) {
                hatau = 1.0/(sum*q11)
        }else{
                hatau = 19.0/q11
        }
	
        hainf = 1.0/(1 + exp((v + 78.0)/6.0))

                :"ma" potassium transient current
        vd = 0
        alpha = 0
        beta = 0 
        sum = 0.37 + exp((v + 35.82)/19.69) + exp(-(v + 79.69)/12.7)
	matau = 1.0/(q11*sum)
        mainf = 1.0/(1 + exp(-(v + 60)/8.5))
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

FUNCTION vtrap1(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap1 = y*(1 - x/y/2)     :vtrap =  y*(1 - x/y/2), technically it should be y*(1 -x/y), but anyways
        }else{
                vtrap1 = x/(exp(x/y) - 1)
        }
}
 
 
UNITSON
