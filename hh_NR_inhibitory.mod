TITLE neurorishika_inhibitory.mod
 
COMMENT


        ***************************************************************************
        *  This file is based on the AL model of locust neurorishika_inhibitory neuron.
        *  The model is based on the Hodgkin-Huxley formalism and includes a calcium-dependent potassium current.
        *  The model was made in neuron by JKG-IISER Pune, India.
        *  ***************************************************************************
  ***************************************************************************
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

? interface

NEURON {
        SUFFIX hh_nrLN
        NONSPECIFIC_CURRENT ica,ik,il,ikl,ikca

        RANGE gcabar,gca,egca,gkbar,gk,egk,gl,el,gkcabar,gkca,ekca
	GLOBAL hinf, ninf, htau, ntau ,mtau ,minf ,mkinf ,mktau , CaCinf , CaCtau

}
 
PARAMETER {
        gcabar = 0.286  (mho/cm2)	<0,1e9>
	egca	= 140 (mV)	
        gkbar = 1.43 (mho/cm2)	<0,1e9>
	egk  = -95 (mV)	
        gl = 0.021 (mho/cm2)	<0,1e9>
        el = -50 (mV)
        gkl = 0.002 (mho/cm2)	<0,1e9>
        ekl = -95 (mV)		
        gkcabar = 0.0358 (mho/cm2)	<0,1e9>
        ekca = -90 (mV)	        	
}
 
STATE {
        m h n mk CaC
}
 
ASSIGNED {
        v (mV)
	celsius (degC)

	gca (mho/cm2)
        ica (mA/cm2)
	gk (mho/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        ikl (mA/cm2)
        gkca (mho/cm2)
        ikca (mA/cm2)
        minf hinf ninf mkinf CaCinf
	htau (ms) ntau (ms) mtau (ms) mktau (ms) CaCtau (ms)
}
 
LOCAL mexp, hexp, nexp , mkexp , CaCexp      
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
	
        gca = gcabar*m*m*h
	ica = gca*(v - egca)
        gk = gkbar*n*n*n*n
	ik = gk*(v - egk)      
        il = gl*(v - el)
        ikl = gkl*(v - ekl)
        gkca = gkcabar*mk
        ikca  = gkca*(v - ekca)*1.28
        
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
	mk = mkinf
        CaC = CaCinf        
}

? states
DERIVATIVE states {  
        rates(v)
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
        m' = (minf -m)/mtau            
        CaC' = (CaCinf -CaC)/CaCtau - 0.0002*ica
        mk' = ((0.01*CaC/(0.01*CaC + 0.02)) - mk )*1.28*(0.02 +0.01*CaC )
}
 
LOCAL q10 , q11


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
		      
        LOCAL  alpha, beta, sum ,vd
        TABLE minf, hinf, htau, ninf, ntau, mtau, minf, mkinf, mktau, CaCinf, CaCtau DEPEND celsius FROM -100 TO 100 WITH 200


UNITSOFF
        q10 = 3^((celsius - 36)/10)
        q11 = 2.3^((26 - 23)/10)
        

               :"m" calcium activation system
        vd = 0
        alpha = 0
        beta =  0
        sum = (1.0 + exp(-(v + 20)/6.5))
        minf = 1.0/sum
        mtau = 1.5

                :"h" calcium inactivation system
        vd = 0 
        alpha = 0.0
        beta = 0 
        sum = 0
	hinf = 1.0/( 1 + exp((v + 25)/12.0) ) 
        htau = 0.3*exp((v - 40)/13) + 0.002*exp((60 - v)/29.0)



                :"n" potassium activation system
        vd = v + 50
        alpha = 0.02*vtrap((15- vd),5) 
        beta = 0.5 * exp((10 -vd)/40.0)
	sum = alpha + beta
        ntau = 1/(q10*sum)
        ninf = alpha/sum



                :"mk" calcium dependent potassium current
        vd = 0
        alpha = 0.01
        beta = 0.02
        sum = 0.0
	mktau = 0.0
        mkinf = 0.0



                :"CaC" calcium concentration
        vd = 0
        alpha = 0
        beta = 0 
        sum = 0
	CaCtau = 150
        CaCinf = 0.00024              :0.00024        
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap =  y*(1 - x/y/2)    :vtrap =  y*(1 - x/y/2), technically it should be y*(1 -x/y), but anyways
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
 
UNITSON
