from neuron import h
from neuron.units import ms, mV
import numpy as np

class Cell:
    def __init__(self, gid, x, y, z, theta):
        self._gid = gid
        self._setup_morphology()
        self.all = self.soma.wholetree()
        self._setup_biophysics()
        self.x = self.y = self.z = 0  # <-- NEW
        # h.define_shape()
        # self._rotate_z(theta)  # <-- NEW
        # self._set_position(x, y, z)  # <-- NEW
        # everything below here in this method is NEW
        self._spike_detector = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        self.spike_times = h.Vector()
        self._spike_detector.record(self.spike_times)

        self._ncs = []

    def __repr__(self):
        return "{}[{}]".format(self.name, self._gid)


    def _set_position(self, x, y, z):
        for sec in self.all:
            for i in range(sec.n3d()):
                sec.pt3dchange(
                    i,
                    x - self.x + sec.x3d(i),
                    y - self.y + sec.y3d(i),
                    z - self.z + sec.z3d(i),
                    sec.diam3d(i),
                )
        self.x, self.y, self.z = x, y, z

    def _rotate_z(self, theta):
        """Rotate the cell about the Z axis."""
        for sec in self.all:
            for i in range(sec.n3d()):
                x = sec.x3d(i)
                y = sec.y3d(i)
                c = h.cos(theta)
                s = h.sin(theta)
                xprime = x * c - y * s
                yprime = x * s + y * c
                sec.pt3dchange(i, xprime, yprime, sec.z3d(i), sec.diam3d(i))



#ball and stick PN for cell

class PN(Cell):
    name = "PN"

    def _setup_morphology(self):
        self.soma = h.Section(name="soma", cell=self)
        #self.dend = h.Section(name="dend", cell=self)
        #self.dend.connect(self.soma)
        self.soma.L = 10.0/3.14
        self.soma.diam = 10.0


    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 100  # Axial resistance in Ohm * cm
            sec.cm = 1  # Membrane capacitance in micro Farads / cm^2
        self.soma.insert("hh_nrPN")
        for seg in self.soma:
            seg.hh_nrPN.gnabar = 7.15 # Sodium conductance in S/cm2
            seg.hh_nrPN.gkbar = 1.43  # Potassium conductance in S/cm2
            seg.hh_nrPN.gl = 0.021  # Leak conductance in S/cm2
            seg.hh_nrPN.el = -55  # Reversal potential in mV


#ball and stick for LN cell

class LN(Cell):
    name = "LN"

    def _setup_morphology(self):
        self.soma = h.Section(name="soma", cell=self)
        #self.dend = h.Section(name="dend", cell=self)
        #self.dend.connect(self.soma)
        self.soma.L = 10.0/3.14
        self.soma.diam = 10.0
        #self.dend.L = 100
        #self.dend.diam = 1

    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 100  # Axial resistance in Ohm * cm
            sec.cm = 1  # Membrane capacitance in micro Farads / cm^2
        self.soma.insert("hh_nrLN2")
        for seg in self.soma:
            seg.hh_nrLN2.gcabar = 0.286 # Sodium conductance in S/cm2
            seg.hh_nrLN2.gkbar = 1.43  # Potassium conductance in S/cm2
            seg.hh_nrLN2.gl = 0.021  # Leak conductance in S/cm2
            seg.hh_nrLN2.el = -50  # Reversal potential in mV
            seg.hh_nrLN2.gkcabar = 4.5*1e-2 # Ca dependemt K+ conductance in S/cm2

# make matrix for 10 LN and 30 PN name to map properly
def create_n_BallAndStick(cell_type,n, r):
    """n = number of cells; r = radius of circle"""
    cells = []
    for i in range(n):
        theta = i * 2 * h.PI / n
        cells.append(cell_type(i, h.cos(theta) * r, h.sin(theta) * r, 0, theta))
    return cells


class Ring_LN_PN :
    def __init__(self,nPN = 10, nLN = 10 ,sim_dur = 1000*ms):
        self.nPN = nPN
        self.nLN = nLN
        self.sim_dur = sim_dur
        self._create_cells()
        self._stimulate()
        self._insert_noise()
        self._connect_cells()


    def _create_cells(self) :
        self.PNcells = create_n_BallAndStick(PN,self.nPN, self.nPN*500000)

        self.LNcells = create_n_BallAndStick(LN,self.nLN, self.nLN*5000)

    def _stimulate(self, active_PN_frac = 1.0, active_LN_frac = 1.0
    ,PN_amp = 5,LN_amp = 5):

        self.iclamp_PN = []
        self.active_PN_list = []
        for i in range(self.nPN) :

            if (np.random.rand() < active_PN_frac) :
                self.active_PN_list = self.active_PN_list + ['active_PN' + str(i)]
                self.iclamp_PN = self.iclamp_PN + ['iclamp_PN' + str(i)]
                self.iclamp_PN[i] = h.IClamp(self.PNcells[i].soma(0.5))
                self.iclamp_PN[i].delay = 0.0
                self.iclamp_PN[i].dur = self.sim_dur
                self.iclamp_PN[i].amp = PN_amp

            else :
                self.iclamp_PN = self.iclamp_PN + ['iclamp_PN' + str(i)]
                self.iclamp_PN[i] = h.IClamp(self.PNcells[i].soma(0.5))
                self.iclamp_PN[i].delay = 0.0
                self.iclamp_PN[i].dur = self.sim_dur
                self.iclamp_PN[i].amp = 0


        self.iclamp_LN = []
        self.active_LN_list = []
        for i in range(self.nLN) :

            if (np.random.rand() < active_LN_frac) :
                self.active_LN_list = self.active_LN_list + ['active_LN' + str(i)]
                self.iclamp_LN = self.iclamp_LN + ['iclamp_LN' + str(i)]
                self.iclamp_LN[i] = h.IClamp(self.LNcells[i].soma(0.5))
                self.iclamp_LN[i].delay = 0.0
                self.iclamp_LN[i].dur = self.sim_dur
                self.iclamp_LN[i].amp = LN_amp

            else :
                self.iclamp_LN = self.iclamp_LN + ['iclamp_LN' + str(i)]
                self.iclamp_LN[i] = h.IClamp(self.LNcells[i].soma(0.5))
                self.iclamp_LN[i].delay = 0.0
                self.iclamp_LN[i].dur = self.sim_dur
                self.iclamp_LN[i].amp = 0

    def _stimulate_manually(self,active_PN_mat, active_LN_mat, PN_amp_list, LN_amp_list):

        self.iclamp_PN = []
        self.active_PN_list = []
        self.iclamp_LN = []
        self.active_LN_list = []

        for i in range(self.nPN) :

            if (active_PN_mat[i] == 1) :
                self.iclamp_PN = self.iclamp_PN + ['iclamp_PN' + str(i)]
                self.iclamp_PN[i] = h.IClamp(self.PNcells[i].soma(0.5))
                self.iclamp_PN[i].delay = 0.0
                self.iclamp_PN[i].dur = self.sim_dur
                self.iclamp_PN[i].amp = PN_amp_list[i]
                self.active_PN_list = self.active_PN_list + ['active_PN' + str(i)]
                
            else :
                self.iclamp_PN = self.iclamp_PN + ['iclamp_PN' + str(i)]
                self.iclamp_PN[i] = h.IClamp(self.PNcells[i].soma(0.5))
                self.iclamp_PN[i].delay = 0.0
                self.iclamp_PN[i].dur = self.sim_dur
                self.iclamp_PN[i].amp = 0

        for i in range(self.nLN) :

            if (active_LN_mat[i] == 1) :
                self.iclamp_LN = self.iclamp_LN + ['iclamp_LN' + str(i)]
                self.iclamp_LN[i] = h.IClamp(self.LNcells[i].soma(0.5))
                self.iclamp_LN[i].delay = 0.0
                self.iclamp_LN[i].dur = self.sim_dur
                self.iclamp_LN[i].amp = LN_amp_list[i]
                self.active_LN_list = self.active_LN_list + ['active_LN' + str(i)]

            else :
                self.iclamp_LN = self.iclamp_LN + ['iclamp_LN' + str(i)]
                self.iclamp_LN[i] = h.IClamp(self.LNcells[i].soma(0.5))
                self.iclamp_LN[i].delay = 0.0
                self.iclamp_LN[i].dur = self.sim_dur
                self.iclamp_LN[i].amp = 0

    def _stim_manual_arb(self,active_PN_mat, active_LN_mat,PN_arb_inp,LN_arb_inp):
        
        self.iclamp_PN = []
        self.iclamp_amp_PN = []
        self.inp_t1_PN = []
        self.active_PN_list = []
        self.iclamp_LN = []
        self.iclamp_amp_LN = []
        self.inp_t1_LN = []
        self.active_LN_list = []

        for i in range(self.nPN) :

            if (active_PN_mat[i] == 1) :
                self.iclamp_PN = self.iclamp_PN + ['iclamp_PN' + str(i)]
                self.iclamp_amp_PN = self.iclamp_amp_PN + ['iclamp_amp_PN' + str(i)]    
                self.iclamp_PN[i] = h.IClamp(self.PNcells[i].soma(0.5))
                self.iclamp_PN[i].delay = 0.0
                self.iclamp_PN[i].dur = 0.0
                self.iclamp_PN[i].amp = 0.0

                #insert the arbitrary input using vector play
                arb_extra = (self.sim_dur,0,PN_arb_inp[i])
                self.iclamp_amp_PN[i] = h.Vector(PN_arb_inp[i])
                self.inp_t1_PN = self.inp_t1_PN + ['inp_t1_PN' + str(i)]
                inp_t1 = np.arange(0,arb_extra[0],h.dt)  
                self.inp_t1_PN[i] = h.Vector(inp_t1)
                self.iclamp_PN[i].dur = arb_extra[0]
                self.iclamp_amp_PN[i].play(self.iclamp_PN[i]._ref_amp,self.inp_t1_PN[i],True)


            else :
                self.iclamp_PN = self.iclamp_PN + ['iclamp_PN' + str(i)]
                self.iclamp_PN[i] = h.IClamp(self.PNcells[i].soma(0.5))
                self.iclamp_PN[i].delay = 0.0
                self.iclamp_PN[i].dur = 0
                self.iclamp_PN[i].amp = 0

        for i in range(self.nLN) :

            if (active_LN_mat[i] == 1) :
                self.iclamp_LN = self.iclamp_LN + ['iclamp_LN' + str(i)]
                self.iclamp_amp_LN = self.iclamp_amp_LN + ['iclamp_amp_LN' + str(i)]    
                self.iclamp_LN[i] = h.IClamp(self.LNcells[i].soma(0.5))
                self.iclamp_LN[i] = h.IClamp(self.LNcells[i].soma(0.5))
                self.iclamp_LN[i].delay = 0.0
                self.iclamp_LN[i].dur = 0.0
                self.iclamp_LN[i].amp = 0.0

                #insert the arbitrary input using vector play
                arb_extra = (self.sim_dur,0,LN_arb_inp[i])
                self.iclamp_amp_LN[i] = h.Vector(LN_arb_inp[i])
                self.inp_t1_LN = self.inp_t1_LN + ['inp_t1_LN' + str(i)]
                inp_t1 = np.arange(0,arb_extra[0],h.dt)
                self.inp_t1_LN[i] = h.Vector(inp_t1)
                self.iclamp_LN[i].dur = arb_extra[0]
                self.iclamp_amp_LN[i].play(self.iclamp_LN[i]._ref_amp,self.inp_t1_LN[i],True)

            else :
                self.iclamp_LN = self.iclamp_LN + ['iclamp_LN' + str(i)]
                self.iclamp_LN[i] = h.IClamp(self.LNcells[i].soma(0.5))
                self.iclamp_LN[i].delay = 0.0
                self.iclamp_LN[i].dur = 0
                self.iclamp_LN[i].amp = 0


    def _insert_noise(self,nus_amp_PN = 0.5 ,nus_amp_LN = 0.5   ,nus_dur = 100*ms):


        self.noise_iclamp_PN = []
        self.noise_amp_PN = []
        self.noise_t_PN = []
        for i in range(self.nPN) :
            self.noise_iclamp_PN = self.noise_iclamp_PN + ['noise_iclamp_PN' + str(i)]
            self.noise_amp_PN = self.noise_amp_PN + ['noise_amp_PN' + str(i)]
            self.noise_t_PN = self.noise_t_PN + ['noise_t_PN' + str(i)]
            self.noise_iclamp_PN[i] = h.IClamp(self.PNcells[i].soma(0.5))
            self.noise_iclamp_PN[i].delay = 0.0
            self.noise_iclamp_PN[i].dur = 0.0
            self.noise_iclamp_PN[i].amp = 0.0


            noise_extra = (nus_dur,0,nus_amp_PN)
            noise_t1 = np.arange(0,noise_extra[0],h.dt)
            self.noise_amp_PN[i] =h.Vector(np.random.normal(noise_extra[1],noise_extra[2],noise_t1.shape))
            self.noise_t_PN[i] = h.Vector(noise_t1)
            self.noise_iclamp_PN[i].dur = noise_extra[0]
            self.noise_amp_PN[i].play(self.noise_iclamp_PN[i]._ref_amp,self.noise_t_PN[i],True)


        self.noise_iclamp_LN = []
        self.noise_amp_LN = []
        self.noise_t_LN = []
        for i in range(self.nLN) :
            self.noise_iclamp_LN = self.noise_iclamp_LN + ['noise_iclamp_LN' + str(i)]
            self.noise_amp_LN = self.noise_amp_LN + ['noise_amp_LN' + str(i)]
            self.noise_t_LN = self.noise_t_LN + ['noise_t_LN' + str(i)]
            self.noise_iclamp_LN[i] = h.IClamp(self.LNcells[i].soma(0.5))
            self.noise_iclamp_LN[i].delay = 0.0
            self.noise_iclamp_LN[i].dur = 0.0
            self.noise_iclamp_LN[i].amp = 0.0


            noise_extra = (nus_dur,0,nus_amp_LN)
            noise_t1 = np.arange(0,noise_extra[0],h.dt)
            self.noise_amp_LN[i] =h.Vector(np.random.normal(noise_extra[1],noise_extra[2],noise_t1.shape))
            self.noise_t_LN[i] = h.Vector(noise_t1)
            self.noise_iclamp_LN[i].dur = noise_extra[0]
            self.noise_amp_LN[i].play(self.noise_iclamp_LN[i]._ref_amp,self.noise_t_LN[i],True)



    def _connect_cells(self, pI2I = 1.0 , pI2E = 1.0 , pE2I = 1.0 , pE2E = 1.0,
    gI2I = 0.4 , gI2E_f = 0.4 , gI2E_s = 0.035 , gE2I = 0.30 , gE2E = 0.35):
        
        # Inhibitory to Inhibitory
        self.syn_LN2LN = []
        k=0
        for i in range(self.nLN) :

            for j in range(self.nLN) :
                if ((np.random.rand() <= pI2I) and (i !=j) ) :
                    self.syn_LN2LN = self.syn_LN2LN + ['syn_LN2LN' + str(k)]
                    self.syn_LN2LN[k] = h.nr_fGABA(self.LNcells[i].soma(0.5))
                    self.syn_LN2LN[k].gsyn = gI2I
                    h.setpointer(self.LNcells[j].soma(1)._ref_v, 'vpre', self.syn_LN2LN[k])
                    k = k + 1

        # Excitatory to Excitatory

        self.nc_E2E = []
        self.syn_PN2PN = []
        k=0
        for i in range(self.nPN) :
            for j in range(self.nPN) :
                if ((np.random.rand() <= pE2E) and (i !=j) ) :
                    self.syn_PN2PN = self.syn_PN2PN + ['syn_PN2PN' + str(k)]
                    self.nc_E2E = self.nc_E2E + ['nc_E2E' + str(k)]

                    self.syn_PN2PN[k] = h.ExpSyn(self.PNcells[i].soma(0.5))
                    self.syn_PN2PN[k].e = -0 * mV            # Nernst resting potential for a excitatory synapse
                    self.syn_PN2PN[k].tau = 6
                    self.nc_E2E[k] = h.NetCon(self.PNcells[j].soma(0.5)._ref_v, self.syn_PN2PN[k], sec = self.PNcells[j].soma)
                    self.nc_E2E[k].weight[0] =  gE2E
                    self.nc_E2E[k].delay = 0.0
                    self.nc_E2E[k].threshold = 30
                    k = k+1

        # Excitatory to Inhibitory
        self.nc_E2I = []
        self.syn_PN2LN = []


        k=0
        for i in range(self.nPN) :

            for j in range(self.nLN) :
                if (np.random.rand() < pE2I ) :
                    self.syn_PN2LN = self.syn_PN2LN + ['syn_PN2LN' + str(k)]
                    self.nc_E2I = self.nc_E2I + ['nc_E2I' + str(k)]

                    self.syn_PN2LN[k] = h.ExpSyn(self.LNcells[j].soma(0.5))
                    self.syn_PN2LN[k].e = -0 * mV            # Nernst resting potential for a Inhibitory synapse
                    self.syn_PN2LN[k].tau = 6
                    self.nc_E2I[k] = h.NetCon(self.PNcells[j].soma(0.5)._ref_v, self.syn_PN2LN[k], sec = self.PNcells[j].soma)
                    self.nc_E2I[k].weight[0] = gE2I
                    self.nc_E2I[k].delay = 0.0
                    self.nc_E2I[k].threshold = 35
                    k = k+1


        # Inhibitory  to Excitatory
        self.syn_LN2PN = []
        self.syn_LN2PN_slow = []

        k=0
        for i in range(self.nLN) :

            for j in range(self.nPN) :
                if (np.random.rand() < pI2E) :
                    self.syn_LN2PN = self.syn_LN2PN + ['syn_LN2PN' + str(k)]
                    #inihibitory synapses
                    self.syn_LN2PN[k] = h.nr_fGABA(self.PNcells[j].soma(0.5))
                    self.syn_LN2PN[k].gsyn = gI2E_f
                    h.setpointer(self.LNcells[i].soma(1)._ref_v, 'vpre', self.syn_LN2PN[k])

                    self.syn_LN2PN_slow = self.syn_LN2PN_slow + ['syn_LN2PN_slow' + str(k)]
                    #inihibitory synapses slow
                    self.syn_LN2PN_slow[k] = h.nr_sGABA(self.LNcells[i].soma(0.5))
                    self.syn_LN2PN_slow[k].gsyn = gI2E_s
                    #inhb_syn1s.gsyn =
                    h.setpointer(self.LNcells[i].soma(1)._ref_v, 'vpre', self.syn_LN2PN_slow[k])
                    k = k+1

    def _connect_cells_manually(self,connection_matrix, weight_matrix, connect_type, connect_index):
        self.syn_LN2LN = []
        self.syn_PN2PN = []
        self.syn_LN2PN = []
        self.syn_LN2PN_slow = []        
        self.syn_PN2LN = []
        self.nc_E2E = []
        self.nc_E2I = []
        self.nc_LN2LN = []
        self.nc_LN2PN = []
        self.nc_PN2LN = []
        self.nc_PN2PN = []

        
        for i in range(np.shape(connection_matrix)[0]):
            for j in range(np.shape(connection_matrix)[1]):

                if connection_matrix[i,j] == 1 :
                    k = 0
                    if connect_type[i,j] == "E2E":
                        i1 = connect_index[i,j][0]
                        i2 = connect_index[i,j][1]
                        self.syn_PN2PN = self.syn_PN2PN + ['syn_PN2PN' + str(k)]
                        self.nc_E2E = self.nc_E2E + ['nc_E2E' + str(k)]

                        self.syn_PN2PN[k] = h.ExpSyn(self.PNcells[i1].soma(0.5))
                        self.syn_PN2PN[k].e = -0 * mV
                        self.syn_PN2PN[k].tau = 6
                        self.nc_E2E[k] = h.NetCon(self.PNcells[i2].soma(0.5)._ref_v, self.syn_PN2PN[k], sec = self.PNcells[i2].soma)
                        self.nc_E2E[k].weight[0] =  weight_matrix[i,j]
                        self.nc_E2E[k].delay = 0.0
                        self.nc_E2E[k].threshold = 30
                        k = k+1
                    
                    if connect_type[i,j] == "I2I":
                        i1 = connect_index[i,j][0]
                        i2 = connect_index[i,j][1]
                        self.syn_LN2LN = self.syn_LN2LN + ['syn_LN2LN' + str(k)]
                        self.syn_LN2LN[k] = h.nr_fGABA(self.LNcells[i1].soma(0.5))
                        self.syn_LN2LN[k].gsyn = weight_matrix[i,j]
                        h.setpointer(self.LNcells[i2].soma(1)._ref_v, 'vpre', self.syn_LN2LN[k])
                        k = k + 1


                    if connect_type[i,j] == "E2I":
                        i1 = connect_index[i,j][0]
                        i2 = connect_index[i,j][1]
                        self.syn_PN2LN = self.syn_PN2LN + ['syn_PN2LN' + str(k)]
                        self.nc_E2I = self.nc_E2I + ['nc_E2I' + str(k)]

                        self.syn_PN2LN[k] = h.ExpSyn(self.LNcells[i2].soma(0.5))
                        self.syn_PN2LN[k].e = -0 * mV
                        self.syn_PN2LN[k].tau = 6
                        self.nc_E2I[k] = h.NetCon(self.PNcells[i1].soma(0.5)._ref_v, self.syn_PN2LN[k], sec = self.PNcells[i1].soma)
                        self.nc_E2I[k].weight[0] = weight_matrix[i,j]
                        self.nc_E2I[k].delay = 0.0
                        self.nc_E2I[k].threshold = 35
                        k = k+1

                    if connect_type[i,j] == "I2E":  

                        i1 = connect_index[i,j][0]
                        i2 = connect_index[i,j][1]
                        self.syn_LN2PN = self.syn_LN2PN + ['syn_LN2PN' + str(k)]
                        self.syn_LN2PN[k] = h.nr_fGABA(self.PNcells[i2].soma(0.5))
                        self.syn_LN2PN[k].gsyn = weight_matrix[i,j][0]
                        h.setpointer(self.LNcells[i1].soma(1)._ref_v, 'vpre', self.syn_LN2PN[k])
                        
                        self.syn_LN2PN_slow = self.syn_LN2PN_slow + ['syn_LN2PN_slow' + str(k)]
                        self.syn_LN2PN_slow[k] = h.nr_sGABA(self.PNcells[i2].soma(0.5))
                        self.syn_LN2PN_slow[k].gsyn = weight_matrix[i,j][1]
                        h.setpointer(self.LNcells[i1].soma(1)._ref_v, 'vpre', self.syn_LN2PN_slow[k])
                        k = k+1
