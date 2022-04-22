from qubitconfig import envelope_pulse

import inspect
import time
import copy
import json
import numpy as np
import re

class Qubit:
    """
    Todo: How standard should we make this? Maybe require 
    certain params (e.g. drive_freq/freq and read_freq)
    """
    def __init__(self, **paradict):
        self.phase={}
        for k, v in paradict.items():
            setattr(self, k, v)
            self.phase[k]=0 #todo: why is this here?

        self.qdrv=[] #todo: what are these?
        self.rdrv=[]
        self.read=[]
        self.mark=[]
    def add_phase(self, freqname, phase):
        self.phase[freqname]=(self.phase[freqname]+phase)%(2*np.pi)
    def add_time_phase(self, freqname, time):
        self.phase[freqname]=(self.phase[freqname]+2*np.pi*getattr(self, freqname)*time)%(2*np.pi)

class QChip:
    """
    Class containing configuration information for a chip described by a qchip json file
    (see examples/qubitcfg.json), as well as routines for calculating pulse envelopes. 

    Attributes
    ----------
        gates : dict
            dictionary of gate objects, with keys corresponding to those in json config
        qubits : dict
            dictionary of qubit objects
        gate_dict : dict
            dictionary of gate config information; key: gatename, value: list of pulse dicts
        paradict : dict
            config dictionary
    """
    def __init__(self, paradict):#, "gate1":Gate}):
        if isinstance(paradict, str):
            with open(paradict) as f:
                paradict = json.load(paradict)

        self.qubits={}
        for k, v in paradict['Qubits'].items():
            self.qubits.update({k:Qubit(**v)})

        gatesdict = {}
        for gatename in paradict['Gates']:
            gatesdict[gatename]=self._gatedictexpand(paradict['Gates'], gatename)

        self.gates = {}
        for k, pulselist in gatesdict.items():
            self.gates.update({k:Gate(pulselist, chip=self, name=k)})


    def _gatedictexpand(self, gatesdict, gatename):
        pulseout=[]
        for p in gatesdict[gatename]:
            if isinstance(p, str):
                if p in gatesdict:
                    pulseout.extend(self._gatedictexpand(gatesdict, p))
                else:
                    raise Exception('Error gate description: %s'%p)
            else:
                pulseout.append(p)
        return pulseout

    def save(self, wfilename):
        fwrite=open(wfilename, 'w')
        json.dump(self.paradict, fwrite, indent=4)
        print('c_qchip.updatecfg', wfilename)
        fwrite.close()

    def get_qubit_freq(self, freqname):
        if isinstance(freqname, str):
            m=re.match('(?P<qname>\S+)\.(?P<fname>\S+)', freqname)
            if m:
                q=m.group('qname')
                f=m.group('fname')
                if q in self.qubits:
                    if hasattr(self.qubits[q], f):
                        return getattr(self.qubits[q], f)

                    else:
                        raise AttributeError("Qubit %s doesn't have %s"%(str(q), str(f)))
                else:
                    raise AttributeError("Qubit %s not found"%str(q))

            else:
                raise Exception("%s does not match qubit.freqname format"%(str(freqname)))

        elif isinstance(freqname, float) or isinstance(freqname, int):
            return freqname #todo: do we actually want this to return int or should we do typechecking beforehand?

        else:
            raise TypeError("%s is not a string or number"%str(freqname))

    def getdest(self, destname):
        [q, d]=destname.split('.')
        return getattr(self.qubits[q], d)

    @property
    def qubit_dict(self):
        qdict = {}
        for name, obj in self.qubits.items():
            qdict[name] = vars(obj)
        return qdict

    @property
    def gate_dict(self):
        gdict = {}
        for name, obj in self.gates.items():
            gdict[name] = obj.paradict
        return gdict

    @property
    def paradict(self):
        return {'Qubits': self.qubit_dict, 'Gates': self.gate_dict}

class Gate:
    def __init__(self, pulse_list, chip, name):
        self.chip=chip
        self.name=name
        self.pulses=[]
        self._norm = 1.0

        for paradict in pulse_list:
            self.pulses.append(GatePulse(gate=self, paradict=paradict))

        self._tstart = None #todo: should this be 0?

    @property
    def norm(self, dt, FS):
        if len(self.pulses)==2: #todo: why hardcoded to 2
            p0 = self.pulses[0].pulseval(dt, self.pulses[0].fcarrier, mod=False).nvval(dt)[1]
            p1 = self.pulses[1].pulseval(dt, self.pulses[1].fcarrier, mod=False).nvval(dt)[1]
            self._norm = np.sum(p0*p1)*FS
        return self._norm

    @property
    def tstart(self):
        return self._tstart

    @tstart.setter
    def tstart(self, tstart):
        self._tstart=tstart
        for pulse in self.pulses:
            pulse.gate_tstart = tstart #todo: how to handle gate pulse tstarts

    @property
    def tlength(self):
        return max([p.t0+p.twidth for p in self.pulses]) - min([p.t0 for p in self.pulses])

    @property
    def tend(self):
        return max([pulse.tend() for pulse in self.pulses])#self.tstart+self.t0+self.tlength()

    @property
    def paradict(self):
        cfg = []
        for pulse in self.pulses:
            cfg.append(pulse.paradict)
        return cfg

    def pcalc(self, dt=0, padd=0, freq=None): #todo: what is this?
        return np.array([p.pcarrier+2*np.pi*(freq  if freq else p.fcarrier)*(dt+p.t0)+padd for p  in  self.pulses])

class GatePulse:
    def __init__(self, gate, paradict):
        '''
        t0: pulse start time relative to the gate start time
        tstart: pulse start time relative to circuit or circuits start time, this will be directly used for the trig_t
        twidth: pulse env function parameter for the pulse width
        patch: patch time added before the pulse, during this time, the pulse envelope value is 0
        tlength: pulse length including patch time, so that include the zeros at the begining, tlength=twidth+patch*dt
        '''
        self.gate=gate #todo: do we need this reference?
        self.chip=gate.chip
        self.dest=None #what is this/do we need it here?
        self.amp = 0
        self.twidth = 0
        self.t0 = 0
        for k in paradict:
            if k in ['amp', 'twidth', 't0', 'pcarrier', 'dest', 'fcarrier']: #do we want this restriction? what about other alphas?
                try:
                    v=eval(str(paradict[k]))
                except Exception as e:
                    v=paradict[k]
                setattr(self, k, v)
            elif k in ['env']: #should we require an envelope?
                setattr(self, k, Envelop(paradict[k]))

        self.gate_tstart = 0
        self.npatch=0 #todo: what are these; do we need them here?
        self.tpatch=0

    @property
    def tstart(self):
        return self.gate_tstart + self.t0

    def pulseval(self, dt, flo=0, mod=True):
        fcarrier = self.chip.get_qubit_freq(self.fcarrier)
        fif = 0 if self.fcarrier==0 else fcarrier - flo
        if self.dest is not None:
            tv=self.env.env_val(dt=dt, twidth=self.twidth, fif=0*fif, pini=self.pcarrier, amp=self.amp, mod=mod)
            tv.tstart(self.tstart)
        else:
            tv=None
        return tv

    @property
    def tend(self):
        return self.tstart + self.twidth

    @property
    def paradict(self):
        #cfg = vars(self)
        cfg = {}
        cfg['amp'] = self.amp
        cfg['twidth'] = self.twidth
        cfg['t0'] = self.t0
        cfg['pcarrier'] = self.pcarrier
        cfg['dest'] = self.dest
        cfg['fcarrier'] = self.fcarrier
        if hasattr(self, 'env'):
            cfg['env'] = self.env.env_desc
        return cfg

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        ampeq=False
        twidtheq=False
        desteq=False
        enveq=False
        if hasattr(self, 'env'):
            if hasattr(other, 'env'):
                if self.env==other.env:
                    enveq=True
        if hasattr(self, 'amp'):
            if hasattr(other, 'amp'):
                if self.amp==other.amp:
                    ampeq=True
        if hasattr(self, 'twidth'):
            if hasattr(other, 'twidth'):
                if self.twidth==other.twidth:
                    twidtheq=True
        if hasattr(self, 'dest'):
            if hasattr(other, 'dest'):
                if self.dest==other.dest:
                    desteq=True
#        return all([self.amp==other.amp, self.twidth==other.twidth, self.dest==other.dest, self.env==other.env])
        return all([ampeq, twidtheq, desteq, enveq])

    def __hash__(self):
        #return hash(str(self.paradict))
        return hash(str({k:self.paradict[k] for k in ("amp", "twidth", "t0", "dest") if k in self.paradict}))

    def timeoverlap(self, other):
        overlap=False
        if self.tstart<other.tstart:
            overlap=self.tend() < other.start
        else:
            overlap=other.tend() < self.tstart
        return overlap

class Envelop:
    def __init__(self, env_desc):
        if not isinstance(env_desc, list):
            env_desc=[env_desc]
        self.env_desc=copy.deepcopy(env_desc)

    def env_val(self, dt, twidth, fif=0, pini=0, amp=1.0, mod=True):
        #todo: what is vbase? and why would we have more than one env per pulse
        vbase=None 
        ti=None
        twidth=round(twidth, 10)
        for env in self.env_desc:
            if vbase is None:
                ti, vbase = getattr(envelope_pulse, env['env_func'])(dt=dt, twidth=twidth, **env['paradict'])
            else:
                ti1, vbasenew = (getattr(envelope_pulse, env['env_func'])(dt=dt, twidth=twidth, **env['paradict']))
                if any(ti!=ti1):
                    print('different time!!?')
                vbase = vbase*vbasenew
        if mod:
            vlo = np.cos(2*np.pi*fif*ti+pini)
        else:
            vlo = 1

        val=amp*vlo*vbase
        tv=c_tv(ti, val)
        return tv
    def __eq__(self, other):
        return sorted(self.env_desc)==sorted(other.env_desc)

class c_tv:
    #todo: what is this class?
    def __init__(self, t, val):
        self.tv={}
        if (len(t)==len(val)):
            for i in range(len(t)):
                self.append(t[i], val[i])
#        print t
#        print self.tv
        self.val=val
    def append(self, t, val):
        if t in self.tv:
            self.tv[t]+=val
        else:
            self.tv[t]=val
    def add(self, tval):
        for it in tval.tv:
            self.append(it, tval.tv[it])
    def tstart(self, tstart):
        newtv=dict((t+tstart, val) for (t, val) in self.tv.items())
        self.tv=newtv
    def tvval(self):
        self.t=[]
        self.val=[]
        for i in sorted(self.tv):
            self.t.append(i)
            self.val.append(self.tv[i])
        return(np.array(self.t), np.array(self.val))
    def nvval(self, dt):
        t, v=self.tvval()
        n=np.array([int(round(it/dt, 10)) for it in t])

#        return (np.array(np.round(t/dt)).astype(int), v)
        return (n, v)
