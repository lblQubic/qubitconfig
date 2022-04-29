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
    def __init__(self, freq=None, readfreq=None, **kwargs):
        self.freq = freq
        self.readfreq = readfreq
        for k, v in kwargs.items():
            setattr(self, k, v)

class QChip:
    """
    Class containing configuration information for a chip described by a qchip json file
    (see examples/qubitcfg.json), as well as routines for calculating pulse envelopes. 

    Attributes
    ----------
        gates : dict
            dictionary of Gate objects, with keys corresponding to those in json config
        qubits : dict
            dictionary of Qubit objects
        gate_dict : dict
            dictionary of gate config information; key: gatename, value: list of pulse dicts
        cfg_dict : dict
            config dictionary
    """
    def __init__(self, cfg_dict):#, "gate1":Gate}):
        if isinstance(cfg_dict, str):
            with open(cfg_dict) as f:
                cfg_dict = json.load(cfg_dict)

        self.qubits={}
        for k, v in cfg_dict['Qubits'].items():
            self.qubits.update({k:Qubit(**v)})

        self.gates = {}
        for k, pulselist in cfg_dict['Gates'].items():
            self.gates.update({k:Gate(pulselist, chip=self, name=k)})

    def save(self, wfilename):
        fwrite=open(wfilename, 'w')
        json.dump(self.cfg_dict, fwrite, indent=4)
        #print('c_qchip.updatecfg', wfilename)
        fwrite.close()

    def get_qubit_freq(self, freqname):
        """
        Get qubit frequency (could be read, drive, etc) from
        a qubit in the self.qubits dictionary.

        Parameters
        ----------
            freqname : str
                Frequency to get; format should be <QubitName>.<freqname>
        Returns
        -------
            int
                qubit frequency in Hz
        """
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
            gdict[name] = obj.cfg_dict
        return gdict

    @property
    def cfg_dict(self):
        return {'Qubits': self.qubit_dict, 'Gates': self.gate_dict}

class Gate:
    def __init__(self, contents, chip, name):
        self.chip=chip
        self.name=name
        self.contents=[]
        self._norm = 1.0

        for item_dict in contents:
            if 'gate' in item_dict.keys():
                self.contents.append(item_dict)
            else:
                self.contents.append(GatePulse(gate=self, cfg_dict=item_dict))

    @property
    def norm(self, dt, FS):
        if len(self.pulses)==2: #todo: why hardcoded to 2
            p0 = self.pulses[0].pulseval(dt, self.pulses[0].fcarrier, mod=False).nvval(dt)[1]
            p1 = self.pulses[1].pulseval(dt, self.pulses[1].fcarrier, mod=False).nvval(dt)[1]
            self._norm = np.sum(p0*p1)*FS
        return self._norm

    @property
    def tlength(self):
        return max([p.t0+p.twidth for p in self.get_pulses()]) - min([p.t0 for p in self.get_pulses()])

    @property
    def cfg_dict(self):
        cfg = []
        for item in self.contents:
            if isinstance(item, GatePulse):
                cfg.append(item.cfg_dict)
            else:
                cfg.append(item)
        return cfg

    def get_pulses(self, gate_t0=0):
        pulselist = []
        for item in self.contents:
            if isinstance(item, GatePulse):
                itemcpy = copy.deepcopy(item)
                itemcpy.t0 += gate_t0
                pulselist.append(itemcpy)
            else:
                pulselist.extend(self.chip.gates[item['gate']].get_pulses(gate_t0))
        return pulselist


    def pcalc(self, dt=0, padd=0, freq=None): #todo: what is this?
        return np.array([p.pcarrier+2*np.pi*(freq  if freq else p.fcarrier)*(dt+p.t0)+padd for p  in  self.pulses])

class GatePulse:
    def __init__(self, gate, cfg_dict):
        '''
        t0: pulse start time relative to the gate start time
        twidth: pulse env function parameter for the pulse width
        patch: patch time added before the pulse, during this time, the pulse envelope value is 0
        tlength: pulse length including patch time, so that include the zeros at the begining, tlength=twidth+patch*dt
        '''
        #self.gate=gate #todo: do we need this reference?
        #self.chip=gate.chip
        self.dest=None 
        self.amp = 0
        self.twidth = 0
        self.t0 = 0
        for k in cfg_dict:
            if k in ['amp', 'twidth', 't0', 'pcarrier', 'dest', 'fcarrier']: #do we want this restriction? what about other alphas?
                try:
                    v=eval(str(cfg_dict[k]))
                except Exception as e:
                    v=cfg_dict[k]
                setattr(self, k, v)
            elif k in ['env']: #should we require an envelope?
                setattr(self, k, Envelop(cfg_dict[k]))

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
    def cfg_dict(self):
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
        #return hash(str(self.cfg_dict))
        return hash(str({k:self.cfg_dict[k] for k in ("amp", "twidth", "t0", "dest") if k in self.cfg_dict}))

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
                ti, vbase = getattr(envelope_pulse, env['env_func'])(dt=dt, twidth=twidth, **env['cfg_dict'])
            else:
                ti1, vbasenew = (getattr(envelope_pulse, env['env_func'])(dt=dt, twidth=twidth, **env['cfg_dict']))
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
