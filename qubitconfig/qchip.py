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
    def __init__(self, cfg_dict):
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

    def add_gate(self, name, gate):
        if isinstance(gate, Gate):
            gate = copy.deepcopy(gate)
            gate.chip = self
            self.gates[name] = gate
        elif isinstance(gate, list):
            self.gates[name] = Gate(gate, self, name)

class Gate:
    def __init__(self, contents, chip, name):
        self.chip=chip
        self.name=name
        self.contents=[]

        for item_dict in contents:
            if 'gate' in item_dict.keys():
                self.contents.append(item_dict)
            else:
                self.contents.append(GatePulse(gate=self, chip=chip **item_dict))
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
                itemcpy.gate = self
                itemcpy.chip = self.chip
                pulselist.append(itemcpy)
            else:
                pulselist.extend(self.chip.gates[item['gate']].get_pulses(gate_t0))
        return pulselist

    def dereference_gates(self):
        self.contents = self.get_pulses()


    def pcalc(self, dt=0, padd=0, freq=None): #todo: what is this?
        return np.array([p.pcarrier+2*np.pi*(freq  if freq else p.fcarrier)*(dt+p.t0)+padd for p  in  self.pulses])

class GatePulse:
    """
    Class describing a calibrated quantum gate. 

    Attributes
    ----------

    Methods
    -------
    """
    def __init__(self, pcarrier, fcarrier, dest=None, amp=None, t0=None, twidth=None, env_desc=None, gate=None, chip=None):
        '''
        t0: pulse start time relative to the gate start time
        twidth: pulse env function parameter for the pulse width
        '''
        self.dest = dest 
        self.fcarrier = fcarrier
        self.pcarrier = pcarrier
        self.chip = chip
        self.gate = gate 
        if env_desc is not None: 
            self.env = Envelope(env_desc)
        if amp is not None: 
            self.amp = amp
        if t0 is not None: 
            self.t0 = t0
        if twidth is not None: 
            self.twidth = twidth

    def get_if_samples(self, dt, flo=0):
        pass
        #if isinstance(fcarrier, str):
        #    if self.chip is not None:
        #        fcarrier = self.chip.get_qubit_freq(self.fcarrier)
        #    else:
        #        raise Exception('Cannot resolve refernce to {}; add attribute for parent chip'.format(fcarrier))
        #fif = 0 if self.fcarrier==0 else fcarrier - flo
        #if self.dest is not None:
        #    tv=self.env.env_val(dt=dt, twidth=self.twidth, fif=0*fif, pini=self.pcarrier, amp=self.amp, mod=mod)
        #    tv.tstart(self.tstart)
        #else:
        #    tv=None
        #return tv
    def get_env_samples(self, dt):
        return self.env.get_samples(dt, self.twidth, self.amp)

    @property
    def tend(self):
        return self.tstart + self.twidth

    @property
    def cfg_dict(self):
        #cfg = vars(self)
        cfg = {}
        cfg['fcarrier'] = self.fcarrier
        cfg['pcarrier'] = self.pcarrier
        if hasattr(self, 'dest'):
            cfg['dest'] = self.dest
        if hasattr(self, 't0'):
            cfg['t0'] = self.t0
        if hasattr(self, 'amp'):
            cfg['amp'] = self.amp
        if hasattr(self, 'env'):
            cfg['env'] = self.env.env_desc
        return cfg

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        #todo: why don't we compare fcarrier and pcarrier here
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


class Envelope:
    def __init__(self, env_desc):
        if not isinstance(env_desc, list):
            env_desc=[env_desc]
        self.env_desc=copy.deepcopy(env_desc)

    def get_samples(self, dt, twidth, amp=1.0):
        samples = None
        tlist = None
        twidth = round(twidth, 10) #todo: why do we round this?

        for env in self.env_desc:
            ti, vali = getattr(envelope_pulse, env['env_func'])(dt=dt, twidth=twidth, **env['cfg_dict'])
            if samples:
                samples *= vali
                assert np.all(ti==tlist)
            else:
                samples = vali
                tlist = ti

        samples *= amp

        return tlist, samples

    def __eq__(self, other):
        return sorted(self.env_desc)==sorted(other.env_desc)

