"""
Zero-dimensional reactors.
"""

import _cantera
from Numeric import array, zeros
import types


class ReactorBase:
    """Base class for reactors and reservoirs.

    Classes Reactor and Reservoir derive from a common base class
    ReactorBase. They have the same set of methods, which are all
    inherited from ReactorBase.

    (This is not quite true in the corresponding classes in the
    Cantera C++ kernel. There class Reactor defines some methods that
    class Reservoir doesn't. These are used internally by the
    ReactorNet instance that integrates the system of ODEs describing
    the network to evaluate the portion of the ODE system associated
    with that reactor.)
    """

    def __init__(self, name = '', contents = None,
                 volume = 1.0, energy = 'on',
                 type = -1, verbose = 0):
        """
        See class 'Reactor' for a description of the constructor parameters.
        The 'type' parameter specifies whether a Reactor (type = 1) or
        Reservoir (type = 2) will be created.
        """
        self.__reactor_id = _cantera.reactor_new(type)
        self._type = type
        self._inlets = []
        self._outlets = []
        self._walls = []
        self._reservoirs = []
        self._name = name
        self._verbose = verbose
        self.insert(contents)
        self._setInitialVolume(volume)
        self._setEnergy(energy)
        if self._verbose:
            print 'Created '+self._name
            print '    Volume = ',volume,' m^3'
            if energy <> 'on':
                print '    Temperature will be held constant'
            print '    Initial State:'
            print contents


    def __del__(self):
        """Delete the reactor instance."""
        if self._verbose:
            print 'Deleting '+self._name
        _cantera.reactor_del(self.__reactor_id)

    def __str__(self):
        s = self._name
        s += ':\n     Volume = '+`self.volume()`
        if self._contents:
            s += "\n"+`self._contents`
        return s
    
    def __repr__(self):
        s = self._name
        s += ':\n     Volume = '+`self.volume()`        
        if self._contents:
            s += ": \n"+`self._contents`
        return s
        
    def name(self):
        """The name of the reactor."""
        return self._name
    
    def reactor_id(self):
        """The integer index used to access the kernel reactor
        object. For internal use."""
        return self.__reactor_id
    
    def insert(self, contents):
        """
        Insert 'contents' into the reactor. Sets the objects used to compute
        thermodynamic properties and kinetic rates.
        """
        # store a reference to contents so that it will live as long
        # as this object
        self._contents = contents
        if contents:
            _cantera.reactor_setThermoMgr(self.__reactor_id, contents._phase_id)
            _cantera.reactor_setKineticsMgr(self.__reactor_id, contents.ckin)

        
##     def setInitialTime(self, T0):
##         """Deprecated.
##         Set the initial time. Restarts integration from this time
##         using the current state as the initial condition. Default: 0.0 s"""
##         _cantera.reactor_setInitialTime(self.__reactor_id, T0)

    def _setInitialVolume(self, V0):
        """Set the initial reactor volume. """
        _cantera.reactor_setInitialVolume(self.__reactor_id, V0)

    def _setEnergy(self, eflag):
        """Turn the energy equation on or off. If the argument is the
        string 'off' or the number 0, the energy equation is disabled,
        and the reactor temperature is held constant at its initial
        value."""
        ie = 1
        if eflag == 'off' or eflag == 0:
            ie = 0
        if self._verbose:
            if ie:
                print 'enabling energy equation for reactor',self._name
            else:
                print 'disabling energy equation for reactor',self._name                
        _cantera.reactor_setEnergy(self.__reactor_id, ie)

    def temperature(self):
        """The temperature in the reactor [K]."""
        return _cantera.reactor_temperature(self.__reactor_id)

    def density(self):
        """The density of the fluid in the reactor [kg/m^3]."""
        return _cantera.reactor_density(self.__reactor_id)

    def volume(self):
        """The total reactor volume [m^3]. The volume may change with time
        if non-rigid walls are installed on the reactor."""
        return _cantera.reactor_volume(self.__reactor_id)            

    def time(self):
        """Deprecated. The current time [s]."""
        return _cantera.reactor_time(self.__reactor_id)

    def mass(self):
        """The total mass of fluid in the reactor [kg]."""
        return _cantera.reactor_mass(self.__reactor_id)

    def enthalpy_mass(self):
        """The specific enthalpy of the fluid in the reactor [J/kg]."""
        return _cantera.reactor_enthalpy_mass(self.__reactor_id)

    def intEnergy_mass(self):
        """The specific internal energy of the fluid in the reactor [J/kg]."""
        return _cantera.reactor_intEnergy_mass(self.__reactor_id)

    def pressure(self):
        """The pressure in the reactor [Pa]."""
        return _cantera.reactor_pressure(self.__reactor_id)            

    def advance(self, time):
        """Deprecated.
        Advance the state of the reactor in time from the current
        time to time 'time'. Note: this method is deprecated. See
        class ReactorNet."""
        return _cantera.reactor_advance(self.__reactor_id, time)

    def step(self, time):
        """Deprecated.
        Take one internal time step from the current time toward
        time 'time'. Note: this method is deprecated. See class
        ReactorNet."""
        return _cantera.reactor_step(self.__reactor_id, time)    
    
    def massFraction(self, s):
        """The mass fraction of species s, specified either by name or
        index number.
        >>> y1 = r.massFraction(7)
        ___0.02
        >>> y2 = r.massFraction('CH3O')
        ___0.02
        """
        if type(s) == types.StringType:
            kk = self._contents.speciesIndex(s)
        else:
            kk = s
        return _cantera.reactor_massFraction(self.__reactor_id, kk)

    def massFractions(self):
        """Return an array of the species mass fractions."""
        nsp = self._contents.nSpecies()
        y = zeros(nsp,'d')
        for k in range(nsp):
            y[k] = self.massFraction(k)
        return y

    def moleFractions(self):
        """Return an array of the species mole fractions."""
        y = self.massFractions()
        self._contents.setMassFractions(y)
        return self._contents.moleFractions()

    def moleFraction(self, s):
        """The mole fraction of species s, specified either by name or
        index number.
        >>> x1 = r.moleFraction(9)
        ___0.00012
        >>> x2 = r.moleFraction('CH3')
        ___0.00012
        """        
        if type(s) == types.StringType:
            kk = self._contents.speciesIndex(s)
        else:
            kk = s
        x = self.moleFractions()
        return x[kk]
        
    def inlets(self):
        """Return the list of flow devices installed on inlets to this reactor.
        This method can be used to access information about the flows entering
        the reactor:
        >>> for n in r.inlets():
        ...    print n.name(), n.massFlowRate()
        See: MassFlowController, Valve, PressureController.
        """        
        return self._inlets

    def outlets(self):
        """Return the list of flow devices installed on outlets
        on this reactor.
        >>> for o in r.outlets():
        ...    print o.name(), o.massFlowRate()
        See: MassFlowController, Valve, PressureController.        
        """
        return self._outlets

    def walls(self):
        """Return the list of walls installed on this reactor.
        >>> for w in r.walls():
        ...    print w.name()
        See: Wall.
        """
        return self._walls
    
    def _addInlet(self, inlet, other):
        """For internal use. Store a reference to 'inlet'
        so that it will not be deleted before this object."""        
        self._inlets.append(inlet)
        if self._type == 1 and other._type == 2:
            self._reservoirs.append(other)        

    def _addOutlet(self, outlet, other):
        """For internal use. Store a reference to 'outlet'
        so that it will not be deleted before this object."""        
        self._outlets.append(outlet)
        if self._type == 1 and other._type == 2:
            self._reservoirs.append(other)        

    def _addWall(self, wall, other):
        """For internal use. Store a reference to 'wall'
        so that it will not be deleted before this object."""
        self._walls.append(wall)
        if self._type == 1 and other._type == 2:
            self._reservoirs.append(other)

    def syncContents(self):
        """Set the state of the object representing the reactor contents
        to the current reactor state.
        >>> r = Reactor(gas)
        >>> (statements that change the state of object 'gas')
        >>> r.syncContents()
        After this statement, the state of object 'gas' is synchronized
        with the reactor state.
        See 'contents'.
        """
        self._contents.setState_TRY(self.temperature(),
                                   self.density(),
                                   self.massFractions())
        
    def contents(self):
        """Return an object representing the reactor contents, after first
        synchronizing its state with the current reactor state. This method
        is useful when some property of the fluid in the reactor is
        needed that is not provided by a method of class Reactor.
        >>> r = Reactor(gas)
        >>> (statements that change the state of object 'gas')
        >>> c = r.contents()
        >>> print c.gibbs_mole(), c.chemPotentials()
        
        Note that after calling method 'contents', object 'c'
        references the same underlying kernel object as object 'gas'
        does. Therefore, all properties of 'c' and 'gas' are
        identical. (Remember that Python objects are really C
        pointers; at the C level, both point to the same data
        structure.)
        It is also allowed to write
        >>> gas = r.contents()
        """
        self.syncContents()
        return self._contents
    

_reactorcount = 0
_reservoircount = 0

class Reactor(ReactorBase):
    """
    Zero-dimensional reactors. Instances of class Reactor represent
    zero-dimensional reactors. By default, they are closed (no inlets
    or outlets), have fixed volume, and have adiabatic, chemically-intert
    walls. These properties may all be changed by adding appropriate
    components.
    See classes 'Wall', 'MassFlowController', and 'Valve'.
    """
    def __init__(self, contents = None, name = '',
                 volume = 1.0, energy = 'on',
                 verbose = 0):
        """
        contents - Reactor contents. If not specified, the reactor is
        initially empty. In this case, call method 'insert' to specify
        the contents.

        name - Used only to identify this reactor in output. If not
        specified, defaults to 'Reactor_n', where n is an integer
        assigned in the order Reactor objects are created.

        volume - Initial reactor volume. Defaults to 1 m^3.

        energy - Set to 'on' or 'off'. If set to 'off', the energy
        equation is not solved, and the temperature is held at its
        initial value. The default in 'on'.

        verbose - if set to a non-zero value, additional diagnostic
        information will be printed.

        Some examples showing how to create Reactor objects are shown below.
        >>> gas = GRI30()
        >>> r1 = Reactor(gas)
        This is equivalent to:
        >>> r1 = Reactor()
        >>> r1.insert(gas)
        Arguments may be specified using keywords in any order:
        >>> r2 = Reactor(contents = gas, energy = 'off',
        ...              name = 'isothermal_reactor')
        >>> r3 = Reactor(contents = gas, name = 'adiabatic_reactor')
        Here's an array of reactors:
        >>> reactor_array = [Reactor(), Reactor(gas), Reactor(Air())]
        """
        global _reactorcount
        if name == '':
            name = 'Reactor_'+`_reactorcount`
        _reactorcount += 1
        ReactorBase.__init__(self, contents = contents, name = name,
                             volume = volume, energy = energy,
                             verbose = verbose, type = 1)
            

class Reservoir(ReactorBase):
    """
    A reservoir is a reactor with a constant state. The temperature,
    pressure, and chemical composition in a reservoir never change from
    their initial values.
    """
    def __init__(self, contents = None, name = '', verbose = 0):
        """
        contents - Reservoir contents. If not specified, the reservoir is
        initially empty. In this case, call method insert to specify
        the contents.

        name - Used only to identify this reservoir in output. If not
        specified, defaults to 'Reservoir_n', where n is an integer
        assigned in the order Reservoir objects are created.

        verbose - if set to a non-zero value, additional diagnostic
        information will be printed.

        Some examples showing how to create Reservoir objects are shown below.
        >>> gas = GRI30()
        >>> res1 = Reservoir(gas)
        This is equivalent to:
        >>> res1 = Reactor()
        >>> res1.insert(gas)
        Arguments may be specified using keywords in any order:
        >>> res2 = Reservoir(contents = Air(), 
        ...                  name = 'environment')
        >>> res3 = Reservoir(contents = gas, name = 'upstream_state')
        """        
        global _reservoircount
        if name == '':
            name = 'Reservoir_'+`_reservoircount`
        _reservoircount += 1
        ReactorBase.__init__(self, contents = contents,
                             name = name, verbose = verbose, type = 2)
            
    def advance(self, time):
        """Deprecated. Do nothing."""
        pass




#------------------ FlowDevice ---------------------------------

class FlowDevice:
    """
    Base class for devices that regulate the flow rate in a fluid line.
    """
    def __init__(self, type, name, verbose):
        """
        Create a new instance of type 'type'
        """
        self._name = name
        self._verbose = verbose
        self.__fdev_id = _cantera.flowdev_new(type)

    def __del__(self):
        """
        Delete the instance.
        """
        if self._verbose:
            print 'deleting '+self._name
        _cantera.flowdev_del(self.__fdev_id)
        
    def name(self):
        """The name specified when initially constructed."""
        return self._name
    
    def ready(self):
        """
        Deprecated. Returns true if the device is ready to use.
        """
        return _cantera.flowdev_ready(self.__fdev_id)

    def massFlowRate(self, time = -999.0):
        """Mass flow rate (kg/s). """
        return _cantera.flowdev_massFlowRate(self.__fdev_id, time)

    def install(self, upstream, downstream):
        """
        Install the device between the upstream and downstream
        reactors or reservoirs.  
        >>> f.install(upstream = reactor1, downstream = reservoir2)
        """
        if self._verbose:
            print
            print self._name+': installing between '+upstream.name()+' and '+downstream.name()
        upstream._addOutlet(self, downstream)
        downstream._addInlet(self, upstream)
        _cantera.flowdev_install(self.__fdev_id, upstream.reactor_id(),
                                  downstream.reactor_id())
    def _setParameters(self, c):
        params = array(c,'d')
        n = len(params)
        return _cantera.flowdev_setParameters(self.__fdev_id, n, params)    

    def setFunction(self, f):
        _cantera.flowdev_setFunction(self.__fdev_id, f.func_id())

    def flowdev_id(self):
        return self.__fdev_id
        
_mfccount = 0

class MassFlowController(FlowDevice):
    
    """Mass flow controllers. A mass flow controller maintains a
    specified mass flow rate independent of upstream and downstream
    conditions. The equation used to compute the mass flow rate is
    \f[
    \dot m = \max(\dot m_0, 0.0),
    \f] where \f$ \dot m_0 \f$ is either
    a constant value or a function of time. Note that if \f$\dot m_0 <
    0\f$, the mass flow rate will be set to zero, since reversal of
    the flow direction is not allowed.
    
    Unlike a real mass flow controller, a MassFlowController object
    will maintain the flow even if the downstream pressure is greater
    than the upstream pressure.  This allows simple implementation of
    loops, in which exhaust gas from a reactor is fed back into it
    through an inlet. But note that this capability should be used
    with caution, since no account is taken of the work required to do
    this.

    A mass flow controller is assumed to be adiabatic, non-reactive,
    and have negligible volume, so that it is internally always in
    steady-state even if the upstream and downstream reactors are
    not. The fluid enthalpy, chemical composition, and mass flow rate
    are constant across a mass flow controller, and the pressure
    difference equals the difference in pressure between the upstream
    and downstream reactors.

    Examples:
    
    >>> mfc1 = MassFlowController(upstream = res1, downstream = reactr,
    ...                           name = 'fuel_mfc', mdot = 0.1)
    >>> air_mdot = Gaussian(A = 0.1, t0 = 2.0, FWHM = 0.1)
    >>> mfc2 = MassFlowController(upstream = res2, downstream = reactr,
    ...                           name = 'air_mfc', mdot = air_mdot)
    
    """
    def __init__(self, upstream=None,
                 downstream=None,
                 name='',
                 verbose=0, mdot = 0.0):
        """
        upstream - upstream reactor or reservoir.

        downstream - downstream reactor or reservoir.

        name - name used to identify the mass flow controller in output.
        If no name is specified, it defaults to 'MFC_n', where n is an
        integer assigned in the order the MassFlowController object
        was created.

        mdot - Mass flow rate [kg/s]. This mass flow rate, which may
        be a constant of a function of time, will be maintained,
        independent of unstream and downstream conditions, unless
        reset by calling method 'set'.

        verbose - if set to a positive integer, additional diagnostic
        information will be printed.
        
        """
        global _mfccount
        if name == '':
            name = 'MFC_'+`_mfccount`
        _mfccount += 1
        FlowDevice.__init__(self,1,name,verbose)
        if upstream and downstream:
            self.install(upstream, downstream)
        if mdot:
            self.set(mdot = mdot)

    def _setMassFlowRate(self, mdot):
        """Set or reset the mass flow rate to 'mdot' [kg/s].
        """
        if self._verbose:
            print self._name+': setting mdot to '+`mdot`+' kg/s'
        if type(mdot) == types.InstanceType:
            self.setFunction(mdot)
        else:
            _cantera.flowdev_setMassFlowRate(self.flowdev_id(), mdot)


    def set(self, mdot = 0.0):
        """Set the mass flow rate [kg/s]. May be called at any time to
        change the mass flow rate to a new value, or to a new function
        of time.

        >>> mfc.set(mdot = 0.2)
        """
        self._setMassFlowRate(mdot)


_valvecount = 0

class Valve(FlowDevice):
    """Valves. In Cantera, a Valve object is a flow devices with mass
    flow rate that is a function of the pressure drop across it. The default behavior
    is linear:
    \f[ \dot m = K_v (P_1 - P_2) \f]
    if \f$ P_1 > P_2. \f$
    Otherwise,
    \f$ \dot m = 0 \f$.
    However, an arbitrary function \f$ F\f$ can also be specified, such that
    \f[
    \dot m = F(P_1 - P_2).
    \f]
    if \f$ P_1 > P_2, \f$
    or  \f$ \dot m = 0 \f$ otherwise.    
    It is never possible for the flow to reverse
    and go from the downstream to the upstream reactor/reservoir through
    a line containing a Valve object. 

    'Valve' objects are often used between an upstream reactor and a
    downstream reactor or reservoir to maintain them both at nearly the
    same pressure. By setting the constant \f$ K_v \f$ to a
    sufficiently large value, very small pressure differences will
    result in flow between the reactors that counteracts the pressure
    difference.

    A Valve is assumed to be adiabatic, non-reactive, and have
    negligible internal volume, so that it is internally always in
    steady-state even if the upstream and downstream reactors are
    not. The fluid enthalpy, chemical composition, and mass flow rate
    are constant across a Valve, and the pressure difference equals
    the difference in pressure between the upstream and downstream
    reactors.
    
    """
    def __init__(self, upstream=None, downstream=None,
                 name='', Kv = 0.0, mdot0 = 0.0, verbose=0):
        """
        upstream - upstream reactor or reservoir.
        
        downstream - downstream reactor or reservoir.
        
        name - name used to identify the valve in output.
        If no name is specified, it defaults to 'Valve_n', where n is an
        integer assigned in the order the Valve object
        was created.

        Kv - the constant in the mass flow rate equation.

        verbose - if set to a positive integer, additional diagnostic
        information will be printed.
        
        """        
        global _valvecount
        if name == '':
            name = 'Valve_'+`_valvecount`
        _valvecount += 1
        FlowDevice.__init__(self,3,name,verbose)
        if upstream and downstream:
            self.install(upstream, downstream)
        self.setValveCoeff(Kv, mdot0)


    def setValveCoeff(self, Kv = -1.0):
        """Set or reset the valve coefficient \f$ K_v \f$."""
        vv = zeros(1,'d')
        vv[0] = Kv
        if self._verbose:
            print
            print self._name+': setting valve coefficient to '+`Kv`+' kg/Pa-s'
        self._setParameters(vv)

    def _setValveCharacteristic(self, f):
        """Set or reset the valve characteristics.
        """
        if type(f) == types.InstanceType:
            self.setFunction(f)
        else:
            raise CanteraError("Wrong type for valve characteristic function.")

    def set(self, Kv = -1.0, F = None):
        """Set or reset valve properties. All keywords are optional.

        Kv - constant in linear mass flow rate equation.

        F - function of \f$\Delta P\f$.
        """
        if F:
            self.setFunction(F)
        if Kv > 0.0:
            self.setValveCoeff(Kv)
            


_pccount = 0

class PressureController(FlowDevice):

    """ A PressureController is designed to be used in conjunction
    with another 'master' flow controller, typically a
    MassFlowController. The master flow controller is installed on the
    inlet of the reactor, and the corresponding PressureController is
    installed on on outlet of the reactor. The PressureController mass
    flow rate is equal to the master mass flow rate, plus a
    small correction dependent on the pressure difference:
    \f[
    \dot m = \dot m_{\rm master} + K_v(P_1 - P_2).
    \f]
    """

    def __init__(self, upstream=None, downstream=None,
                 name='', master = None, Kv = 0.0, verbose=0):
        """
        upstream - upstream reactor or reservoir.
        
        downstream - downstream reactor or reservoir.
        
        name - name used to identify the pressure controller in
        output.  If no name is specified, it defaults to
        'PressureController_n', where n is an integer assigned in the
        order the PressureController object was created.

        Kv - the constant in the mass flow rate equation.

        verbose - if set to a positive integer, additional diagnostic
        information will be printed.
        
        """        
        global _pccount
        if name == '':
            name = 'PressureController_'+`_pccount`
        _pccount += 1
        FlowDevice.__init__(self,2,name,verbose)
        if upstream and downstream:
            self.install(upstream, downstream)
        self.setPressureCoeff(Kv)
        self.setMaster(master)


    def setPressureCoeff(self, Kv):
        """Set or reset the pressure coefficient \f$ K_v \f$."""
        vv = zeros(1,'d')
        vv[0] = Kv
        if self._verbose:
            print
            print self._name+': setting pressure coefficient to '+`Kv`+' kg/Pa-s'
        self._setParameters(vv)

    def setMaster(self, master):
        """Set the master flow controller."""
        _cantera.flowdev_setMaster(self.flowdev_id(),
                                   master.flowdev_id())
            
    def set(self, Kv = -1.0, master = None):
        if master:
            self.setMaster(master)
        if Kv > 0.0:
            self.setPressureCoeff(Kv)
            

            
#------------- Wall ---------------------------

_wallcount = 0

class Wall:
    """
    Reactor walls.
    
    A Wall separates two reactors, or a reactor and a reservoir. A
    wall has a finite area, may conduct or radiate heat between the
    two reactors on either side, and may move like a piston.

    Walls are stateless objects in Cantera, meaning that no
    differential equation is integrated to determine any wall
    property. Since it is the wall (piston) velocity that enters the
    energy equation, this means that it is the velocity, not the
    acceleration or displacement, that is specified. The wall
    velocity is computed from
    \f[
    v = K(P_{\\rm left} - P_{\\rm right}) + v_0(t),
    \f]
    where $K$ is a non-negative constant, and \f$v_0(t)$ is a
    specified function of time. The velocity is positive if the wall is
    moving to the right.

    The heat flux through the wall is computed from
    \f[
    q = U(T_{\\rm left} - T_{\\rm right}) + \epsilon\sigma (T_{\\rm left}^4
    - T_{\\rm right}^4) + q_0(t),
    \f]
    where \f$ U \f$ is the overall heat transfer coefficient for
    conduction/convection, and \f$ \\epsilon \f$ is the emissivity.
    The function \f$ q_0(t)$ is a specified function of time.
    The heat flux is positive when heat flows from the reactor on the left
    to the reactor on the right.

    A heterogeneous reaction mechanism may be specified for one or
    both of the wall surfaces. The mechanism object (typically an
    instance of class Interface) must be constructed so that it is
    properly linked to the object representing the fluid in the
    reactor the surface in question faces. The surface temperature on
    each side is taken to be equal to the temperature of the reactor
    it faces.
    
    """
    def __init__(self, left, right, name = '',
                 A = 1.0, K = 0.0, U = 0.0,
                 Q = None, velocity = None,
                 kinetics = [None, None]):
        """
        Constructor arguments:
        
        left - Reactor or reservoir on the left. Required.

        right - Reactor or reservoir on the right. Required.

        name - Name string.
        If omitted, the name is 'Wall_n', where 'n' is an integer
        assigned in the order walls are created.

        A - Wall area [m^2]. Defaults to 1.0 m^2.

        K - Wall expansion rate parameter [m/s/Pa]. Defaults to 0.0.

        U - Overall heat transfer coefficient [W/m^2]. Defaults to 0.0
        (adiabbatic wall).

        Q - Heat flux function \f$ q_0(t) \f$ [W/m^2]. Optional. Default:
        \f$ q_0(t) = 0.0 \f$.

        velocity - Wall velocity function \f$ v_0(t) \f$ [m/s].
        Default: \f$ v_0(t) = 0.0 \f$.

        kinetics - Surface reaction mechanisms for the left-facing and
        right-facing surface, respectively. These must be instances of
        class Kinetics, or of a class derived from Kinetics, such as
        Interface. If chemistry occurs on only one side, enter 'None'
        for the non-reactive side.

        """
        typ = 0
        self.__wall_id = _cantera.wall_new(typ)

        global _wallcount
        if name == '':
            _nm = 'Wall_'+`_wallcount`
        else:
            _nm = name
        _wallcount += 1
        
        if left and right:
            self.install(left, right)
        else:
            raise CanteraError('both left and right reactors must be specified.')
        self.setArea(A)
        self.setExpansionRateCoeff(K)
        self.setVelocity(velocity)        
        self.setHeatTransferCoeff(U)
        self.setHeatFlux(Q)

        self.setKinetics(kinetics[0],kinetics[1])

    def __del__(self):
        """ Delete the Wall instance. This method is called
        automatically when no Python object stores a reference to this
        Wall. Since reactors and reserviors store references to all
        Walls installed on them, this method will only be called after
        the reactors/reservoirs have been deleted.  """
        
        _cantera.wall_del(self.__wall_id)
        
    def ready(self):
        """
        Return 1 if the wall instance is ready for use, 0 otherwise. Deprecated.
        """
        return _cantera.wall_ready(self.__wall_id)

    def area(self):
        """
        The wall area (m^2).
        """
        return _cantera.wall_area(self.__wall_id)

    def setArea(self, a):
        """
        Set the area (m^2). The wall area may be changed manually at any time during a simulation.
        """
        _cantera.wall_setArea(self.__wall_id, a)

    def setThermalResistance(self, rth):
        """Deprecated."""
        return _cantera.wall_setThermalResistance(self.__wall_id, rth)

    def setHeatTransferCoeff(self, u):
        """
        Set the overall heat transfer coefficient [W/m^2/K]
        """
        return _cantera.wall_setHeatTransferCoeff(self.__wall_id, u)

    def setEmissivity(self, epsilon):
        """
        Set the emissivity.
        """
        _cantera.wall_setEmissivity(self.__wall_id, epsilon)
        
    
    def setHeatFlux(self, qfunc):
        """
        Specify the time-dependent heat flux function [W/m2].
        'qfunc' must be a functor (an instance of a subclass of Cantera.Func1).
        See: Func1.
        """
        n = 0
        if qfunc: n = qfunc.func_id()
        return _cantera.wall_setHeatFlux(self.__wall_id, n)

    def setExpansionRateCoeff(self, k):
        """Set the coefficient K that determines the expansion rate
        resulting from a unit pressure drop."""
        _cantera.wall_setExpansionRateCoeff(self.__wall_id, k)        
        
    def setVelocity(self, vfunc):
        """
        Specify the velocity function [m/s]. 'vfunc' must
        be a functor (an instance of a subclass of Cantera.Func1)
        See: Func1.
        """
        n = 0
        if vfunc: n = vfunc.func_id()
        _cantera.wall_setVelocity(self.__wall_id, n)

    def vdot(self):
        """Rate of volume change [m^3]. A positive value corresponds
        to the left-hand reactor volume increasing, and the right-hand
        reactor volume decreasing."""
        return _cantera.wall_vdot(self.__wall_id)

    def velocity(self):
        return self.vdot()/self.area()

    def heatFlowRate(self):
        """Rate of heat flow through the wall. A positive value
        corresponds to heat flowing from the left-hand reactor to the
        right-hand one."""
        return _cantera.wall_Q(self.__wall_id)

    def heatFlux(self):
        return self.heatFlowRate()/self.area()
    
    def install(self, left, right):
        left._addWall(self, right)
        right._addWall(self, left)
        _cantera.wall_install(self.__wall_id, left.reactor_id(),
                               right.reactor_id())

    def setKinetics(self, left, right):
        """Specify surface reaction mechanisms for the left and right sides of the wall."""
        ileft = 0
        iright = 0
        if left:
            ileft = left.kin_index()
        if right:
            iright = right.kin_index()
        _cantera.wall_setkinetics(self.__wall_id, ileft, iright)
                                  
    def set(self, **p):
        """Set various wall parameters: 'A', 'U', 'K', 'Q'. 'velocity'.
        These have the same meanings as in the constructor.
        """
        for item in p.keys():
            if item == 'A' or item == 'area':
                self.setArea(p[item])
            elif item == 'R':
                self.setThermalResistance(p[item])
            elif item == 'U':
                self.setHeatTransferCoeff(p[item])                
            elif item == 'K':
                self.setExpansionRateCoeff(p[item])
            elif item == 'Q':
                self.setHeatFlux(p[item])
            elif item == 'velocity':
                self.setVelocity(p[item])
            else:
                raise 'unknown parameter: ',item
                

class ReactorNet:
    
    """Networks of reactors. ReactorNet objects are used to
    simultaneously advance the state of a set of coupled reactors.

    Example:

    >>> r1 = Reactor(gas1)
    >>> r2 = Reactor(gas2)
    >>> <... install walls, inlets, outlets, etc...>

    >>> reactor_network = ReactorNet([r1, r2])
    >>> reactor_network.advance(time)
    
    """


    def __init__(self, reactorlist = None):
        """
        Create a new ReactorNet instance. If a list of reactors is supplied,
        these will be added to the network.
        """
        self._reactors = []
        self.__reactornet_id = _cantera.reactornet_new()
        if reactorlist:
            for r in reactorlist:
                self.add(r)


    def __del__(self):
        """Delete the reactor network instance. The reactors in the
        network are not deleted."""
        _cantera.reactornet_del(self.__reactornet_id)


    def reactornet_id(self):
        """ The integer index used to access the
        kernel reactornet object. For internal use.  """
        return self.__reactornet_id

    
    def add(self, reactor):
        """
        Add a reactor to the network.
        """
        self._reactors.append(reactor)
        _cantera.reactornet_addreactor(self.__reactornet_id,
                                       reactor.reactor_id())

        
    def setInitialTime(self, t0):
        """Set the initial time. Restarts integration from this time
        using the current state as the initial condition. Default: 0.0 s"""
        _cantera.reactornet_setInitialTime(self.__reactornet_id, t0)

    def advance(self, time):
        """Advance the state of the reactor network in time from the current
        time to time 'time'."""
        return _cantera.reactornet_advance(self.__reactornet_id, time)

    def step(self, time):
        """Take a single internal time step toward time 'time'.
        The time after taking the step is returned."""
        return _cantera.reactornet_step(self.__reactornet_id, time)    

    def reactors(self):
        return self._reactors
    
