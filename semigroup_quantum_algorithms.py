import QArithmetic
import math
from qiskit import QuantumRegister, QuantumCircuit


class groverPhase():

    def __init__(self, initH=False):
        self.initH = initH  # Ensamblado con el H inicial?

    def assembleSubCirc(func):

        def inner(self, qCircuit, *args, **kwargs):
            registerClass = type(QuantumRegister(1, 'dummy'))
            dummyCircuitSize = 0
            wireList = []
            wires = []
            for var in args:
                if type(var) == type(list()):
                    for subVar in var:
                        wires.append(subVar)
                elif type(var) == registerClass:
                    wires.append(var)

            for wire in wires:
                dummyCircuitSize += len(wire)
                wireList += list(wire)
            dummyCircuit = QuantumCircuit(*wires, name=f'{func.__name__}')
            if 'dagger' in kwargs.keys():
                dagger = kwargs['dagger']
            else:
                dagger = False

            if dagger:
                qCircuit.append(func(self, dummyCircuit, *args).inverse(), wireList)
            else:
                qCircuit.append(func(self, dummyCircuit, *args), wireList)

        return inner

    @assembleSubCirc
    def linearCombination(self, qCircuit, wiresOfGen, numberOfGen, wiresOfLambda,
                          wiresOfLinear, dagger=False):
        '''operador de combinación lineal'''

        dummy = QuantumCircuit(*wiresOfGen, *wiresOfLambda, *wiresOfLinear,
                               name='linearCombination')

        num = len(wiresOfGen[0])
        for i in range(numberOfGen):

            QArithmetic.mult(dummy, wiresOfGen[i], wiresOfLambda[i],
                             wiresOfLinear[i][:2*num], num)

        for i in range(numberOfGen-1):

            num = len(wiresOfLinear[i])
            QArithmetic.add(dummy, wiresOfLinear[i], wiresOfLinear[i+1], num)

        return dummy

    @assembleSubCirc
    def substractSought(self, qCircuit, wiresOfSought, wiresOfLinearCom, dagger=False):

        QArithmetic.sub(qCircuit, wiresOfSought, wiresOfLinearCom, len(wiresOfSought)-1)
        return qCircuit

    @assembleSubCirc
    def isZero(self, qCircuit, LinCtrl, target, *tCtrl):

        qCircuit.x(LinCtrl)
        qCircuit.mct([*LinCtrl, *tCtrl], target)

        return qCircuit

    @assembleSubCirc
    def difussor(self, qCircuit, *lambdaWires):
        '''Crea un difusor en los hilos implicados'''
        nlw = len(lambdaWires)
        qCircuit.h(range(nlw))
        qCircuit.x(range(nlw))

        qCircuit.h(0)
        qCircuit.mct(list(range(1, nlw)), 0)
        qCircuit.h(0)

        qCircuit.x(range(nlw))
        qCircuit.h(range(nlw))

        return qCircuit

    def induceSuperposition(self, qCircuit, wires):

        if not self.initH:

            for i in range(len(wires)):
                qCircuit.h(wires[i])
            self.initH = True

    def semigroupMembershipOracle(self, qCircuit, tControl,  wiresOfGenerators,
                                  numberOfGenerators, wiresOfLambda, wiresOfLinCom,
                                  wiresOfSought, threadPhaseKickback):

        self.linearCombination(qCircuit, wiresOfGenerators, numberOfGenerators,
                               wiresOfLambda, wiresOfLinCom)
        self.substractSought(qCircuit, wiresOfSought, wiresOfLinCom[-1])

        self.isZero(qCircuit,  wiresOfLinCom[-1], threadPhaseKickback, tControl)

        self.substractSought(qCircuit, wiresOfSought, wiresOfLinCom[-1], dagger=True)
        self.linearCombination(qCircuit, wiresOfGenerators, numberOfGenerators,
                               wiresOfLambda, wiresOfLinCom, dagger=True)


class semigroup(groverPhase):
    '''El numero opt. de iter (sol simple) es sqrt(2^n)'''

    def __init__(self, sizeOfLambda=None, *generators):
        super().__init__()
        self.aproximationQBits = 1
        self.soughtElement = 0
        self.generators = generators
        self.numberOfGenerators = len(generators)
        if sizeOfLambda is None:
            self.sizeOfLambda = max([math.ceil(math.log2(gen)) for gen in generators])
        else:
            self.sizeOfLambda = sizeOfLambda

    def dec2binQR(self, qReg, number):  # ****

        '''Toma un valor decimal y un registro, introduce el equiv. binairo en el
        registro.'''

        binExp = bin(number)  # −> 0bxxxxx...

        if len(qReg) < len(binExp)-2:
            print(f"error al convertir en binario{len(qReg)} and {len(binExp)-2}")
            return
        i = 0
        # De derecha a izquierda
        while binExp[-(i+1)] != 'b':
            if binExp[-(i+1)] == '1':
                self.circ.x(qReg[i])
            i += 1
        del(i, binExp)

    def setUpWires(self):

        # Tienen que coincidir para poder multiplicarlos.
        sizeOfGenerators = self.sizeOfLambda

        if self.aproximationQBits != 0:
            self.wiresQCounting = []
            for i in range(self.aproximationQBits):

                self.wiresQCounting.append(QuantumRegister(1, f't_{i}'))

        self.wiresOfGenerators = [QuantumRegister(sizeOfGenerators, f's{i}')
                                  for i in range(self.numberOfGenerators)]

        self.wiresOfLambda = [QuantumRegister(self.sizeOfLambda, f'lambda{i}')
                              for i in range(self.numberOfGenerators)]

        self.wiresOfLinCom = []

        for j in range(self.numberOfGenerators-1):

            self.wiresOfLinCom.append(QuantumRegister(self.sizeOfLambda*2 + j,
                                                      f'linearCom{j}'))

        self.wiresOfLinCom.append(
            QuantumRegister(self.sizeOfLambda*2 + self.numberOfGenerators,
                            f'linearCom{self.numberOfGenerators-1}'))

        # Same size to substract
        self.wiresOfSought = QuantumRegister(self.sizeOfLambda*2 +
                                             self.numberOfGenerators, 'a')

        self.threadPhaseKickback = QuantumRegister(1, 'b')

        self.circ = QuantumCircuit(
            *self.wiresQCounting, *self.wiresOfGenerators, *self.wiresOfLambda,
            *self.wiresOfLinCom, self.wiresOfSought, self.threadPhaseKickback,
            name='Semigroup Membership Algorithm')

    def setUpPhaseKickback(self):

        self.circ.x(self.threadPhaseKickback)
        self.circ.h(self.threadPhaseKickback)

    def setUpValues(self):

        for i in range(self.numberOfGenerators):
            self.dec2binQR(self.wiresOfGenerators[i], self.generators[i])
        self.dec2binQR(self.wiresOfSought, self.soughtElement)

    def quantumCountingInit(self):
        # Esta parte es común al inicio de todos los algoritmos
        self.setUpWires()
        self.setUpValues()
        self.setUpPhaseKickback()
        for i in range(self.aproximationQBits):

            self.circ.h(self.wiresQCounting[i])
        self.circ.barrier()

        self.induceSuperposition(self.circ, self.wiresOfLambda)

        self.circ.barrier()

    def quantumCountingParam(self):
        # definimos el elemento que vamos a buscar
        try:
            self.soughtElement = int(input('Introduce elemento a comprobar: '))
            self.aproximationQBits = int(input('Bits de aproximación: '))
            if self.aproximationQBits <= 0:
                raise ValueError
        except ValueError:
            raise ValueError

    def semigroupMembershipAlgorithm(self):

        try:
            self.quantumCountingParam()
            self.quantumCountingInit()
            for i in range(self.aproximationQBits):
                for j in range(i+1):
                    self.circ.barrier()
                    self.semigroupMembershipOracle(
                        self.circ, self.wiresQCounting[i], self.wiresOfGenerators,
                        self.numberOfGenerators, self.wiresOfLambda,
                        self.wiresOfLinCom[:self.numberOfGenerators],
                        self.wiresOfSought,
                        self.threadPhaseKickback)

                    self.circ.barrier()  # Si no pones las barreras las cajas se mezclan

                    self.difussor(self.circ, *self.wiresOfLambda)

                    self.circ.barrier()

            if self.aproximationQBits != 0:
                pass
              #  self.circ.h(self.wiresQCounting)  # QFT == H para qubits = 1
        except ValueError:
            print('Algo salió mal.')

test1 = semigroup(None, 7, 3)
test1.semigroupMembershipAlgorithm()
print(test1.circ)