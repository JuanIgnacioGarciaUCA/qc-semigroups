import QArithmetic
import math
from qiskit import QuantumRegister, QuantumCircuit

mult1 = QuantumRegister(3, 'mult1')
mult2 = QuantumRegister(3, 'mult2')
result = QuantumRegister(6, 'result')

qc1 = QuantumCircuit(mult1, mult2, result)
QArithmetic.mult(qc1, mult1, mult2, result, 3)
print(qc1)

qc2 = QuantumCircuit(mult1, mult2, result)
qc2.append(qc1, [*mult1, *mult2, *result])
# print(qc2)
