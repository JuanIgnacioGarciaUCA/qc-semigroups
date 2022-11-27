import QArithmetic
import math
from qiskit import QuantumRegister, QuantumCircuit

mult1 = QuantumRegister(2, 'mult1')
mult2 = QuantumRegister(2, 'mult2')
result = QuantumRegister(4, 'result')

qc1 = QuantumCircuit(mult1, mult2, result)
QArithmetic.mult(qc1, mult1, mult2, result, 2)
# print(qc1.inverse())

qc2 = QuantumCircuit(mult1, mult2, result)
qc2.append(qc1, [*mult1, *mult2, *result])
qc2.append(qc1.inverse(), [*mult1, *mult2, *result])

# print(qc2)

print([x for x in mult1])
xx = [1, 2, 3]
xx.append(2)
print(xx)
