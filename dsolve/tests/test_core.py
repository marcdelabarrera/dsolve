from atoms import Variable, Parameter, E
import unittest

class TestVariable(unittest.TestCase):
    def test_initialize(self):
        self.assertEqual(str(Variable('x_{t}')), 'x_{t}')
        self.assertEqual(str(Variable('E_{t}[x_{t+1}]')), 'E_{t}[x_{t+1}]')
        self.assertEqual(str(Variable('E_{t}[x_{i,t+1}]')), 'E_{t}[x_{i,t+1}]')

    def test_expectation(self):
        self.assertEqual(str(E(Variable('x_{t+1}'),'t')),'E_{t}[x_{t+1}]')
        self.assertEqual(str(E(Variable('x^{p}_{t+1}'),'t-1')),'E_{t-1}[x^{p}_{t+1}]')

if __name__ == '__main__':
    unittest.main()