import json
import unittest
import test.mdoodz


class Test(unittest.TestCase):
    model_name = 'ShearTemplate'
    model: test.mdoodz.model.MdoodzModel

    @classmethod
    def setUpClass(cls):
        """Compiling C code"""
        cls.model = test.mdoodz.model.MdoodzModel(cls.model_name)
        cls.model.compile()

    def test_iteration_count(self):
        """Iteration number count should be 2 if pwlv = 0"""
        self.model.run()

        results_hdf5 = self.model.get_results()
        number_of_steps = results_hdf5.get('Iterations/NumberSteps')[0]
        self.assertEqual(number_of_steps, 2)


if __name__ == '__main__':
    unittest.main()
