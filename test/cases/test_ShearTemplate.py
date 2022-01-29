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
        """Should be converged with one step if pwlv = 0"""
        self.model.run()

        results_hdf5 = self.model.get_results()
        number_of_steps = results_hdf5.get('Iterations/NumberSteps')[0]
        self.assertEqual(number_of_steps, 1)


if __name__ == '__main__':
    unittest.main()
