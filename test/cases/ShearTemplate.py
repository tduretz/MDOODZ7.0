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

    def test_string(self):
        """FIRST TEST"""
        a = 'some'
        b = 'some'
        self.assertEqual(a, b)

    def test_boolean(self):
        """SECOND TEST"""
        a = True
        b = False
        self.assertEqual(a, b)


if __name__ == '__main__':
    unittest.main()
