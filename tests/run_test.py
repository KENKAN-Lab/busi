import unittest


def get_suite():

    return unittest.defaultTestLoader.discover("tests", pattern="test_*.py")


if __name__ == "__main__":

    suite = get_suite()
    unittest.TextTestRunner().run(suite)
