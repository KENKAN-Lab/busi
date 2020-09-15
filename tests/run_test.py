import unittest


def get_suite():

    test_suite = unittest.TestSuite()

    all_test_suite = unittest.defaultTestLoader.discover("tests", pattern="test_*.py")

    for ts in all_test_suite:
        test_suite.addTest(ts)

    return test_suite


if __name__ == "__main__":

    suite = get_suite()
    unittest.TextTestRunner().run(suite)
