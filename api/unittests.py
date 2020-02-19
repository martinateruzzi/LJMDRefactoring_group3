import unittest
from ljmd import *

class Test(unittest.TestCase):
	md = Ljmd("argon_003.inp")

	def test_force(self):
		