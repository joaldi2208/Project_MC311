import unittest
import importlib
from Hashed_Keys import HashedKey
        
class Hashing(unittest.TestCase):
    def test_uniqueSubstrings(self):
        # test algorithm
        smiles1 = HashedKey("ABC")
        uniqueSubstring_1 = smiles1.createUniqueSubstrings()
        smiles2 = HashedKey("AABC")
        uniqueSubstring_2 = smiles2.createUniqueSubstrings()
        self.assertEqual(uniqueSubstring_1, ["A","AB","ABC","B","BC","C"])
        self.assertEqual(uniqueSubstring_2, ["A","AA","AAB","AABC","AB","ABC","B","BC","C"])
