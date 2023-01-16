from openbabel import openbabel
from collections import Counter
import re 

def smiles2ConnectionTable(smiles):
    """convert a SMILES string to a connection table."""

    # maybe in a class
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "ct")
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, smiles)
    outMDL = obConversion.WriteString(mol)
    

    MDL = outMDL.strip()
    MDL_first_line = MDL[:8]
    n_rows_connection_table = re.search(r"[0-9]+ ([0-9]*)", MDL_first_line).group(1)

    lines_MDL_file = MDL.split("\n")

    lines = []
    for line in lines_MDL_file:
        line = line.strip()
        line = re.split(r'\s+', line)
        lines.append(line)

    connection_table = lines[-int(n_rows_connection_table):]
    return connection_table





def getBondPairs(connection_table):
    bond_pairs = []
    for row in connection_table:
        bond_pair = (row[0], row[1])
        bond_pairs.append(bond_pair)
    return bond_pairs



def getBondNumbers(connection_table):
    """get the number of bonds for each atom"""

    cnt = Counter()
    for row in connection_table:

        cnt[row[0]] += 1
        cnt[row[1]] += 1

    n_unique_numbers = len(set(cnt.values()))
    return cnt, n_unique_numbers
                           

def transferNumber(bond_pairs, bond_numbers, n_unique_numbers1):
    """transfers the bond number until it starts to fluctuate"""

    updated_bond_numbers = dict.fromkeys(list(bond_numbers),0)
    for atom1, atom2 in bond_pairs:
        updated_bond_numbers[atom1] += bond_numbers[atom2]
        updated_bond_numbers[atom2] += bond_numbers[atom1]

    n_unique_numbers = len(set(updated_bond_numbers.values()))
    #print(updated_bond_numbers)
    print(n_unique_numbers)

    if n_unique_numbers > n_unique_numbers1:
        transferNumber(bond_pairs, updated_bond_numbers, n_unique_numbers)
    else:
        return updated_bond_numbers
                


if __name__=='__main__':

    #connection_table = smiles2ConnectionTable("O[C@@H](Cc1ccccc1)C(O)=O")
    connection_table = smiles2ConnectionTable("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
    print(len("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"))
    #print(connection_table)
    bond_pairs = getBondPairs(connection_table)
    bond_numbers, n_unique_numbers1 = getBondNumbers(connection_table)
    #print(bond_pairs)
    #print(bond_numbers["2"])
    updated_bond_numbers = transferNumber(bond_pairs, bond_numbers, n_unique_numbers1)
    # calc morgan symmetrie

