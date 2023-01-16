from open_babel import smiles2ConnectionTable
from open_babel import getBondOrder
from open_babel import getBondPairs




def identify_last_character (table):
    lines = table.split("\n")
    character_list = []
    for line in lines:
        character_list.append(line[-1])
    return character_list

def find_atoms (character_list):
    numbers = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
    atoms = []
    for character in character_list:
        if character not in numbers:
            atoms.append(character)
    return atoms

def get_index_hetero_atoms (atoms_list):
    hetero_atoms = []
    for index, atom in enumerate(atoms_list):
        if atom != "C":
            hetero_atoms.append(index+1)
    return hetero_atoms

def check_bond_type (bond_list):
    bonds = []
    for index, bond in enumerate(bond_list):
        if bond != "1":
            bonds.append(index)
    return bonds

def get_binding_partners (bonds, bond_pairs):
    partners = []
    for bond in bonds:
        partners.append(bond_pairs[bond][0])
        partners.append(bond_pairs[bond][1])
    return partners

def find_acetal (bondpairs, heteroatom_index):
    possible_acetals = []
    for pairs in bondpairs:
        for index in  heteroatom_index:
            if str(index) in pairs:
                possible_acetals.append(pairs[0])
                possible_acetals.append(pairs[1])
    return possible_acetals





if __name__=='__main__':

    connection_table, MDL = smiles2ConnectionTable("O[C@@H](Cc1ccccc1)C(O)=O")

    #print(MDL)
    #find_atoms(identify_last_character(MDL))

    BondPairs = getBondPairs(connection_table)
    BondOrder = getBondOrder(connection_table)

    
    print(BondOrder)
    print(BondPairs)

    higher_bonds = check_bond_type(BondOrder)
    partners_higher_bonds = get_binding_partners(higher_bonds, BondPairs)
    last_character = identify_last_character (MDL)
    atoms_list = find_atoms(last_character)
    index_heteroatoms = get_index_hetero_atoms (atoms_list)


    possible_acetals = find_acetal(BondPairs, index_heteroatoms)

    print(possible_acetals)
    print(index_heteroatoms)